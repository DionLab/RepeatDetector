#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include "pb_hdf5.h"
#include "pfProfile.h"
#include "pfCompute.h"
#include "pfOutput.h"
#include "pfSequence.h"
#include "pfStatistics.h"
#include <smmintrin.h>

char OutputFileName[256] = "Results";
struct IO_Wrapper OptionsDefault = { .OptimalOnly = false };
struct IO_Data OptionsData = {
	.BaseFileName = OutputFileName,
	.IsBinary = false,
	.Separate = false
};
#ifdef USE_GD
struct IO_PNG  OptionsPNG  = {
	.BaseFileName = OutputFileName,
	.PixelWidth = 1,
	.PixelHeight = 1,
	.ScaleMin = 0,
	.ScaleMax = 0,
	.x = { 0, 0},
	.y = { 0, 0},
	.GivenRegion = false,
	.GivenScale = false,
	.ScaleOnImage = false,
	.WithReverse = false,
	.compress = {1,1}
};
#endif
#ifdef USE_PDF
struct IO_PDF OptionsPDF  = { 
	.BaseFileName = OutputFileName,
	.WithPaths = true, 
	.WithProfileScoreGraph = false,
	.WithSequenceScoreGraph = false,
	.WithReverse = false,
	.GivenRegion = false,
	.x = {0,0},
	.y = {0,0}
};
#endif

OutputMethod Output = DataOutput;
_Bool PacBioAveragedPerHole = false;
const void * OutputMethodOptions = &OptionsDefault;
_Bool OutputVerbose = false;
const _Bool CompleteCycleoOnly = false;
int PlotAllSequence = 0;
unsigned int OutputPrintWidth = 100;

int align(const PacBio_t * const restrict PB, const struct Profile * const restrict prf,
					struct Alignment * * Alignments, const unsigned int HoleNumber)
{
	if (!(PB->Content & BASECALLING)) return ERROR; 
	HoleData_t HReg = HOLE_DATA_INIT;
	int res = getHoleRegions(PB, &HReg, HoleNumber);
	if (res < 0) return res;
	
	/* Get the HQRegion from HReg starting for the end to speed things up */
	unsigned int count = HReg.nRegions;
	int HQStart, HQStop;
	int LargestSequence = 0;
	{
		int HQindex = 0;
		_Bool HasValidHQ = true;
		HoleRegion_t * const HR = HReg.Regions;
		/* WARNING: That assumes inserts are always prior to HQregion !! */
		for (int r=0; r<count; r++) {
			//fprintf(stderr, "%i : %u - %u\n", HR[r].type, HR[r].start, HR[r].stop);
			if (HR[r].type == Insert) {
				const int lstop = HR[r].stop;
				LargestSequence = (LargestSequence > lstop) ? LargestSequence : lstop;
			}
			else if (HR[r].type == HQRegion) {
				if (HR[r].stop == 0) {
					HasValidHQ = false;
				}
				HQindex = r;
			}
		}
		
		HQStart = HR[HQindex].start;
		HQStop  = HR[HQindex].stop;
		
		if (HR[HQindex].quality < PacBioHQThreshold) {
			HQStop = 0;
		}
		
		if (PacBioKeepInvalid && !HasValidHQ) {
			HQStart = 0;
			HQStop  = LargestSequence;
		}
	}
	if (OutputVerbose) fprintf(stderr, "Hole %u HQ region set to [%i,%i]\n", HoleNumber, HQStart, HQStop);
	

	HReg.Regions[0].type  = Insert;
	HReg.Regions[0].start = 0 /*HQStart*/;
	HReg.Regions[0].stop  = HQStop;
	HReg.Regions[1].type  = HQRegion;
	HReg.Regions[1].start = 0 /*HQStart*/;
	HReg.Regions[1].stop  = HQStop;
	if (OutputVerbose) fprintf(stderr, "alignHoleReads treating hole %u as a stream of length %i\n",
	                           HoleNumber, HQStop);
	count = HReg.nRegions = 2;
	LargestSequence++;
	
	HoleReads_t HReads = HOLE_READS_INIT;
	res = getHoleReads_HR(PB, &HReads, &HReg);
	if (res < 0) return res;
	
	unsigned char RevComp_Alphabet_Mapping[ALPHABET_SIZE+2] __attribute__((aligned(16)));

	const size_t prfLength = prf->Length;
	const size_t WorkSize  = prf->Length + 1;
	
	if (OutputVerbose) fprintf(stderr, "alignHoleReads on hole %u with cutoff set to %i\n", HoleNumber, prf->CutOff);
		
	/*************************************************************************/
  /*                          ALLOCATE MEMORY                              */
  /*************************************************************************/
	
	union lScores * restrict matrix   = _mm_malloc(WorkSize*LargestSequence*sizeof(union lScores), 64);
	union lScores * restrict rvmatrix = _mm_malloc(WorkSize*LargestSequence*sizeof(union lScores), 64);
	int * restrict WORK     = _mm_malloc(2*WorkSize*4*sizeof(int)+63,64);
	unsigned char * restrict SequenceIndex = (unsigned char*) malloc(LargestSequence*sizeof(unsigned char));
	if ( rvmatrix == NULL || matrix == NULL || WORK == NULL || SequenceIndex == NULL ) {
		res = -1; 
		goto FIN;
	}
	
	PFSequence PFSeq;
	PFSeq.ProfileIndex = SequenceIndex;

	/* Builds RevComp index */
	memcpy(RevComp_Alphabet_Mapping, prf->Alphabet_Mapping, (ALPHABET_SIZE+2)*sizeof(unsigned char));
	RevComp_Alphabet_Mapping[(int) 'A' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'T' - (int) 'A' ];
	RevComp_Alphabet_Mapping[(int) 'C' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'G' - (int) 'A' ];
	RevComp_Alphabet_Mapping[(int) 'G' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'C' - (int) 'A' ];
	RevComp_Alphabet_Mapping[(int) 'T' - (int) 'A'] =  prf->Alphabet_Mapping[(int) 'A' - (int) 'A' ];
	
	const size_t SeqLength = HReg.Regions[0].stop - HReg.Regions[0].start;
	char * const restrict seq =  &(HReads.Stream[HReg.Regions[0].start]);
	
	/* Copy sequence to local space index */ 
	memcpy(PFSeq.ProfileIndex, seq, SeqLength);
	PFSeq.Length = SeqLength;
	
	/* Translate into indices */
	TranslateSequenceToIndex(&PFSeq, prf->Alphabet_Mapping);
	
	/* Build matrix */
	Repeat_sse41.BuildMatrix(prf, PFSeq.ProfileIndex, matrix, WORK, NULL, 0, SeqLength);
	
	struct Alignment * als;
	const int c1 = Repeat_sse41.GetAlignments(matrix, prf, &als, SeqLength);
	if (c1 < 0) {res = -1; goto FIN2;}
	
	/* Copy sequence to local space index */ 
	memcpy(PFSeq.ProfileIndex, seq, SeqLength);
	
	/* Translate into indices */
	TranslateSequenceToIndex(&PFSeq, RevComp_Alphabet_Mapping);
	ReverseTranslatedSequence(&PFSeq);
	
	/* Build reverse matrix */
	Repeat_sse41.BuildMatrix(prf, PFSeq.ProfileIndex, rvmatrix, WORK, NULL, 0, SeqLength);
	struct Alignment * rvals;
	const int c2 = Repeat_sse41.GetAlignments(rvmatrix, prf, &rvals, SeqLength);
	if (c2 < 0) {res = -1; goto FIN3;}
	
	*Alignments = malloc((c1+c2)*sizeof(struct Alignment));
	if (*Alignments) {
		for (int i=0; i<c1; i++) (*Alignments)[i] = als[i];
		for (int i=0; i<c2; i++) {
			(*Alignments)[c1+i] = rvals[i];
			(*Alignments)[c1+i].Matrix.row.Begin = (int) SeqLength - rvals[i].Matrix.row.End;
			(*Alignments)[c1+i].Matrix.row.End   = (int) SeqLength - rvals[i].Matrix.row.Begin;
		}
		res = c1+c2;
	}
	else {
		res = -1;
	}
	
	FIN4:
		if (c1) free(rvals);
	FIN3:
		if (c2) free(als);
	FIN2:
	if (WORK) _mm_free(WORK);
	if (matrix) _mm_free(matrix);
	if (rvmatrix) _mm_free(rvmatrix);
	free(SequenceIndex);
	FIN:
	
	return res;
}

int main(int argc, char *argv[])
{
	if (argc > 4) {
		fprintf(stderr, "Usage: %s [pac bio HDF5 file] [hole number] [Window Mean Size]\n", argv[0]);
		exit(1);
	}
	
	/* Loading PacBio Trace File */
	
	fprintf(stderr, "Testing trace file %s\n", argv[1]);
	PacBio_t * const restrict PBTrc = PacBioOpen(argv[1]);
	if (PBTrc == NULL) {
		fprintf(stderr, "Error opening Pac Bio file %s\n", argv[1]);
		exit(1);
	}
	
	fprintf(stderr,"Base file name : %s\n"
				 "Directory      : %s\n"
				 "Parts          : %u\n",
				PBTrc->BaseFileName, PBTrc->Directory, PBTrc->nParts);
	for (unsigned int i=0; i<PBTrc->nParts; i++) {
		fprintf(stderr,"               : %s\n", PBTrc->PartsFileName[i]);
	}
	
	fputs("The given file [set] contains : ", stderr);
	if (PBTrc->Content & CONSENSUS_BASECALLING) fputs("CONSENSUS ", stderr);
	if (PBTrc->Content & BASECALLING) fputs("BASECALLING ", stderr);
	if (PBTrc->Content & PULSE) fputs("PULSES ", stderr);
	if (PBTrc->Content & TRACE) fputs("TRACES ", stderr);
	fputs("\n", stderr);
	
	fprintf(stderr,"Holes          : %u\n", PBTrc->nHoles);
	
	/* Loading PacBio Bas file */
	char BASFileName[256];
	snprintf(BASFileName, 256, "%s/Analysis_Results/%.*sbas.h5", 
					 PBTrc->Directory,
					 (int) (strlen(PBTrc->BaseFileName) - 6),
					 PBTrc->BaseFileName);
	fprintf(stderr,"Looking for %s...\t", BASFileName);
	PacBio_t * const restrict PBBas = PacBioOpen(BASFileName);
	if (PBBas == NULL) {
		fputs("FAILED\n", stderr);
		goto END;
	}
	else 
		fputs("OK\n", stderr);
	
	fprintf(stderr, "Testing Basecalling file %s\n", BASFileName);
	fprintf(stderr,"Base file name : %s\n"
				 "Directory      : %s\n"
				 "Parts          : %u\n",
				PBBas->BaseFileName, PBBas->Directory, PBBas->nParts);
	for (unsigned int i=0; i<PBBas->nParts; i++) {
		fprintf(stderr,"               : %s\n", PBBas->PartsFileName[i]);
	}
	
	fputs("The given file [set] contains : ", stderr);
	if (PBBas->Content & CONSENSUS_BASECALLING) fputs("CONSENSUS ", stderr);
	if (PBBas->Content & BASECALLING) fputs("BASECALLING ", stderr);
	if (PBBas->Content & PULSE) fputs("PULSES ", stderr);
	if (PBBas->Content & TRACE) fputs("TRACES ", stderr);
	fputs("\n", stderr);
	fprintf(stderr,"Holes          : %u\n", PBBas->nHoles);
	
	/* Checking content of both TRc and Bas Pac Bio files */
	if (!(PBTrc->Content & TRACE)) {
		fprintf(stderr, "%s does NOT contains traces...\n", PBTrc->BaseFileName);
		goto END;
	}
	if (!(PBBas->Content & BASECALLING)) {
		fprintf(stderr, "%s does NOT contains basecalling...\n", PBBas->BaseFileName);
		goto END;
	}
	
	if (PopulateHoleStatus(PBTrc) != SUCCESS) {
			fprintf(stderr,"Hole status not populated...\n");
			goto END;
	}
	
	if (PopulateHoleCoordinates(PBTrc) != SUCCESS) {
		fprintf(stderr,"Hole coordinates not populated...\n");
		goto END;
	}
	
	if (PopulateDecodeTable(PBTrc) != SUCCESS) {
		fputs("Error getting Decoding tables and Bias\n", stderr);
		goto END;
	}
	
	if (IndexRegions(PBBas) != SUCCESS) {
		fputs("Error indexing basecalling regions...\n", stderr);
		goto END;
	}
	
	/* Print out decoding tables and offsets */
// 	{
// 		for (unsigned int i=0; i<PBTrc->nParts; i++) {
// 			printf("%s decode table:\n", PBTrc->PartsFileName[i]);
// 			for (int l=0; l<16; l++) {
// 				for (int k=0; k<7; k++) {
// 					printf("%8.2f ", PBTrc->Traces.DecodeTable[i][8*l+k]);
// 				}
// 				printf("%8.2f\n", PBTrc->Traces.DecodeTable[i][8*l+7]);
// 			}
// 			printf("\n%s decode bias: %8.2f\n\n", PBTrc->PartsFileName[i], PBTrc->Traces.Bias[i]);
// 		}
// 	}


	int start, stop;
	if (sscanf(argv[2], "%i:%i", &start, &stop) != 2) {
		const int holeN = atoi(argv[2]);
		
		fprintf(stderr,"Hole %i is at location [%i,%i]\n", holeN, PBTrc->Coordinates[holeN][0], PBTrc->Coordinates[holeN][1]);
		fputs("This hole is of type ", stderr);
		switch(PBTrc->Status[holeN]) {
			case SEQUENCING: fputs("SEQUENCING\n", stderr); break;
			case ANTIHOLE: fputs("ANTIHOLE\n", stderr); break;
			case FIDUCIAL: fputs("FIDUCIAL\n", stderr); break;
			case SUSPECT: fputs("SUSPECT\n", stderr); break;
			case ANTIMIRROR: fputs("ANTIMIRROR\n", stderr); break;
			case FDZMW: fputs("FDZMW\n", stderr); break;
			case FBZMW: fputs("FBZMW\n", stderr); break;
			case ANTIBEAMLET: fputs("ANTIBEAMLET\n", stderr); break;
			case OUTSIDEFOV: fputs("OUTSIDEFOV\n", stderr); break;
		}

		Traces_t * const Indices = getTraceIndex(PBTrc, holeN);
		const size_t nTimes = PBTrc->Traces.Length;
		
		unsigned int Histogram[256];
		
		computeHistogram(Indices->T, Histogram, nTimes);
		
		/* Check the mean */
		{
			float Mean = 0.0f;
			const float invN = 1.0f/(float) nTimes;
			for (size_t i=0; i<256;i++) {
				Mean += (float) Histogram[i] * invN * PBTrc->Traces.DecodeTable[0][i];
			}
			Mean -= PBTrc->Traces.Bias[0];
			fprintf(stderr, "Mean is given by %lf\n", Mean);
		}
		
		FILE* hist = fopen("hist.dat", "w");
		size_t s = 0;
		for (int i=0;i<256;i++) {
			s += Histogram[i];
			fprintf(hist, "%u %lf\t%u\t%zu\n", i, PBTrc->Traces.DecodeTable[0][i], Histogram[i], s);
		}
		fclose(hist);
		
		BaseCallsPulses_t BP;
		if (getBasecallsPulses(PBBas, &BP, holeN) != SUCCESS) {
			fputs("Error getting basecalls and pulses...\n", stderr);
			goto END;
		}
		
		unsigned int Tcount = 0U;
		unsigned int Gcount = 0U;
		unsigned int Acount = 0U;
		unsigned int Ccount = 0U;
		for (size_t i=0; i<BP.nBases; i++) {
			switch (BP.Bases[i]) {
				case 'T': Tcount++; break;
				case 'G': Gcount++; break;
				case 'A': Acount++; break;
				case 'C': Ccount++; break;
			}
		}
		fprintf(stderr, "PacBio Pulses are %u T's, %u G's, %u A's and %u C's\n", Tcount, Gcount, Acount, Ccount);
		
		char * const restrict Bases = (char*) malloc(nTimes*sizeof(char));
		if (Bases == NULL) {
			fputs("Cannot allocate memory for trace bases and types...\n", stderr);
			goto END;
		}
		memset(Bases, 'N', PBTrc->Traces.Length*sizeof(char));
		for (size_t i=0; i<BP.nBases; i++) {
			const unsigned int stop = BP.FramesEnd[i];
			const char TheBase = BP.Bases[i];
			for (unsigned int j=BP.FramesStart[i]; j<stop; j++) {
				Bases[j] = TheBase;
			}
		}
		
		Tcount = 0U;
		Gcount = 0U;
		Acount = 0U;
		Ccount = 0U;
		for (size_t i=0; i<nTimes; i++) {
			switch (Bases[i]) {
				case 'T': Tcount++; break;
				case 'G': Gcount++; break;
				case 'A': Acount++; break;
				case 'C': Ccount++; break;
			}
		}
		fprintf(stderr, "PacBio frames are %u T's, %u G's, %u A's and %u C's\n", Tcount, Gcount, Acount, Ccount);
#define TESTS 100
		unsigned char Thresholds[TESTS];
		float Coefs[TESTS];
		float coef = 0.0f;
		for (int i=0; i<TESTS; i++) {
			Coefs[i] = coef;
			Thresholds[i] = getIQRLimit(Histogram, nTimes, &(PBTrc->Traces.DecodeTable[0][0]), coef);
			coef += 0.1f;
		}
		
		
// 		Tcount = 0U;
// 		Gcount = 0U;
// 		for (size_t i=0; i<nTimes; i++) {
// 			if (Indices->T[i] >= Threshold && Bases[i] == 'T') Tcount++;
// 			if (Indices->T[i] < Threshold && Bases[i] == 'T') Gcount++;
// 		}
// 		fprintf(stderr, "T TP: %u\t FN: %u\n", Tcount, Gcount);
		
		ContengencyTable_t * CT = computeContengency(Indices->T, Bases, Thresholds, 'T', nTimes, TESTS);

		FILE * ROC = fopen("ROC.dat","w");
		for (int i=0; i<TESTS; i++) {
// 			fprintf(stderr, "TP : %10u\tFP : %10u   |   %u\n"
// 		                "FN : %10u\tTN : %10u   |   %u\n"
// 										"--------------------------------\n"
// 										"     %10u      %10u\n",
// 										CT[i].TP, CT[i].FP, CT[i].TP + CT[i].FP,
// 					          CT[i].FN, CT[i].TN, CT[i].FN + CT[i].TN,
// 										CT[i].TP+CT[i].FN, CT[i].FP + CT[i].TN);
			const float inv1 = (CT[i].TP+CT[i].FN) == 0 ? 1.0f : 1.0f/((float)(CT[i].TP+CT[i].FN));
			const float inv2 = (CT[i].FP+CT[i].TN) == 0 ? 1.0f : 1.0f/((float)(CT[i].FP+CT[i].TN));
			fprintf(ROC,"%lf\t%lf\t%lf\n", Coefs[i], 1.0f - (float)CT[i].FP*inv2, (float)CT[i].TP*inv1);
		}
		fclose(ROC);
		
		goto END;
		
		
		
		__m128* const restrict __tr = DecodeTracePositionIndex(PBTrc, Indices, holeN);
		
		
		
		
		
		
		FILE* out = fopen("Traces.dat","w");
		if (out != NULL) {
			float(* const restrict ftmp)[4] = (float (* const restrict)[4]) __tr; 
			char * const restrict Types  = (char*) malloc(PBTrc->Traces.Length*sizeof(char));
			if ( Types == NULL) {
				fputs("Cannot allocate memory for trace bases and types...\n", stderr);
				goto END;
			}
			memset(Types, 'S', PBTrc->Traces.Length*sizeof(char));
			
			HoleData_t HReg = HOLE_DATA_INIT;
			int res = getHoleRegions(PBBas, &HReg, holeN);
			if (res < 0) {
				fputs("Error getting Hole regions...\n", stderr);
				goto END;
			}
			
			/* Identify adapters */
			HoleRegion_t * restrict Regions = HReg.Regions;
			for (unsigned int r=0; r<HReg.nRegions; r++) {
				if (Regions[r].type == Adapter) {
					const unsigned int FramesStart = BP.FramesStart[Regions[r].start];
					const unsigned int FramesEnd   = BP.FramesEnd[Regions[r].stop];
					fprintf(stderr, "Setting type for adapter from [%u,%u] with quality %i\n",
									Regions[r].start, Regions[r].stop, Regions[r].quality);
					for (unsigned int k=FramesStart; k<FramesEnd; k++) Types[k] = 'A';
				}
			}
			
			/* Identify CAG repeats */
			fputs("Looking for CAG repeat profile cag.prf\t", stderr);
			FILE * dummy = fopen("cag.prf","r");
			if (dummy) {
				fputs("FOUND\n",stderr);
				fclose(dummy);
				
				struct Profile * const prf = (struct Profile *) _mm_malloc(sizeof(struct Profile), 16);
			  if (prf == NULL) {
			      fputs("Unable to allocate memory for the profile structure\n", stderr);
			      exit(1);
			  }
			  OutputVerbose = true;
				
			  /* Read the profile */
				{
				  const int ProfileCount = ReadProfile("cag.prf", prf, true, false);
				  if (ProfileCount != 1) {
				    fputs("Error found reading profile or multiprofile or PATTERN.\n", stderr);
				    exit(1);
				  }
				}
				
				/* Compute the alignments */
				struct Alignment * Alignments;
				prf->CutOff = 400;
				const int AlnCount = align(PBBas, prf, &Alignments, holeN);
					
				for (int l=0; l<AlnCount; l++) {
					fprintf(stderr, "Setting type for alignment of cag from [%u,%u] with score %i\n",
									Alignments[l].Matrix.row.Begin,
									Alignments[l].Matrix.row.End,
									Alignments[l].Score);
					const unsigned int FramesStart = BP.FramesStart[Alignments[l].Matrix.row.Begin];
					const unsigned int FramesEnd   = BP.FramesEnd[Alignments[l].Matrix.row.End];
					for (unsigned int k=FramesStart; k<FramesEnd; k++) Types[k] = 'R';
				}
				
				char fname[32];
				snprintf(fname, 32, "DNA_%i.txt", holeN);
				dummy = fopen(fname, "w");
				if (dummy) {
					int OutWidth = 0;
					size_t k=0;
					char previous = '\t';
					register char c;
					
					for(; k<nTimes; k++) {
						c = Bases[k];
						if (c != previous) {
							if (c != 'N') {
								if (OutWidth == 100) {OutWidth = 0; fputc('\n', dummy);}
								const char C = (Types[k] == 'S') ? c = c - 'A' + 'a' : c; 
								fputc(C, dummy);
								OutWidth++;
							}
							previous = c;
						}
					}
					
					fclose(dummy);
				}
				
				/* free memory */
				free(Alignments);
				FreeProfile(prf, true);
			}
			else {
				fputs("NOT FOUND\n",stderr);
			}
		
			for (size_t i=0; i<PBTrc->Traces.Length; i++) {
				fprintf(out, "%7lu\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%c\t%c\n",
				        i, ftmp[i][0], ftmp[i][1], ftmp[i][2], ftmp[i][3], Bases[i], Types[i]);
			}
			fclose(out);
			fprintf(stderr, "Hole %i traces exported to Traces.dat\n", holeN);
		
			out = fopen("StairTraces.dat","w");
			if (out) {
				union { float f; unsigned int u;} Ones = {.u=-1};
				const __m128 __Ones = _mm_load1_ps(&(Ones.f));
				const __m128 __Zero = _mm_setzero_ps();
				__m128 __Base         = __tr[0];
				__m128 __previousDiff = _mm_sub_ps(__tr[1], __tr[0]);
				
// 				{
// 					__m128 V1 = _mm_shuffle_ps(__Base, __Base, 0x0);
// 					__m128 V2 = _mm_shuffle_ps(__previousDiff, __previousDiff, 0x0);
// 					__m128 V3 = _mm_shuffle_ps(__Zero, __Zero, 0x0);
// 					__m128 __Values   = _mm_blend_ps(__tr[0], V1, 0x2);
// 					__Values   = _mm_blend_ps(__Values, V2, 0x4);
// 					__Values   = _mm_blend_ps(__Values, V3, 0x8);
// 					_mm_store_ps(&ftmp[0][0], __Values);
// 				}
				const size_t Ni = PBTrc->Traces.Length-1;
				for (size_t i=1; i<Ni; i++) {
					__m128 __isEqual  = _mm_cmpeq_ps(__tr[i+1],__tr[i]);
					__m128 __nextDiff = _mm_sub_ps(__tr[i+1],__tr[i]);
					__nextDiff        = _mm_blendv_ps(__nextDiff, __previousDiff, __isEqual);
					__m128 __diffSign = _mm_xor_ps(__previousDiff, __nextDiff);
					__m128 __mask     = _mm_cmpgt_ps(__diffSign, __Zero);
					__m128 __Values   = _mm_blendv_ps(__tr[i], __Base, __mask);
					__Base            = _mm_blendv_ps(__Base, __tr[i], _mm_xor_ps(__mask, __Ones));
					
// 					__m128 V1 = _mm_shuffle_ps(__Base, __Base, 0x0);
// 					__m128 V2 = _mm_shuffle_ps(__previousDiff, __previousDiff, 0x0);
// 					__m128 V3 = _mm_shuffle_ps(__nextDiff, __nextDiff, 0x0);
// 					__Values   = _mm_blend_ps(__Values, V1, 0x2);
// 					__Values   = _mm_blend_ps(__Values, V2, 0x4);
// 					__Values   = _mm_blend_ps(__Values, V3, 0x8);
					_mm_store_ps(&ftmp[i][0], __Values);
					
					__previousDiff = __nextDiff; 
				}
				
				for (size_t i=0; i<PBTrc->Traces.Length; i++) {
					fprintf(out, "%7lu\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%c\t%c\n",
					        i, ftmp[i][0], ftmp[i][1], ftmp[i][2], ftmp[i][3], Bases[i], Types[i]);
				}
				fclose(out);
			}
			
// 			_mm_free(__tr);
			
			
			if (argc == 4) {
				const size_t WindowSize = (size_t) atoi(argv[3]);
				char fname[32];
				snprintf(fname,32,"Mean_%u_%zu.dat", holeN, WindowSize);
				out = fopen(fname, "w");
				if (out) {
					const float one = 1.0f;
					float values[4] __attribute__((aligned(16)));
					const int TotalOut = 
					fprintf(stderr, "Exporting Mean window of size %zu to %s\n", WindowSize, fname);
					const __m128i __itmp = _mm_set1_epi32(WindowSize);
					const __m128 __ftmp = _mm_cvtepi32_ps(__itmp);
					const __m128 __invN = _mm_div_ps(_mm_load1_ps(&one),__ftmp);
					
					register __m128 __Sum = _mm_setzero_ps();
					register size_t k;
					for (k=0; k<WindowSize; k++) __Sum = _mm_add_ps(__Sum, __tr[k]);
					
					__m128 __res = _mm_mul_ps(__Sum, __invN);
					_mm_store_ps(&values[0], __res);
					fprintf(out, "%lf\t%lf\t%lf\t%lf\n", values[0], values[1], values[2], values[3]);
					__m128 __first = __tr[0];
					for (;k<nTimes;k++) {
						__m128 __diff = _mm_sub_ps(__tr[k], __first);
						__Sum = _mm_add_ps(__Sum, __diff);
						__first = __tr[k-WindowSize];
						_mm_store_ps(&values[0], _mm_mul_ps(__Sum, __invN));
						fprintf(out, "%lf\t%lf\t%lf\t%lf\n", values[0], values[1], values[2], values[3]);
					}
					
					fclose(out);
				}
			}
			
			unsigned int Stats[256][256][26];
			
			unsigned char * restrict IndicesPtr = Indices->memory;
			unsigned int * restrict Accounting = (unsigned int * const ) __tr;
			memset(Accounting, 0, 4*nTimes*sizeof(unsigned int));
			for (int i=0; i<4; i++) {
				memset(Stats, 0, 256*256*26*sizeof(unsigned int));
				register unsigned char previous = IndicesPtr[0]; 
				unsigned int previousBase;
				switch(Bases[i]) {
					case 'T': previousBase = 0; break;
					case 'G': previousBase = 1; break;
					case 'A': previousBase = 2; break;
					case 'C': previousBase = 3; break;
					case 'N': previousBase = 4; break;
				}
				for (size_t t=1; t<nTimes; t++) {
					const unsigned char Id = IndicesPtr[t];
					Stats[previous][Id][0]++;
					unsigned int currentBase;
					switch(Bases[t]) {
						case 'T': currentBase = 0; break;
						case 'G': currentBase = 1; break;
						case 'A': currentBase = 2; break;
						case 'C': currentBase = 3; break;
						case 'N': currentBase = 4; break;
					}
					
					Stats[previous][Id][1+previousBase*5+currentBase]++;
					previous = Id;
					previousBase = currentBase;
				}
				char FName[16];
				snprintf(FName, 16, "Stat_%i.dat", i);
				out = fopen(FName, "w");
				if (out != NULL) {
					for (int l=0; l<256; l++) {
						for (int m=0;m<256;m++) {
							fprintf(out, "%u %u ", l, m);
							for (int n=0; n<26; n++) fprintf(out, "%u ", Stats[l][m][n]);
							fputc('\n', out);
						}
						fputc('\n', out);
					}
					fclose(out);
				}
				
				/* Compute cumulative sum */
				unsigned int * const restrict Cumulative = (unsigned int *const) &Stats[0][0][0]; 
				memset(Cumulative, 0, 256*sizeof(unsigned int *));
				for (size_t t=0; t<nTimes; t++) Cumulative[IndicesPtr[t]]++;
				unsigned int Sum = 0U;
				for(size_t i=0; i<256; i++) {
					Sum += Cumulative[i];
					Cumulative[i] = Sum;
				}
				
				for (size_t t=0; t<nTimes; t++) Accounting[t] = Cumulative[IndicesPtr[t]];
				
				Accounting += nTimes;
				IndicesPtr += Indices->stride;
			}
			
			out = fopen("Traces_cumulative.dat", "w");
			if (out != NULL) {
				unsigned int * restrict Accounting = (unsigned int * const ) __tr;
				for (size_t i=0; i<PBTrc->Traces.Length; i++) {
					fprintf(out, "%7lu\t%10u\t%10u\t%10u\t%10u\t%c\n",
					        i, Accounting[i], Accounting[i+nTimes], Accounting[i+2*nTimes], Accounting[i+3*nTimes], Bases[i]);
				}
				fclose(out);
				fprintf(stderr, "Hole %i cumulative traces exported to Traces_cumulative.dat\n", holeN);
			}
			
			_mm_free(__tr);
			free(Bases);
			free(Indices);
		}
	}
		
	
	
	if (PacBioClose(PBTrc) != SUCCESS) {
		fprintf(stderr, "Error closing Pac Bio Trace file...\n");
	}
	if (PacBioClose(PBBas) != SUCCESS) {
		fprintf(stderr, "Error closing Pac Bio Bas file...\n");
	}
	
	exit(0);
	
	END:
	if (PacBioClose(PBTrc) != SUCCESS) {
		fprintf(stderr, "Error closing Pac Bio Traces...\n");
	}
	if (PacBioClose(PBBas) != SUCCESS) {
		fprintf(stderr, "Error closing Pac Bio Basecalls...\n");
	}
	exit(1);
}

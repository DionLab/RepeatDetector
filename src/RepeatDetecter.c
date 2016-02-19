/***************************************************************************************************
                                        PACBIO REPEAT
 ***************************************************************************************************
  Apr 1, 2016 pfrepeat.c
 ***************************************************************************************************
 (C) 2016 Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@isb-sib.ch)
 ***************************************************************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#ifdef __USE_WINAPI__
# include <windows.h>
#else
# include <pthread.h>
#endif
#ifdef HAVE_ALLOCA_H
# include <alloca.h>
#endif
#include <stdint.h>
#include <getopt.h>
#ifdef __USE_AFFINITY__
# include <unistd.h>
# include <sched.h>
#endif
#include <pthread.h>
#include <nmmintrin.h>
#include "Version.h"
#include "pfVersion.h"
#include "pfProfile.h"
#ifdef PRF_CORE_PCRE
#include "pfRegexp.h"
#endif
#include "pfInput.h"
#include "pfOutput.h"
#include "pfDispatch.h"
#include "system.h"

#define HEADER \
         "%----------------------------------------------------------------------------%\n"\
	       "|                          RepeatDecoder v" PROJECT_VERSION "                       |\n"\
	       "|                        built with libPRF v" PRF_VERSION "                     |\n"\
	       "%----------------------------------------------------------------------------%\n"\
	       "| Built on " __DATE__ " at " __TIME__ ".                                          |\n"
				 
static SystemInfo System;
#ifdef __USE_AFFINITY__
Affinity_Mask_t * Thread_masks[2] = {0,0};						/* Define variables to hold thread affinity mask */
unsigned int Thread_count[2] = {0,0};
pthread_attr_t * restrict threads_attr = NULL;
#endif

#define TEST_IF_AVAILABLE(TYPE, option) if (OutputType.Type != TYPE) {\
	fprintf(stderr,"Invalid option for this output type: " option ".\n");\
	exit(1);\
}

static const char opt_to_test[] = "P:ac:H:I:i:t:hs0:o:Vw>:!:X:FD:em,Q"
#ifdef __USE_AFFINITY__
	"012:3"
#endif
# ifdef PRF_OUTPUT_PDF
	"nN"
#endif
#ifdef PRF_OUTPUT_GRAPHICS
	"W:H:S4:5:6:"
#endif
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
	"jJK:yp:L:qYTM"
#endif
#if defined(PRF_OUTPUT_PDF) || defined(PRF_OUTPUT_GRAPHICS)
	"kr:"
#endif
#ifdef PRF_CORE_PCRE
	"R:z"
#endif
;

static const struct option long_options[] =
{
  /*
	 * These options set a flag.
	 */

  /*
	 * These options don't set a flag. We distinguish them by their indices.
	 */
	{"help",               		no_argument,       	0,	'h'},
	{"sse2",									no_argument,				0,	's'},
	{"verbose",								no_argument,				0,	'V'},
	/* Profile */
	{"prf",										required_argument,	0,	'P'},
	{"cutoff",								required_argument,	0,	'c'},
	{"optimal", 							no_argument,				0,	'a'},
	{"bestof",								no_argument,				0,	','},
	/* Sequence */
	{"sequence-clip",					required_argument,	0,	'D'},
	{"score-range",						required_argument,	0,	'!'},
	{"cycle-range",						required_argument,	0,	'>'},
	{"use-cycle",							no_argument,				0,	'm'},
	{"with-revcomp",					no_argument,				0,	'e'},
	{"not",										no_argument,				0,	'Q'},
#ifdef PRF_CORE_PCRE
	/* Regex */
	{"regex",									required_argument,	0,	'R'},
	{"cycle-count",						no_argument,				0,	'z'},
#endif
	/* Database indexing options */
	{"create-index-database",	required_argument,	0,	'I'},
	{"use-index-database",		required_argument,	0,	'i'},
	{"stream-fasta",					no_argument,				0,	'.'},
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
	/* PacBio */
	{"subreads",							no_argument,				0,	'j'},
	{"consensus",							no_argument,				0,	'J'},
	{"polymerase",						no_argument,				0,	'x'},
	{"best-subread",					no_argument,				0,	'M'},
	{"invalid-zmw",						no_argument,				0,	'O'},
	{"zmw-number",						required_argument,	0,	'K'},
	{"not-only-sequencing",		no_argument,				0,	'y'},
	{"min-filter-score",			required_argument,	0,	'p'},
	{"max-filter-score",			required_argument,	0,	'L'},
	{"keep-invalid-zmw",			no_argument,				0,	'q'},
	{"discard-filter",				no_argument,				0,	'Y'},
	{"test-filter",						no_argument,				0,	'T'},
#endif
	/* Others */
	/* SMP options*/
	{"nthreads",							required_argument,	0,	't'},
#ifdef __USE_AFFINITY__
	{"no-affinity",						no_argument,				0,	'0'},
	{"split", 								no_argument,				0,	'1'},
	{"thread-affinity",				required_argument,	0,	'2'},
	{"no-shared-core",				no_argument,				0,	'3'},
#endif
	/* Print output methods*/
	{"output",								required_argument,	0,	'o'},
	{"output-name",						required_argument,	0,	'u'},
	{"best",									no_argument,				0,	'A'},
	{"with-reverse",					no_argument,				0,	'k'},
#ifdef PRF_OUTPUT_GRAPHICS
	{"png-pixel-width",				required_argument,	0,	'W'},
	{"png-pixel-height",			required_argument,	0,	'H'},
	{"scale-on-image",				no_argument,				0,	'S'},
	{"scale-min",							required_argument,	0,	'4'},
	{"scale-max",							required_argument,	0,	'5'},
	{"compress",							required_argument,	0,	'6'},
#endif
#ifdef PRF_OUTPUT_PDF
	{"with-path",							no_argument,				0, 	'n'},
	{"with-graphs",						no_argument,				0,	'N'},
	{"region",								required_argument,	0,	'r'},
	{"whole-sequence",				no_argument,				0,	'w'},
#endif
	/* Print output methods*/
	{"output-method",					required_argument,	0,	'o'},
	{"output-width",					required_argument,	0,	'X'},
	{"source", 								no_argument,				0,	'F'},
	{"before",								no_argument,				0,	'b'},
	{"after",									no_argument,				0,	'B'},
	{"in-between",						no_argument,				0,	'C'},
	{0, 0, 0, 0}
};

char OutputFileName[256] = "Results";

const void *OutputMethodOptions;
const _Bool CompleteCycleoOnly = false;
int PlotAllSequence = 0;
_Bool WithReverseComplement = false;

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
#ifndef PRF_CORE_PCRE
	" RepeatDecoder [--prf file] [options] [Database, - ]\n"
#else
	" RepeatDecoder [--prf file|--regex {pattern}] [options] [Database, - ]\n"
#endif
	" Options:\n"
	"  Sequence\n"
	"    --sequence-clip [b:e]          : constrain alignement to start after b\n"
	"                                     and to end before e\n"
	"    --score-range [a:b]            : contrain score within given range\n"
	"    --cycle-range [a:b]            : contrain repeat cycle within given range\n"
	"    --not                          : inverse above criteria\n"
	"    --with-revcomp                 : test also reverse complement\n"
	"    --bestof                       : when used with --with-revcomp and --optimal,\n"
	"                                     output only the best of both orientation\n"
	"                                     this option is for text output only, in\n"
	"                                     histogram mode when reverse complement is\n"
	"                                     required, best of is implicit.\n\n" 
	"  Profile\n"
	"    --cutoff                  [-c] : score cutoff value\n"
	"    --optimal                 [-a] : output only optimal alignment per sequence\n\n"
#ifdef PRF_CORE_PCRE
	"  Regex\n"
	"   --cycle-count                   : account cycle of regex pattern till\n"
	"                                     histogram size (--histogram-bins)\n\n"
#endif
	"  Database\n"
	"   FASTA\n"
	"    Fasta format allows to be piped with the '-' symbol as database file name\n"
	"     --create-index-database  [-I] : output indices to given file\n"
	"     --use-index-database     [-i] : use indices stored in given file\n"
	"     --stream-fasta                : stream the file rather than indexing\n\n"
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
#if !defined(PRF_INPUT_HDF5)
	"   PacBio BAM\n"
#elif !defined(PRF_INPUT_PBBAM)
	"   PacBio HDF5\n"
#else
	"   PacBio HDF5 or BAM\n"
#endif
	"     --subreads                    : work on subreads (default)\n"
	"     --consensus                   : work on consensus\n"
	"     --polymerase                  : work on entire polymerase data\n"
	"     --invalid-zmw                 : work on invalid zmw only\n"
	"     --best-subread                : keep only best subread of zmw\n"
	"     --zmw-number i                : limit to zmw i\n"
	"     --not-only-sequencing         : compute on any zmw type\n"
	"     --min-filter-score i          : minimal High Quality region score i to keep data\n"
	"     --max-filter-score i          : maximum High Quality region score i to keep data\n"
	"     --keep-invalid-zmw            : keep zmw removed by PacBio\n"
	"     --discard-filter              : does not filter or trim according to HQ region\n"
	"     --test-filter                 : output regions based upon the above filtering\n\n"
#endif
#ifdef PRF_INPUT_RANDOM
	"   Randomly generated sequences\n"
	"     Random[N,L,S]                 : provide this as database name where\n"
	"                                     N is the number, L the sequence length\n"
	"                                     and S the generator seed\n\n"
#endif
	" Optimizations\n"
// 	"   --sse2                     [-s] : enforces SSE 2 only instruction set\n"
	"   --nthreads                 [-t] : max number of threads to use\n"
#ifdef __USE_AFFINITY__
	"   --no-affinity                   : disable CPU affinity file\n"
	"   --thread-affinity               : file containing thread mask,\n"
	"                                     one row for one thread\n"
	"   --no-shared-core                : Prevent core resource sharing\n"
	"   --split                         : if both SSE 2 & 4.1 are available,\n"
	"                                     split half-half using linked resources\n\n"
#else
	"\n"
#endif
	"  Output\n"
	"    --output-name                  : specify output file name\n"
	"    --output type             [-o] : output in format type\n"
	"    --output-width            [-X] : text output width\n"
	"   type is:\n"
#ifdef PRF_OUTPUT_FORMAT_TEST
	"    - test                         : compare alignment with sequence header\n"
#endif
	"    - dat                          : ASCII matrix dump\n"
	"    - tab                          : tab separated matrix dump\n"
	"    - histogram                    : generate an histogram of number of cycles\n"
	"      --score-range                :\n"
	"      --cycle-range                :\n"
	"      --use-cycle                  : in case both score and cycle range are\n"
	"                                   : provided, use cycle rather than score\n"
	"    - density                      : generate a 2D histogram, score vs cycles\n"
	"      --score-range and --cycle-range are mandatory\n"
#ifdef PRF_OUTPUT_FORMAT_XPSA
	"    - xPSA\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_TSV
	"    - TSV\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_FASEARCH
	"    - FAsearch\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_INCMATCH
	"    - IncMatch\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_INTERPRO
	"    - Interpro\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_PSMAKER
	"    - PSMaker\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_PFSCAN
	"    - PfScan\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_ONELINE
	"    - OneLine\n"
#endif	
#ifdef PRF_OUTPUT_FORMAT_SIMPLE
	"    - Simple\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_FASTA
	"    - FASTA                        : FASTA file format\n"
	"      --source                     : output entire source sequence\n"
	"      --before\n"
	"      --after\n"
	"      --in-between\n"
#endif
#ifdef PRF_OUTPUT_FORMAT_HIGHLIGHT
	"    - Highlight                    : Highlighted FASTA output\n"
#endif
#ifdef PRF_OUTPUT_PDF
	"    - pdf                          : pdf outout\n"
	"      --whole-sequence        [-w] : does not limit output to alignment\n"
	"      --with-path                  : does output paths\n"
	"      --with-graphs                : does output side graphs\n"
	"      --region [a,b][c,d]     [-r] : specify region to output\n"
#endif
#ifdef PRF_OUTPUT_GRAPHICS
	"    - png                          : png output\n"
	"      --png-pixel-width       [-W] : set pixel width\n"
	"      --png-pixel-height      [-H] : set pixel height\n"
	"      --scale-on-image        [-S] : print scale on data image\n"
	"      --scale-min                  : lower scale score limit\n"
	"      --scale-max                  : upper scale score limit\n"
	"      --region [x0,y0][x1,y1]      : specify region to output, x=profile, y=sequence\n"
	"      --compress [x,y]             : output a pixel per rectangle\n"
	"                                     size given by [x,y]\n"
#endif
	"  Other\n"
	"    --verbose                 [-V] : verbose on stderr\n"
	"    --help                    [-h] : output command help\n\n"
	"This is version " PROJECT_VERSION ".\n"
	"It uses libPRF version " PRF_VERSION ".\n",
	stream);
  exit(1);
}

int main (int argc, char *argv[])
{
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL STRUCTURES
  ////////////////////////////////////////////////////////////////////////////////////////////////

  struct Profile * restrict prf = NULL; 		/* Profile */
  struct Profile * * restrict prfs = NULL; 	/* Profiles */
  FASTAStructure FASTA = { 0 };							/* Sequence Database File */
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
  PacBio_t * restrict PBS = NULL;						/* PacBio HDF5 datafile */
#endif
  struct timeval _t0, _t1;									/* Timing structures */
	struct RegEx * restrict regex = NULL;			/* Regex structure */

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL DATA
  ////////////////////////////////////////////////////////////////////////////////////////////////

  size_t nCPUs=0;										/* number of threads */
  int res;
  int Cutoff = -1;									/* Default Coverage cutoff from command line, if not zero then enforces that value */
  char * ProfileFile = NULL;				/* Profile file */
  char * DB = NULL;									/* FASTA sequence file */
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
  int isPacBio = 0;									/* Are we crunching HDF5 file from PacBio ?*/
  PacBioDispatchOptions_t PacBioDispatchOptions = PB_DEFAULT_OPTIONS;
	unsigned int PacBioZmWNumber;
#endif
	enum ComplementaryMode CplmMode = NONE;	/* Complementary mode to use (DNA, IUPAC,etc...)*/
  size_t * shares = 0;
  size_t MaxProfileLength = 0;						/* Maxiumum profile length when several are given */
  struct ThreadArrayData *threads_array_arg = NULL;
  enum Version ComputeVersion = SSE41;		/* Trigger SSE version to use for filter and alignment */
  OutputType_t OutputType = OUTPUT_TYPE_INIT; OutputMethodOptions = &(OutputType.Specific);
	_Bool FastaStream = false;
#ifdef PRF_INPUT_RANDOM
	RandomData_t RD = { .seed = 0U, .Length = 0U, .N = 0UL, .nCAGSuffix=0U };	/* Random data generator parameters */
#endif
#ifdef PRF_CORE_PCRE
	_Bool RegexInFile = false;					/* Are regex provided on stdin or in file */
	char * RegexSource = NULL;					/* Pointer to the regex source argument */
#endif
#ifdef __USE_AFFINITY__
  char buffer[128] __attribute__((aligned(16)));	/* buffer to read affinity file mask */
  _Bool noAffinity = false;												/* disable use of cpu affinity */
  _Bool split = false;
  _Bool noSharedCore = false;					/* Prevent hyperthreading or AMD compute unit to share resources */
  _Bool GivenAffinityFile = false;		/* File holding a mask for each thread */
  char * AffinityMaskFileName;				/* Name of affinity mask file provided by option m */
#endif

	////////////////////////////////////////////////////////////////////////////////////////////////
  // INITIALIZE VALUES
  ////////////////////////////////////////////////////////////////////////////////////////////////
 
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // SYSTEM ARCHITECTURE ANALYSIS
  ////////////////////////////////////////////////////////////////////////////////////////////////
  getSystemInfo(&System);

  /* Check for minimum requirement */
  if (!(System.Extensions & MM_SSE2)) {
      fputs("pfrepeat requires at least a CPU capable of SSE 2.\n", stderr);
      exit(1);
  }

  /* Allow fast SSE 4.1 extensions ? */
  if (System.Extensions & MM_SSE41) {
      ComputeVersion = SSE41;
  }
  else {
      fputs("Sorry only compute methods based upon SSE 4.1 are available yet\n", stderr);
      exit(1);
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // OPTIONS
  ////////////////////////////////////////////////////////////////////////////////////////////////

	while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
			case 'P':
				OutputType.SearchWith |= PRF;
				break;
			case 'R':
				OutputType.SearchWith |= REGEX;
				break;
			case 'o':
				if (strncmp(optarg, "dat", 3) == 0) {
					OutputType.OutputFct = DataOutput;
					OutputType.Type = DATA;
				}
				else if (strncmp(optarg, "tab", 3) == 0) {
					OutputType.OutputFct = TabOutput;
					OutputType.Type = DATA;
				}
				else if (strncmp(optarg, "histogram", 9) == 0) {
					OutputType.Type = HISTOGRAM;
				}
				else if (strncmp(optarg, "density", 7) == 0) {
					OutputType.Type = DENSITY;
				}
#ifdef PRF_OUTPUT_PDF
				else if (strncmp(optarg, "pdf", 3) == 0) {
					OutputType.Type = PDF;
					OutputType.OutputFct = PDFOutput;
				}
#endif
#ifdef PRF_OUTPUT_GRAPHICS
				else if (strncmp(optarg, "png", 3) == 0) {
					OutputType.Type = PNG;
					OutputType.OutputFct = PNGOutput;
				} 
#endif
#ifdef PRF_OUTPUT_FORMAT_TEST
				else if (strncmp(optarg, "test", 4) == 0) {
					OutputType.Type = TEST;
				}
#endif
				else {
					OutputType.Type = TEXT;
				}
				if (OutputType.Type != TEXT) memset(&OutputType.Specific, 0, sizeof(union dummy));
				break;
			case 'V':
				OutputVerbose = true;
				break;
      case 'h':
				Usage(stderr);
				break;
    }
  }
  
	if (OutputVerbose) {
   fputs(HEADER
#ifdef __USE_MMAP__
	       "| Using Linux kernel MMAP function.                                          |\n"
#endif
				 "%----------------------------------------------------------------------------%\n"
      ,stderr);
    printSystemInfo(stdout, &System);
#ifdef USE_32BIT_FORMAT
    fputs("Using 32 bit format integer for scores\n", stderr);
#endif
    if (ComputeVersion == SSE2 && (System.Extensions & MM_SSE41)) {
			fputs("Enforcing SSE 2...\n", stderr);
    }
  }
  
  optind = 1;
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
			/****************************************** Profile  *****************************************/
			case 'P':
				ProfileFile = optarg;
				break;
			case 'c':
				Cutoff = atoi(optarg);
				break;
			/****************************************** Sequence *****************************************/
			case 'D':
				{
					if (sscanf(optarg, "[%hu:%hu]", &OutputType.BorderClip[0], &OutputType.BorderClip[1]) != 2) {
						fprintf(stderr, "Unable to identify sequence clip: %s\n", optarg);
						exit(1);
					}
				}
				OutputType.Constrains |= BORDER_CLIP;
				break;
			case 'e':
				OutputType.Constrains |= WITH_REVERSE;
				break;
			case '!':
				if (sscanf(optarg, "[%i:%i]", &(OutputType.ScoreRange[0]), &(OutputType.ScoreRange[1])) != 2) {
					fprintf(stderr, "Unable to read range given by %s\n", optarg);
					exit(1);
				}
				if (OutputVerbose) {
					fprintf(stderr, "Score range set to [%i:%i]\n", 
									OutputType.ScoreRange[0], OutputType.ScoreRange[1]);
				}
				OutputType.Constrains |= SCORE_RANGE; 
				break;
			case '>':
				{
					int tmp[2];
					if (sscanf(optarg, "[%i:%i]", &(tmp[0]), &(tmp[1])) != 2) {
						fprintf(stderr, "Unable to read range given by %s\n", optarg);
						exit(1);
					}
					OutputType.CycleRange[0] = tmp[0] < 0 ? -tmp[0] : tmp[0];
					OutputType.CycleRange[1] = tmp[1] < 0 ? -tmp[1] : tmp[1];
					if (OutputVerbose) {
						fprintf(stderr, "Cycle range set to [%i:%i]\n", 
										OutputType.CycleRange[0], OutputType.CycleRange[1]);
					}
					OutputType.Constrains |= CYCLE_RANGE; 
				}
				break;
			case 'Q':
				OutputType.Constrains |= INVERSE_SELECTION;
				break;
			case 'm':
				TEST_IF_AVAILABLE(HISTOGRAM, "--use-cycle");
				OutputType.Specific.Histogram.CycleRatherThanScore = true;
				break;
#ifdef __USE_AFFINITY__
      case '0':
        noAffinity = true;
        break;
      case '3':
				noSharedCore = true;
				break;
      case '1':
				if (ComputeVersion == SSE41) {
					split = true;
				} else {
					fputs("Split not possible without SSE 4.1\n", stderr);
					exit(1);
				}
				break;
      case '2':
				GivenAffinityFile = true;
				AffinityMaskFileName = optarg;
				break;
#endif
      case 'V':
				break;
      case 't':
				nCPUs = (size_t) atoi(optarg);
				break;
      case 's':
				ComputeVersion = SSE2;
				break;
      case 'I':
				FASTA.Options |= DoFastaIndexExport;
				FASTA.indexFileName = optarg;
				break;
      case 'i':
				FASTA.Options |= DoFastaIndexImport;
				FASTA.indexFileName = optarg;
				break;
			case '.':
				FastaStream = true;
				break;
			case 'X':
				OutputPrintWidth = (unsigned int) atoi(optarg);
				break;
				/************************************** TEXT specific options ******************************/
			case 'a':
				TEST_IF_AVAILABLE(TEXT, "--optimal");
				OutputType.Specific.Text.OptimalOnly = true;
				break;
			case ',':
				TEST_IF_AVAILABLE(TEXT, "--bestof");
				OutputType.Specific.Text.BestOfStdAndRevComp = true;
				break;
			case 'F':
				TEST_IF_AVAILABLE(TEXT, "--source");
				OutputType.Specific.Text.Range = Wrapper_Source;
				break;
			case 'b':
				TEST_IF_AVAILABLE(TEXT, "--before");
				OutputType.Specific.Text.Range = Wrapper_Before;
				break;
			case 'B':
				TEST_IF_AVAILABLE(TEXT, "--after");
				OutputType.Specific.Text.Range = Wrapper_After;
				break;
			case 'C':
				TEST_IF_AVAILABLE(TEXT, "--in-between");
				OutputType.Specific.Text.Range = Wrapper_InBetween;
				break;
			case 'n':
#ifdef PRF_OUTPUT_PDF
				if (OutputType.Type == PDF)	OutputType.Specific.PDF.WithPaths = true;
#endif
				break;
			case 'N':
#ifdef PRF_OUTPUT_PDF
				if (OutputType.Type == PDF)	OutputType.Specific.PDF.WithProfileScoreGraph = true;
#endif
				break;
			case 'k':
#ifdef PRF_OUTPUT_PDF
				if (OutputType.Type == PDF)	OutputType.Specific.PDF.WithReverse = true;
#endif
#ifdef PRF_OUTPUT_GRAPHICS
				if (OutputType.Type == PNG) OutputType.Specific.PNG.WithReverse = true;
#endif
				break;
			case 'w':
				TEST_IF_AVAILABLE(PDF, "--whole-sequence");
				PlotAllSequence = 1;
				break;
			case 'u':
				strncpy(OutputFileName, optarg, 256);
				break;
			case 'o':
				if (OutputVerbose) fprintf(stderr,"Output string is '%s'\n", optarg);
				if (OutputType.Type == TEXT) {
#ifdef PRF_OUTPUT_FORMAT_INTERPRO
 					if (strncasecmp(optarg, "Interpro", 8) == 0 ){
 						OutputType.Specific.Text.Print = PrintInterpro;
 						goto Found;
 					}
#endif
#ifdef PRF_OUTPUT_FORMAT_PSMAKER
					if (strncasecmp(optarg, "PSMaker", 7) == 0 ){
						OutputType.Specific.Text.Print = PrintPSMaker;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_TSV
					if (strncasecmp(optarg, "TSV", 3) == 0 ){
						OutputType.Specific.Text.Print = PrintTSV;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_XPSA
					if (strncasecmp(optarg, "xPSA", 4) == 0 ){
						OutputType.Specific.Text.Print = PrintxPSA;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_FASEARCH
					if (strncasecmp(optarg, "FAsearch", 8) == 0 ){
						OutputType.Specific.Text.Print = Printfasearch;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_ONELINE
					if (strncasecmp(optarg, "OneLine", 4) == 0 ){
						OutputType.Specific.Text.Print = PrintOneLine;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_INCMATCH
					if (strncasecmp(optarg, "IncMatch", 4) == 0 ){
						OutputType.Specific.Text.Print = PrintIncMatch;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_PFSCAN
					if (strncasecmp(optarg, "PfScan", 4) == 0 ){
						OutputType.Specific.Text.Print = PrintPfscan;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_PFSCAN
					if (strncasecmp(optarg, "Simple", 4) == 0 ){
						OutputType.Specific.Text.Print = PrintSimple;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_FASTA
					if (strncasecmp(optarg, "FASTA", 5) == 0){
						OutputType.Specific.Text.Print = PrintFASTA;
						goto Found;
					}
#endif
#ifdef PRF_OUTPUT_FORMAT_HIGHLIGHT
					if (strncasecmp(optarg, "Highlight", 9) == 0) {
						OutputType.OutputFct = HighlightedFASTAOutput;
						goto Found;
					}
#endif

					fprintf(stderr, "Output format %s not recognized!\n", optarg);
					exit(1);
					Found: ;
				}
				break;
#if defined(PRF_INPUT_HDF5) || defined(PRF_INPUT_PBBAM)
			case 'j':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_SUBREADS;
				break;
			case 'J':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_CONSENSUS;
				break;
			case 'x':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_ZMW;
				break;
			case 'O':
				PacBioDispatchOptions.Selection |= PB_DISPATCH_INVALID;
				break;
			case 'K':
				{
					const int itmp = atoi(optarg);
					if (itmp < 0 ) {
						fprintf(stderr, "PacBio hole number must be positive or null (%i)\n",itmp);
						exit(1);
					}
					else {
						PacBioZmWNumber = (unsigned int) itmp;
						PacBioDispatchOptions.ZMW = &PacBioZmWNumber;
						PacBioDispatchOptions.nZMW = 1;
					}
				}
				break;
			case 'p':
				PacBioDispatchOptions.minReadAccuracy = atoi(optarg);
				PacBioDispatchOptions.Selection |= PB_HAS_FILTER_SCORE;
				break;
			case 'L':
				PacBioDispatchOptions.maxReadAccuracy = atoi(optarg);
				PacBioDispatchOptions.Selection |= PB_HAS_FILTER_SCORE;
				break;
			case 'q':
				PacBioDispatchOptions.Selection |= PB_KEEP_INVALID;
				break;
			case 'Y':
				PacBioDispatchOptions.Selection |= PB_DISCARD_FILTER;
				break;
			case 'y':
				PacBioDispatchOptions.Selection |= PB_NOT_ONLY_SEQUENCING;
				break;
			case 'T':
				PacBioDispatchOptions.Selection |= PB_TEST_OUTPUT;
				break;
			case 'M':
				PacBioDispatchOptions.Selection |= PB_BEST_SUBREAD_OF_ZMW;
				break;
#endif
#ifdef PRF_CORE_PCRE
			case 'R':
				RegexSource = optarg;
				break;
			case 'z':
				RegexCycleCount = true;
				break;
#endif
#if defined(PRF_OUTPUT_PDF) || defined(PRF_OUTPUT_GRAPHICS)
			case 'r':
				if (OutputType.Type == PDF) {
					const int count = sscanf(optarg, "[%u,%u][%u,%u]", &OutputType.Specific.PDF.x[0], &OutputType.Specific.PDF.y[0],
						                                                 &OutputType.Specific.PDF.x[1], &OutputType.Specific.PDF.y[1]);
					OutputType.Specific.PDF.GivenRegion = true;
					if (count != 4) {
							fprintf(stderr, "Error reading region range, only read %i element from %s\n", count, optarg);
							exit(1);
					}
				}
				
				break;
#endif
      case 'h':
				Usage(stderr);
				break;
			case '?':
				fprintf(stderr, "Unknown option %s\n", argv[option_index]);
				exit(1);
			default:
				fprintf(stderr,"Option %c is unknown\n", c);
    }
  }

  if (optind >= argc) {
    fputs("Expected arguments after options\n", stderr);
    Usage(stderr);
  }
  else {
    DB = argv[optind];
  }
  
  if (OutputType.Type == HISTOGRAM || OutputType.Type == DENSITY) {
			OutputType.Specific.Histogram.BaseFileName = OutputFileName;
	}
	else if (OutputType.Type == DATA) {
			OutputType.Specific.Data.BaseFileName = OutputFileName;
	}
#if defined(PRF_OUTPUT_PDF) || defined(PRF_OUTPUT_GRAPHICS)
	else if (OutputType.Type == PDF) {
			OutputType.Specific.PDF.BaseFileName = OutputFileName;
	}
#endif
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // COMMAND LINE COHERENCY
  ////////////////////////////////////////////////////////////////////////////////////////////////
  {
		_Bool AllOK = true;
		
		if (OutputType.SearchWith == 0) {
#ifdef PRF_CORE_PCRE
			fputs("Expecting a profile with --prf or a regex with --regex {}.\n", stderr);
#else
			fputs("Expecting a profile with --prf.\n", stderr);
#endif
			AllOK = false;
		}

		if (!DB) {
			fputs("Expected arguments after options\n", stderr);
			AllOK = false;
		}
  
		if (!AllOK) Usage(stderr);
	}
	
#ifdef PRF_CORE_PCRE	
	if ( OutputType.SearchWith == REGEX && OutputType.Type == DENSITY) {
		fprintf(stderr, "Density plot not possible with regex!\n");
		res = 1;
		goto bail;
	}
#endif

	if (OutputType.Type == HISTOGRAM) {
		if ((OutputType.Constrains & CYCLE_RANGE) && !(OutputType.Constrains & SCORE_RANGE))
			OutputType.Specific.Histogram.CycleRatherThanScore = true;
	}

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // INPUT ANALYSIS : PROFILE(S)
  ////////////////////////////////////////////////////////////////////////////////////////////////
  int ProfileCount;
	if (OutputType.SearchWith == PRF) {
	  /* allocates memory for the profile structure */
	  prf = (struct Profile *) _mm_malloc(sizeof(struct Profile), 16);
	  if (prf == NULL) {
	      fputs("Unable to allocate memory for the profile structure\n", stderr);
	      exit(1);
	  }
	  
	  /*
	  * Read the profile and output some infos
	  */
	  gettimeofday(&_t0,0);
	  ProfileCount = ReadProfile(ProfileFile, prf, true, CompleteCycleoOnly);
	  gettimeofday(&_t1,0);
	  if (ProfileCount < 0) {
	    fputs("Error found reading profile.\n", stderr);
	    exit(1);
	  }
	  const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
	  
	  if (ProfileCount > 1) {
	    if (OutputVerbose) fprintf(stderr, "Profiles reading took %lf seconds: %i found\n", T, ProfileCount);
	    if (prf->Type == PF_MATRIX) {
	      /*
	       * Performs some tests to see how if profiles are similar
	       */
	      prfs = (struct Profile * *) alloca(ProfileCount*sizeof(struct Profile *));
	      
	      struct Profile * tmpPrf = prf;
	      prfs[0] = prf;
	      MaxProfileLength = prf->Length;
	      size_t count = 1;
	      while (tmpPrf->next) {
					tmpPrf = tmpPrf->next;
					prfs[count++] = tmpPrf;
					if (tmpPrf->Type == PF_PATTERN) {
						fputs("There is some pattern profile within the given list of profiles, coverage cannot be performed.\n",
						stderr);
						exit(1);
					}
					if ((prf->Alphabet_Length != tmpPrf->Alphabet_Length ) || (strncmp(prf->CABC, tmpPrf->CABC, ALPHABET_SIZE+2) != 0)) {
						fputs("There is some inconsistencies within profile's alphabets, coverage cannot be performed.\n",
						stderr);
						exit(1);
					}
					if ((prf->DisjointData.NDIP[0] != tmpPrf->DisjointData.NDIP[0]) || (prf->DisjointData.NDIP[0] != tmpPrf->DisjointData.NDIP[0])) {
						fputs("There is some inconsistencies within profile's disjointness, coverage cannot be performed.\n",
						stderr);
						exit(1);
					}
					MaxProfileLength = (MaxProfileLength < tmpPrf->Length) ? tmpPrf->Length : MaxProfileLength;
	      }
	      
				if (OutputVerbose) {
					fprintf(stderr,"Output base file names is %s\n", OutputFileName);
					fprintf(stderr,"%i profiles of max length %lu and alphabet size of %lu\n",
					        ProfileCount, MaxProfileLength, prf->Alphabet_Length);
					fputs("Alphabet Mapping\n",stderr);
					for (size_t i=0; i<ALPHABET_SIZE; ++i) {
						fprintf(stderr,"Map %c=%2u  ", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) prf->Alphabet_Mapping[i]);
						if ((i+1) % 8 == 0 ) fputs("\n",stderr);
					}
					fputs("\n\n",stderr);

					fprintf(stderr,"Disjoint set: %i to %i\n", prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1]);
				}  
	    }
	    else {
				fputs("There is some pattern profile within the given list of profiles, coverage cannot be performed.\n",
					stderr);
				exit(1);
	    }
	  }
	  else {
	    /* Do we have a single pattern profile to treat */
	    if (prf->Type == PF_PATTERN) {
				fputs("Pattern profile are not accepted in pfcoverage", stderr);
				exit(1);
	    }
	    else if (OutputVerbose) {
	      fprintf(stderr, "Profile reading took %lf seconds: %i found\n", T, ProfileCount);

	      fprintf(stderr,"Profile %s has length %lu and alphabet size of %lu\n",
		    ProfileFile, prf->Length, prf->Alphabet_Length);

	      fputs("Alphabet Mapping\n",stderr);
	      for (size_t i=0; i<ALPHABET_SIZE; ++i) {
					fprintf(stderr,"Map %c=%2u  ", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) prf->Alphabet_Mapping[i]);
					if ((i+1) % 8 == 0 ) fputs("\n",stderr);
	      }
	      fputs("\n",stderr);

	      fprintf(stderr,"Disjoint set: %i to %i\n", prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1]);
	    }
	  }
	  
	  /* Set the Mode and Level provided we are on a matrix profile */
	  if (prf->Type == PF_MATRIX) {
	    prf->CutOff = Cutoff;
	    if (OutputVerbose) fprintf(stderr, "Cutoff set to %i\n", Cutoff);
	  }
	}
#ifdef PRF_CORE_PCRE
	if (OutputType.SearchWith == REGEX) {
		/* Parse the regex argument to get the real regular expression */
		/* Create a copy of the source on the stack */
		const size_t regexSourceSize = 1+strlen(RegexSource);
		char * ctmp = (char *) alloca(regexSourceSize*sizeof(char));
		for (size_t i=0; i<regexSourceSize; ++i) ctmp[i] = RegexSource[i];
		ctmp[regexSourceSize] = '\0';
		
		/* Check format extract '{' '}' border */
		char * Start;
		char * ptr = ctmp;
		char * const End = &ctmp[regexSourceSize-1];
		_Bool ErrorinParsing = false;
		
		while ( *ptr != '{' && ptr <= End) ++ptr;
		if (ptr != End) {
			Start = ++ptr;
			ptr = End;
			while ( *ptr != '}' && ptr <= End) --ptr;
			if (ptr != Start)
				*ptr = '\0';
			else
				ErrorinParsing = true;
		}
		else
			ErrorinParsing = true;
		
		if (! ErrorinParsing) 
			RegexSource = Start;
		else {
			fprintf(stderr, "Unable to retrieve correctly the regular expression or pattern from %s\n", RegexSource);
			exit(1);
		}
		ProfileCount = 1;
		
		if (CplmMode != NONE) {
			fputs("Sequences complementation is not possible on pattern\n", stderr);
			exit(1);
		}
		
		/* Initialize regex structure */
		regex = (struct RegEx*) malloc(sizeof(struct RegEx));
		if (regex == NULL) {
			fputs("Unable to allocate regex structure\n", stderr);
			goto bail;
		}
		if (RegexCycleCount) {
			const size_t HistogramSize = OutputType.CycleRange[1] - OutputType.CycleRange[0] + 1;
			res = InitCycleRegExFromString(regex, RegexSource, HistogramSize, 16);
		}
		else 
			res = InitRegExFromString(regex, RegexSource, false, 16);
		if (res<0) {
			fprintf(stderr, "Regex initialization return error code %i\n", res);
			exit(1);
		}
	}
#endif

	/* Treat complement alphabet */
	if (CplmMode != NONE) {
		if (ComplementAlphabet(prf, CplmMode)!= 0) {
			fputs("Error complementing alphabet\n", stderr);
			res = 1;
			goto bail;
		}
	}

	if (ProfileCount != 1) {
		fputs("Only one profile is allowed at a time, please correct this.\n", stderr);
		goto bail;
	}
  
	////////////////////////////////////////////////////////////////////////////////////////////////
  // THREADING ANALYSIS
  ////////////////////////////////////////////////////////////////////////////////////////////////
  /*
   * Retrieve number of cores
   */
  nCPUs = (nCPUs == 0) ? (size_t) System.nOverallCores : nCPUs;

#ifdef __USE_AFFINITY__
  if (noAffinity) {
    // -----------------------------------------------------------------------------
    //                        ***  NO AFFINITY ***
    // -----------------------------------------------------------------------------
    if (OutputVerbose) fputs("Thread affinity disabled\n", stderr);
    Thread_count[0] = System.nOverallCores;
    Thread_masks[0] = (Affinity_Mask_t*) malloc(System.nOverallCores*sizeof(Affinity_Mask_t));
    for (size_t thread=0; thread<System.nOverallCores; ++thread) {
      CPU_ZERO(&Thread_masks[0][thread]);
      for (int i=0; i<(int) System.nOverallCores; ++i) CPU_SET_S(i, sizeof(Affinity_Mask_t),&Thread_masks[0][thread].data);
    }
  }
  else if (GivenAffinityFile) {
    // -----------------------------------------------------------------------------
    //                     ***  INPUT FILE HOLDING NASKS ***
    // -----------------------------------------------------------------------------
    if (OutputVerbose)
      fprintf(stderr,"Parsing file %s for affinity mask and number of threads\n", AffinityMaskFileName);
    FILE* in = fopen(AffinityMaskFileName, "r");
    if (in == NULL) {
			fprintf(stderr, "Cannot open thread affinity file %s.\n", optarg);
			exit(1);
    }
    size_t lines = 0;
    while (!feof(in)) {
			int num = fread(buffer, sizeof(char), 64, in);
			for (unsigned int i=0; i<num; i++)
		    if (buffer[i] == '\n') lines++;
    }
    rewind(in);
    if (lines != 0) {
			if (lines > System.nOverallCores) lines = System.nOverallCores;
			Thread_masks[0] = (Affinity_Mask_t*) malloc(lines*sizeof(Affinity_Mask_t));
			for (size_t i=0; i<lines; i++) {
			    if (fscanf(in, "%s\n", buffer) == 1) {
			      const size_t tmp_size = strlen(buffer) - 1;
			      CPU_ZERO(&Thread_masks[0][i]);
			      for (int j=tmp_size; j>=0; j--) {
				  if (buffer[j] != '0') CPU_SET_S(j, sizeof(Affinity_Mask_t), &Thread_masks[0][i].data);
			      }
			    }
			}
			Thread_count[0] = lines;
			if (OutputVerbose) fprintf(stderr,"Found %2lu threads affinity masks.",nCPUs);
    }
    else {
			if (OutputVerbose) printf("Cannot understand cpu mask, keep on normally\n");
    }
    fclose(in);
  }
  else if ( split ) {
    // -----------------------------------------------------------------------------
    //                 ***  HALF SSE 2 HALF SSE 4.1 HYPERTHREADING***
    // -----------------------------------------------------------------------------
    Thread_count[0] = getMasks(&System, -1, -1, 1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
    Thread_count[1] = getMasks(&System, -1, -1, 2, &Thread_masks[1]);
    if (Thread_count[1] == 0) {
      fputs("No potential affinity mask found with hyperthreading !!!\n", stderr);
      exit(0);
    }
    if (OutputVerbose)
      fprintf(stderr, "%u threads will use SSE 4.1 and %u SSE 2\n", Thread_count[0], Thread_count[1]);
  }
  else if (noSharedCore) {
    if (OutputVerbose)
      fputs("No sharing of core resources will be used: Intel Hyperthreading or AMD Compute Unit\n", stderr);
    Thread_count[0] = getMasks(&System, -1, -1, 1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
  }
  else {
    // -----------------------------------------------------------------------------
    //                        *** OPERATING SYSTEM CHOICE ***
    // -----------------------------------------------------------------------------
    Thread_count[0] = getMasks(&System, -1, -1, -1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
  }

  {
    register size_t total = (size_t) (Thread_count[0] + Thread_count[1]);
    if (nCPUs > total) nCPUs = total;
  }
 
  threads_attr = (pthread_attr_t*) alloca(nCPUs*sizeof(pthread_attr_t));
  {
    register const Affinity_Mask_t * current = &Thread_masks[0][0];
    for (size_t i=0; i<nCPUs; ++i) {
      pthread_attr_init(&threads_attr[i]);
      if (i == (size_t) Thread_count[0]) current = &Thread_masks[1][0];
      pthread_attr_setaffinity_np(&threads_attr[i], sizeof(Affinity_Mask_t), &(current->data));
      ++current;
    }
  }
#endif
  
  if (OutputVerbose) {
#ifdef USE_PACBIO 
		if (!(PacBioDispatchOptions.nZMW == 1U))
#endif
			fprintf(stderr, "Job dispatched over %lu cores.\n", nCPUs);
	}
  
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // ANALYSIS BASED ON DATABASE
  ////////////////////////////////////////////////////////////////////////////////////////////////
	/* --------------------------------- RANDOM SEQUENCES --------------------------------------- */
#ifdef PRF_INPUT_RANDOM
	if (strncmp(DB, "Random[", 7) == 0) {
		if (sscanf(&DB[6],"[%zu,%u,%u,%u]", &RD.N, &RD.Length, &RD.nCAGSuffix, &RD.seed) != 4) {
			fprintf(stderr, "Error reading random generator data from %s\n", DB);
			res = 1;
		}
		else {
			if (OutputVerbose) {
				fprintf(stderr,"Random Generator activated for %'zu sequences of length %u bearing %u CAG suffices with seed %u\n",
								RD.N, RD.Length, RD.nCAGSuffix, RD.seed);
			}
			if (nCPUs > RD.N) nCPUs = RD.N;
			fprintf(stderr, "Currently removed!\n");
		}
		goto Input_Found;
	}
#endif
	/* -------------------------------- FASTA SEQUENCES STREAM ---------------------------------- */
#ifdef PRF_INPUT_FASTA
	if (( DB[0] == '-' && DB[1] == '\0') || FastaStream) {
		FILE *fd;
		if (FastaStream) {
			fd = fopen(DB, "r");
			if (fd == NULL) {
				fprintf(stderr, "Error opening file %s\n", DB);
				goto bail;
			}
		}
		else {
			fd = stdin;
		}
		res = dispatchStreamFASTA(prf, fd,
#ifdef PRF_CORE_PCRE
                              regex,
#endif
                              &OutputType, nCPUs);
		if (FastaStream) fclose(fd);
		goto Input_Found;
	}
	/* ---------------------------------  FASTA SEQUENCES --------------------------------------- */
	else if ( isFASTA(DB) ) {
		if (OpenFASTAStructure(DB, &FASTA) != 0) {
			fputs("Error opening FASTA file\n", stderr);
			res = 1;
			goto bail;
		}
				
		/* No more threads than sequences in DB file */
		if (nCPUs > FASTA.SequenceCount) nCPUs = FASTA.SequenceCount;
		
		res = dispatchFASTAFile(prf, &FASTA,
#ifdef PRF_CORE_PCRE
                              regex,
#endif
                              &OutputType, nCPUs);
			
		CloseFASTAStructure(&FASTA);
	
		goto Input_Found;
	}
#endif

	/* -------------------------------  PAC BIO HDF5 SEQUENCES ---------------------------------- */
#ifdef PRF_INPUT_HDF5
	if ( isPacBioH5(DB) ) {
		if (RegexCycleCount) {
			fprintf(stderr, "Regex cycle count not be done on Pacbio data yet\n");
			goto bail;
		}
		PBS = PacBioOpen(DB);
		if (PBS == NULL) {
			fprintf(stderr, "Error opening Pac Bio file %s\n", DB);
			goto bail;
		}
		if (!(PBS->Content & BASECALLING)) {
			fputs("Error Pac Bio file is not containing Basecalling data\n", stderr);
			goto bail;
		}
		
		gettimeofday(&_t0,0);
		const int tres = (PacBioDispatchOptions.Selection & PB_DISPATCH_CONSENSUS) ? IndexConsensus(PBS) : IndexRegions(PBS);
		gettimeofday(&_t1,0);
		if (tres != SUCCESS) {
			fprintf(stderr,"PacBio error indexing regions, code %i\n", tres);
			goto bail;
		}
		
		if (OutputVerbose) {
			fprintf(stderr,
							"Base file name : %s\n"
							"Directory      : %s\n"
							"Parts          : %u\n",
							PBS->BaseFileName, PBS->Directory, PBS->nParts);
			for (unsigned int i=0; i<PBS->nParts; i++) {
				fprintf(stderr, "               : %s\n", PBS->PartsFileName[i]);
			}
			
			fprintf(stderr, "ZMW number    : %u\n", PBS->nHoles);
			const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Indexing holes took %lf seconds.\n", T);
		}
		
		if (PopulateHoleStatus(PBS) != SUCCESS) {
			fprintf(stderr, "ZMW status not populated correctly...\n");
			res = 1;
			goto bail;
		}
		
		if (PopulateHoleCoordinates(PBS) != SUCCESS) {
			fprintf(stderr, "ZMW coordinates not populated correctly...\n");
			res = 1;
			goto bail;
		}
		
		if (PacBioDispatchOptions.ZMW == NULL) {
			if (!(PacBioDispatchOptions.Selection & PB_NOT_ONLY_SEQUENCING)) {
				PacBioDispatchOptions.nZMW = getSequencingHoles(PBS, &(PacBioDispatchOptions.ZMW));
				if (PacBioDispatchOptions.nZMW > 0) {
					if (OutputVerbose) fprintf(stderr, "Found %i sequencing holes\n", PacBioDispatchOptions.nZMW);
				}
				else {
					fprintf(stderr, "getSequencingHoles returned %i\n", PacBioDispatchOptions.nZMW);
					goto bail;
				}
			}
			else {
				PacBioDispatchOptions.nZMW = PBS->nHoles;
				PacBioDispatchOptions.ZMW = (unsigned int*) malloc(PacBioDispatchOptions.nZMW*sizeof(unsigned int));
				if (!PacBioDispatchOptions.ZMW) {
					fprintf(stderr, "Unable to allocate memory for list of sequencing holes\n");
					goto bail;
				}
				for(int i=0; i<PacBioDispatchOptions.nZMW; i++) PacBioDispatchOptions.ZMW[i] = i;
			}
		}
		
		if (nCPUs > PacBioDispatchOptions.nZMW) nCPUs = PacBioDispatchOptions.nZMW;
		res = dispatchPacBio(prf, PBS, 
#ifdef PRF_CORE_PCRE
												 regex,
#endif
												 &OutputType,
												 &PacBioDispatchOptions,
												 nCPUs);
		
		
		if ((PacBioDispatchOptions.nZMW != 1) && PacBioDispatchOptions.ZMW) free(PacBioDispatchOptions.ZMW);
		if (PacBioClose(PBS) != SUCCESS) {
			fprintf(stderr, "Error closing Pac Bio...\n");
		}
		
		goto Input_Found;
	}
#endif
	/* -------------------------------  PAC BIO BAM SEQUENCES ----------------------------------- */
#ifdef PRF_INPUT_PBBAM
	if (isPacBioBAM(DB)) {
		PacBioBAM_t * const PBBAM = OpenPacBioBAM(DB);
		if (PBBAM) {
			res = dispatchPacBioBAM(prf, PBBAM, 
#ifdef PRF_CORE_PCRE
														regex,
#endif
														&OutputType,
														&PacBioDispatchOptions,
														nCPUs);
			ClosePacBioBAM(PBBAM);
		}
		else {
			fprintf(stderr, "Error in OpenPacBioBAM\n");
			res = 1;
		}
		goto Input_Found;
	}
#endif
	/* ------------------------------------  UNDEFINED ------------------------------------------ */
	
	fprintf(stderr, "Database file %s could not be determined\n.", DB);
	res = 1;
	Input_Found:;

	////////////////////////////////////////////////////////////////////////////////////////////////
  // CLEANLY CLOSE
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  /* Free Memory */
bail:;
  if (prf) FreeProfile(prf, true);
#ifdef PRF_CORE_PCRE
	if (regex) free(regex);
#endif
	
#ifdef __USE_AFFINITY__
  if (Thread_masks[0]) free(Thread_masks[0]);
  if (Thread_masks[1]) free(Thread_masks[1]);
#endif
	
	freeSystemInfo(&System);

  return res;
}

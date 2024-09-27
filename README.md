# Repeat Detector
## About
Repeat Detector (RD) is a deterministic profile weighting algorithm for counting repeats in targeted sequencing data. Based on the pfsearch algorithm, it is an alignment free method that uses unaligned FASTA files, along with a preset motif file. 

**New executable GUI for Windows and macOS**

Instructions and the application/executable are located under the macOS and Windows folders, respectively, to avoid dependency and installation issues.

To use on Windows: https://github.com/DionLab/RepeatDetector/tree/main/Windows
To use on macOS: https://github.com/DionLab/RepeatDetector/tree/main/macOS

For Linux/HPC, use the Singularity image repeat_detector.sif, available here: https://zenodo.org/records/13847199/files/repeat_detector.sif?download=1.

Currently, only the restrictive profile is available. The permissive profile, along with other features, will be added in the next update!

## Installation
Repeat Detector is available as a built executable in the 'RepeatDetector-1.0.15eb445-Release-Linux-x86_64' folder.

To build from scratch, the source code of Repeat Detector and its dependency libprf (DionLab/libprf) can be installed using:
```
git clone https://github.com/DionLab/libprf.git
git clone https://github.com/DionLab/RepeatDetector.git
```
WARNING: libprf MUST be installed in the same directory prior to installing RepeatDetector

Once the repository has been cloned, enter the directory (i.e. RepeatDetector) and perform the following commands:

```
mkdir build
cd build
cmake ../
make 
make install
```
HTSLib must also be installed for RepeatDetector.

## General Usage
Repeat Detector requires an input fasta file and a selected profile. When run with no output format selected, RD will output a table of reads, with associated longest repeat tract to stdout. 

Profiles either allow an interruption in the repeat tract (permissive/Annex2) or do not allow an interruption (restrictive/Annex10).

```
RepeatDetector --prf <profile> <fasta-file> <options>

### Options
  Sequence
    --sequence-clip [b:e]          : constrain alignement to start after b (nucleotide position start) and to end before e (nucleotide position end)
    --score-range [a:b]            : contrain score within given range (int)
    --cycle-range [a:b]            : contrain repeat cycle within given range (int)
    --not                          : inverse above criteria
    --with-revcomp                 : test also reverse complement
    --bestof                       : when used with --with-revcomp and --optimal,
                                     output only the best of both orientation
                                     this option is for text output only, in
                                     histogram mode when reverse complement is
                                     required, best of is implicit.

  Profile
    --cutoff                  [-c] : score cutoff value
    --optimal                 [-a] : output only optimal alignment per sequence

  Regex
   --cycle-count                   : account cycle of regex pattern till
                                     histogram size (--histogram-bins)

  Database
   FASTA
    Fasta format allows to be piped with the '-' symbol as database file name
     --create-index-database  [-I] : output indices to given file
     --use-index-database     [-i] : use indices stored in given file
     --stream-fasta                : stream the file rather than indexing

   PacBio HDF5 or BAM
     --subreads                    : work on subreads (default)
     --consensus                   : work on consensus
     --polymerase                  : work on entire polymerase data
     --invalid-zmw                 : work on invalid zmw only
     --best-subread                : keep only best subread of zmw
     --zmw-number i                : limit to zmw i
     --not-only-sequencing         : compute on any zmw type
     --min-filter-score i          : minimal High Quality region score i to keep data
     --max-filter-score i          : maximum High Quality region score i to keep data
     --keep-invalid-zmw            : keep zmw removed by PacBio
     --discard-filter              : does not filter or trim according to HQ region
     --test-filter                 : output regions based upon the above filtering

 Optimizations
   --nthreads                 [-t] : max number of threads to use
   --no-affinity                   : disable CPU affinity file
   --thread-affinity               : file containing thread mask,
                                     one row for one thread
   --no-shared-core                : Prevent core resource sharing
   --split                         : if both SSE 2 & 4.1 are available,
                                     split half-half using linked resources

  Output
    --output-name                  : specify output file name
    --output type             [-o] : output in format type
    --output-width            [-X] : text output width
   type is:
    - dat                          : ASCII matrix dump
    - tab                          : tab separated matrix dump
    - histogram                    : generate an histogram of number of cycles
      --score-range                :
      --cycle-range                :
      --use-cycle                  : in case both score and cycle range are
                                   : provided, use cycle rather than score
    - density                      : generate a 2D histogram, score vs cycles
      --score-range and --cycle-range are mandatory
    - xPSA
    - TSV
    - OneLine
    - FASTA                        : FASTA file format
      --source                     : output entire source sequence
      --before
      --after
      --in-between
    - pdf                          : pdf outout
      --whole-sequence        [-w] : does not limit output to alignment
      --with-path                  : does output paths
      --with-graphs                : does output side graphs
      --region [a,b][c,d]     [-r] : specify region to output
  Other
    --verbose                 [-V] : verbose on stderr
    --help                    [-h] : output command help
 ```
## Repeat Detector Paper
Alysha S Taylor, Dinis Barros, Nastassia Gobet, Thierry Schuepbach, Branduff McAllister, Lorene Aeschbach, Emma L Randall, Evgeniya Trofimenko, Eleanor R Heuchan, Paula Barszcz, Marc Ciosi, Joanne Morgan, Nathaniel J Hafford-Tear, Alice E Davidson, Thomas H Massey, Darren G Monckton, Lesley Jones, REGISTRY Investigators of the European Huntington’s disease network, Ioannis Xenarios, Vincent Dion. Repeat Detector: versatile sizing of expanded tandem repeats and identification of interrupted alleles from targeted DNA sequencing, NAR Genomics and Bioinformatics, Volume 4, Issue 4, December 2022, lqac089, https://doi.org/10.1093/nargab/lqac089

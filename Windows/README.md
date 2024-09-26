# Windows Instructions

1. Download [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/), install it, and ensure it is running.
2. Download the Repeat Detector Docker image from [Zenodo](https://zenodo.org/api/records/12773712/draft/files/repeat-detector.tar) (link to be modified).
3. Copy the GUI files and run them locally from your computer rather than from networked folders.
4. Locate `RDProgramRunner.exe` and double-click it to run.
5. If Windows shows a security warning, click **More info**, then click **Run anyway**.
6. When the app opens, browse and locate the `repeat-detector.tar` file.
7. Next, browse and locate the FASTA file. If you have FASTQ/FASTQ.gz files, use `seqkit fq2fa` to convert them to FASTA format.
8. Once everything is loaded, click **Run RD Program**.
9. It may take a minute or two for the Docker image to load.
10. Histograms will be stored in the same location as the source FASTA file.


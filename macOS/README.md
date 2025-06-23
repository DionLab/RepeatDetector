# macOS Instructions

1. Download Docker Desktop for macOS from [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/), install it, and ensure it is running.
2. Download the Repeat Detector Docker image from [Zenodo](https://zenodo.org/records/13847199/files/repeat-detector.tar?download=1) (link to be modified).
3. Copy the GUI files and run them locally from your computer.
4. Download RDProgramRunner from here: https://zenodo.org/records/15721525. Locate, unzip and double-click 'RDProgramRunner.app' to run.
5. If macOS shows a security warning, go to **System Preferences** > **Security & Privacy**, and allow the app to run.
6. When the app opens, browse and locate the `repeat-detector.tar` file.
7. Next, browse and locate the FASTA file. If you have FASTQ/FASTQ.gz files, use `seqkit fq2fa` to convert them to FASTA format.
8. Once everything is loaded, click **Run RD Program**.
9. It may take a minute or two for the Docker image to load.
10. Histograms will be stored in the same location as the source FASTA file.

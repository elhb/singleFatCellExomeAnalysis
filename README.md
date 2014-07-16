#singleFatCellExomeAnalysis
### Intro
Scripts used for the analysis of single cells and recipient and donor whole blood samples in the study <STUDY>

###The scripts
There are 5 scripts in the analysisScripts folder, the five scripts are the following:   

1. startAnalysis
1. filterVariants
1. removeEmptyReads.py
1. wgaAdapterTrimmer.py  
1. mappingStatsForExcel.py  

###Reproducing the analysis
####Note on Uppmax and Slurm
All scripts were run at Uppmax, Uppsala University's resource of high-performance computers in Uppsala, Sweden (http://www.uppmax.uu.se/). The Uppmax resource uses the SLURM workload manager (http://slurm.schedmd.com/) and the scripts are therefore configure to use this system. it also loads binaries such as bowtie2, fastqc etc into the PATH through "module load" commands.  If these type of systems is not present on you machine you should be able to run the automatically generated sbatch scripts using bash (possibly after some manual editing).

####Dependencies
To run the analysis some additional software is requiered:

1. TrimBWAstyle.pl (from http://wiki.bioinformatics.ucdavis.edu/index.php/TrimBWAstyle.pl)
1. samtools (http://samtools.sourceforge.net/)
1. picard tools (http://picard.sourceforge.net/)
1. bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
1. cutadapt (https://code.google.com/p/cutadapt/)
1. fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
1. GATK (https://www.broadinstitute.org/gatk/)

####Hardcoded variables
Some paths and filenames etc are hardcoded, if you want to reproduce the analysis you need to change the paths in the script (or rename any files downloaded from the SRA to match the scripts as well as changing the locations of executables on your machine).

####Do not hesitate to contact me if you need assistance in running the scripts.

a pipe
==========================

$ analysisScripts/startAnalysis

usage:  [-h] -s <filename.tsv> [-p <path>] [-d <filename>] -id STRING
        [--submitDummies] [--submitCollection] [-A <projectId>]
        [-x <ref.fa location>] [--GATK <path>] [--Picard <path>]
        [--Bundle <path>]

optional arguments:

  -h, --help            show this help message and exit

  -s <filename.tsv>, --samples <filename.tsv>
                        Tab delimited file with sample information.

  -p <path>, --resultsPath <path>
                        Path tu where all results should be stoored (default ./results).

  -d <filename>, --1000G_bams <filename>
                        List of 1000G bam files to use for the recalibration and realignment steps.

  -id STRING, --collectionId STRING
                        Name for the sample collection.

  --submitDummies       Flag for subbmitting scripts for gvcfs creation from 1000Genomes bam files.

  --submitCollection    Flag for subbmitting scripts for genotype gvcfs and vqsr.

  -A <projectId>, --project <projectId>
                        ProjectId (default b2011011).

  -x <ref.fa location>, --reference <ref.fa location>
                        reference location (default "~/singleFatCellExomeAnalysis/references/GATKbundle/human_g1k_v37.fasta").

  --GATK <path>         GATK location (default ~/singleFatCellExomeAnalysis/bin/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar).

  --Picard <path>       Picard location (default ~/bin/picard-tools-1.114).

  --Bundle <path>       GATKbundle location (default ~/singleFatCellExomeAnalysis/references/GATKbundle/ ).

example commandline:
$ analysisScripts/startAnalysis --samples samples_list_example.tsv  --collectionId testSampleCollection --submitCollection -p resultsAreStoredHere -A b2011001 -x referencefile.fa -d 1000_genomes_bams_list_example.tsv  --submitDummies

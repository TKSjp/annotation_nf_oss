# annotation_nf_oss

Nextflow annotation pipeline for Entamoeba genome.

infernal-tblout2gff.pl needs to be copied from https://github.com/nawrockie/jiffy-infernal-hmmer-scripts/

## environment

- This script is written for super computer in SGE system, especially shirokane (https://gc.hgc.jp)
- You need nextflow and conda environment written in config files.

## setup of databases

- before use, set up databases
    - Rfam.raw.cm
        - wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
        - gunzip Rfam.cm.gz
        - wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
        - cmpress Rfam.cm
    - uniref90.dmnd
        - wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
        - gunzip uniref90.fasta.gz
        - diamond makedb --in uniref90.fasta -d uniref90.dmnd
    - augustus/config
        - git clone https://github.com/Gaius-Augustus/Augustus.git -b v3.5.0 (in somewhere)
        - set --AUGUSTUS_CONFIG_PATH to augustus/config

## setup of conda environment

see config and prepare all environments.

## settings

- before run, modify annotation.config
    - mkdir log input
    - git clone https://github.com/TKSjp/annotation_nf script
    - run repeatmasker first (take 8hours) and use masked fasta for annotation
    - modify config
        - set correct paths for the scripts, conda environments, sequence files, and containers (.nf and .config files)
        - rename_path = '!{set full path to the script} rename_fasta_by_length.py'
        - rename_options = '--width 2 --prefix Chr' *this will Chr01. It is depends on contig number*
        - trnascan_modes = ['E', 'O'] // E: eukaryote, B: bacteria, A: archaea, O: other organellar, G: general
        - trnascan_gff = '!{set full path to the script} tRNAScan_SE_to_gff3.py'
        - infernal_cm = '!{set full path to the rfam/Rfam.raw.cm}'
        - infernal_tblout2gff = '!{set full path to the script} infernal-tblout2gff.pl'
        - diamond_db = '!{set full path to the uniref90/uniref90.dmnd}'
        - diamond_options = '--sensitive --evalue 1e-5 --max-target-seqs 10 --query-cover 50 --subject-cover 50 --outfmt 6 qseqid sseqid pident length  mismatch gapopen qstart qend sstart send evalue bitscore stitle'
        - adddiamond2gffproduct = '!{set full path to the script} adddiamond2gffproduct.py'
        - addhypothetical2emptycds = '!{set full path to the script} addhypothetical2emptycds.py'
        - augustus_config_path = '!{set full path to the augustus/config}'
        - braker_protseqfasta = '!{set full path to the braker/protseq.faa depends on species}'
        - braker_species = 'hoge' // if new run you need to change. if rerun keep it.

- rna-seq pipeline config
    - put fastq like (Eh_1.fastq, Eh_2.fastq), Eh.fna, Eh.gff in input directory (Eh is recognized as sample name)
    - rdna   = '/home/hoge/plasmid.fna' // for removing ribosomal RNA reads. they are in plasmid in case of Entamoeba.
    - fnadir = '/home/hoge/input/' // for input fasta dir path (last "/" is mandatory)
    - gffdir = '/home/hoge/input/' // for input gff dir path (last "/" is mandatory)

## run

here, you need log dir, input dir, and script dir which has each pipeline script.

- qsub -N repeatmasker script/repeatmasker.sh
- qsub -N annotation script/annotation.sh
- qsub -N rnaseq script/rnaseq_pe.sh

## output of the pipeline

results are in result directory

- braker
- diamond
- infernal
- rename
- summary *this is the final output*
    - id.{faa,ffn,faa.stat,ffn.stat,fna,gff}
    - version.txt
- trnascan/{E,O}
    - E has iMet and Met while O has only Met (mixed).
    - O has more tRNAs than E (in case of Entamoeba genome)

## output of repeatmasker

## output of rna-seq pipeline

- fastp
- hisat2
- removerrna
- remapping
- remapping2 *use bam in here to look up in IGV*
- featurecount
    - counts_join.tsv *this is the final output*
    - counts.txt
    - counts.txt.summary
    - version.txt

## for submission to NCBI

- You need modifications below:
    - remove tRNA-Sup from gff
    - manual check of ncRNA annotation by infernal (maybe several incorrect ncRNAs need to be removed)
    - "homolog", and "homologue" must be convert to "-like protein"
    - some orphan genes should be removed using removeorphangene.py

## contact

If you have any troubles, please feel free to open an issue.
Tetsuro Kawano-Sugaya (sugaya@tetsu.ro)
https://tetsu.ro
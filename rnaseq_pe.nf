#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process fastp {

    publishDir "result/fastp", mode: 'copy', overwrite: true
    conda '/home/user/miniconda3/envs/fastp'

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_qc_1.fastq.gz"), path("${id}_qc_2.fastq.gz"), emit: trimmed_reads
    path "${id}_report.html"

    shell:
    '''
    read1=!{reads[0]}
    read2=!{reads[1]}
    fastp -w !{task.cpus} -f 10 -F 10 -q 20 -i ${read1} -I ${read2} -o !{id}_qc_1.fastq.gz -O !{id}_qc_2.fastq.gz --html !{id}_report.html
    '''
}

process hisat2 {

    publishDir "result/hisat2/${id}", mode: 'copy', overwrite: true
    conda '/home/user/miniconda3/envs/hisat2'

    input:
    val(reference)
    tuple val(id), path(qc1), path(qc2)

    output:
    tuple val(id), path("${id}.bam"), emit: bam
    path "${id}.bam.bai"
    path "version.txt"

    shell:
    '''
    ref="!{reference}!{id}.fna"
    hisat2-build -p !{task.cpus} ${ref} INDEX
    hisat2 -x INDEX -1 !{qc1} -2 !{qc2} -p !{task.cpus} --max-intronlen 2000 | samtools sort -@ !{task.cpus} -o !{id}.bam && samtools index -@ !{task.cpus} !{id}.bam
    hisat2 --version > version.txt
    '''
}

process removerrna {

    publishDir "result/removerrna/${id}", mode: 'copy', overwrite: true
    conda '/home/user/miniconda3/envs/removerrna'

    input:
    val(reference)
    val(rdna)
    tuple val(id), path(bam)

    output:
    tuple val(id), path("${id}.wordna.bam"), emit: bam    
    path "${id}.tsv"
    path "${id}.bed"

    shell:
    '''
    ref="!{reference}!{id}.fna"
    makeblastdb -in ${ref} -dbtype nucl -out db -parse_seqids
    blastn -db db -query !{rdna} -out !{id}.tsv -outfmt 6
    cut -f2,9,10 !{id}.tsv | awk '{if($2 > $3){print $1 "\t" $3-1 "\t" $2}else{print $1 "\t" $2-1 "\t" $3}}' > !{id}.bed
    intersectBed -abam !{bam} -b !{id}.bed -v > !{id}.wordna.bam
    '''
}

process remapping {

    publishDir "result/remapping/${id}", mode: 'copy', overwrite: true
    conda '/home/user/miniconda3/envs/remapping'

    input:
    tuple val(id), path(bam)

    output:
    path "${id}.tsv"
    path "version.txt"

    shell:
    '''
    samtools view -@ !{task.cpus} -O bam -F4 -F256 -F2048 -o !{id}.mapped.bam   !{bam}
    samtools view -@ !{task.cpus} -O bam -f4 -F256 -F2048 -o !{id}.unmapped.bam !{bam}
    map=`  samtools view -@ !{task.cpus} -c !{id}.mapped.bam`
    unmap=`samtools view -@ !{task.cpus} -c !{id}.unmapped.bam`
    ratio=`echo $map $unmap | awk '{OFMT="%.6f"} {print $1/($1 + $2)}'`
    echo -e "!{id}\t${map}\t${unmap}\t${ratio}" > !{id}.tsv
    sed -i "1i count_genome_id\tcount_mapped\tcount_unmapped\tmappedratio" !{id}.tsv
    samtools --version | grep "samtools" > version.txt
    '''
}

process remapping2 {

    publishDir "result/remapping2/${id}", mode: 'copy', overwrite: true
    conda '/home/user/miniconda3/envs/remapping'

    input:
    tuple val(id), path(bam)

    output:
    path "${id}.tsv"
    path "version.txt"
    tuple val(id), path("${id}.mapped.sort.bam"), emit: bam    
    path("${id}.mapped.sort.bam.bai")

    shell:
    '''
    samtools view -@ !{task.cpus} -O bam -F4 -F256 -F2048 -o !{id}.mapped.bam   !{bam}
    samtools sort -@ !{task.cpus} -o !{id}.mapped.sort.bam !{id}.mapped.bam
    samtools index !{id}.mapped.sort.bam
    samtools view -@ !{task.cpus} -O bam -f4 -F256 -F2048 -o !{id}.unmapped.bam !{bam}
    map=`  samtools view -@ !{task.cpus} -c !{id}.mapped.sort.bam`
    unmap=`samtools view -@ !{task.cpus} -c !{id}.unmapped.bam`
    ratio=`echo $map $unmap | awk '{OFMT="%.6f"} {print $1/($1 + $2)}'`
    echo -e "!{id}\t${map}\t${unmap}\t${ratio}" > !{id}.tsv
    sed -i "1i count_genome_id\tcount_mapped\tcount_unmapped\tmappedratio" !{id}.tsv
    samtools --version | grep "samtools" > version.txt
    '''
}

process featurecount {

    publishDir "result/featurecount/${id}", mode: 'copy', overwrite: true
    conda '/home/user/miniconda3/envs/featurecount'

    input:
    val(gffdir)
    tuple val(id), path(bam)

    output:
    path "counts_join.tsv"
    path "counts.txt"
    path "counts.txt.summary"
    path "version.txt"

    shell:
    '''
    gff="!{gffdir}!{id}.gff"
    featureCounts -T !{task.cpus} -M -O -p -B -C --countReadPairs --fraction -t exon -g Parent -s 0 -a ${gff} -o counts.txt !{bam}
    featureCounts -v &> version.txt
    python !{params.join_product} counts.txt !{gffdir}!{id}.gff > counts_join.tsv 
    '''
}

process tpmcalc {

    publishDir "result/tpmcalc/${id}", mode: 'copy', overwrite: true
    conda '/home/user/miniconda3/envs/tpmcalc'

    input:
    val(gffdir)
    tuple val(id), path(bam)

    output:
    path "${id}.mapped.sort_genes.ent"
    path "${id}.mapped.sort_genes.out"
    path "${id}.mapped.sort_genes.uni"    
    path "${id}.agat.log"    
    path "${id}.gtf"
    path "version.txt"
    path "log.txt"

    shell:
    '''
    gff="!{gffdir}!{id}.gff"
    agat_convert_sp_gff2gtf.pl --gff ${gff} -o !{id}.gtf
    TPMCalculator -g !{id}.gtf -b !{bam} > output.txt 2> log.txt
    TPMCalculator -version > version.txt
    '''
}

workflow {
    fastp(channel.fromFilePairs('input/*_{1,2}.fastq.gz'))
    hisat2(params.fnadir, fastp.out.trimmed_reads)
    removerrna(params.fnadir, params.rdna, hisat2.out.bam)
    remapping(hisat2.out.bam)
    remapping2(removerrna.out.bam)
    //tpmcalc(params.gffdir, remapping2.out.bam)
    featurecount(params.gffdir, remapping2.out.bam)
}


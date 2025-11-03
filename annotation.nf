#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process rename {
    publishDir "result/rename/${id}", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fna)

    output:
    tuple val(id), path("${id}_renamed.fna"), emit:fna

    shell:
    '''
    python !{params.rename_path} --in !{fna} --out !{id}_renamed.fna !{params.rename_options}
    '''
}

process trnascan {
    publishDir "result/trnascan/${id}/${mode}", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fna), val(mode)

    output:
    path "${id}_${mode}.bed"
    path "${id}_${mode}.fasta"
    path "${id}_${mode}.stats"
    tuple val(id), val(mode), path("${id}_${mode}_woheader.gff"), emit:gff
    tuple val(id), val(mode), path("version_tRNAscan-SE.txt"), emit:version

    shell:
    '''
    tRNAscan-SE --thread !{task.cpus} -!{mode} --bed !{id}_!{mode}.bed --fasta !{id}_!{mode}.fasta --stats !{id}_!{mode}.stats !{fna}
    python !{params.trnascan_gff} !{id}_!{mode}.bed > !{id}_!{mode}.gff
    grep -v "^#" !{id}_!{mode}.gff > !{id}_!{mode}_woheader.gff
    tRNAscan-SE /dev/null 2>&1 | sed 1d | head -n1 > version_tRNAscan-SE.txt
    '''
}

process infernal {
    publishDir "result/infernal/${id}", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fna)

    output:
    path "${id}.tblout"
    path "${id}.out"
    path "${id}_infernal.gff"
    path "${id}_infernal_wotrna.gff", emit:gff
    path "version_infernal.txt", emit:version
    
    shell:
    '''
    cmscan --cut_ga --tblout !{id}.tblout --fmt 2 --cpu !{task.cpus}  -o !{id}.out !{params.infernal_cm} !{fna}
    perl !{params.infernal_tblout2gff} --cmscan --fmt2 !{id}.tblout > !{id}_infernal.gff
    grep -v "tRNA" !{id}_infernal.gff > !{id}_infernal_wotrna.gff
    cmscan --devhelp | sed 1d | head -n1 > version_infernal.txt
    '''
}

process braker {
    publishDir "result/braker/${id}", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fna)

    output:
    path "intermediate"
    path "${id}.gff", emit:gff
    tuple val(id), path("${id}.faa"), emit:faa
    path "${id}.ffn", emit:ffn
    path "version_braker.txt", emit:version

    shell:
    '''
    braker.pl --genome=!{fna} --prot_seq=!{params.braker_protseqfasta} --workingdir=intermediate --threads !{task.cpus} --gff3 --species !{params.braker_species} --useexisting --AUGUSTUS_CONFIG_PATH !{params.augustus_config_path}
    cp intermediate/braker.gff3 !{id}.gff
    cp intermediate/braker.aa !{id}.faa
    cp intermediate/braker.codingseq !{id}.ffn
    braker.pl --version > version_braker.txt
    '''
}

process diamond {
    publishDir "result/diamond/${id}", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(faa)

    output:
    path "${id}.tsv"
    path "${id}_best.tsv", emit:tsv
    path "version_diamond.txt", emit:version

    shell:
    '''
    diamond blastp -d !{params.diamond_db} -q !{faa} -o !{id}.tsv --threads !{task.cpus} !{params.diamond_options}
    sort -k1,1 -k12,12nr !{id}.tsv | awk '!seen[$1]++' > !{id}_best.tsv
    diamond --version > version_diamond.txt
    '''
}

process summary {
    publishDir "result/summary/${id}", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fna)
    path(trnascan_gff)
    path(trnascan_version)
    path(infernal_gff)
    path(infernal_version)
    path(braker_gff)
    tuple val(id2), path(braker_faa)
    path(braker_ffn)
    path(braker_version)
    path(diamond_tsv)
    path(diamond_version)

    output:
    path "${id}.fna"
    path "${id}.faa"
    path "${id}.ffn"
    path "${id}.faa.stat"
    path "${id}.ffn.stat"
    path "${id}.gff"
    path "version.txt"

    shell:
    '''
    # gff integration
    cp !{fna} !{id}.fna
    mv !{braker_gff} braker.gff
    seqkit stat !{fna} > !{id}.fna.stat
    python !{params.adddiamond2gffproduct} braker.gff !{diamond_tsv} > braker_product.gff
    python !{params.addhypothetical2emptycds} braker_product.gff | grep -v "^#" > braker_product_hypothetical.gff
    cat braker_product_hypothetical.gff !{infernal_gff} !{trnascan_gff} > braker_infernal_trnascan.gff
    python !{params.fix_cmscan_gff} braker_infernal_trnascan.gff > braker_infernal_trnascan_cmscanfixed.gff
    sed "s/;anticodon=.*//g" braker_infernal_trnascan_cmscanfixed.gff | grep -v "UndetNNN" > trna_fixed.gff

    # extract expected pattern
    mv !{id}.faa !{id}_ori.faa
    seqkit stat !{id}_ori.faa > !{id}_ori.faa.stat
    seqkit grep -srp "^M[ACDEFGHIKLMNPQRSTVWY]+\\*$" !{id}_ori.faa > !{id}.faa
    seqkit stat !{id}.faa > !{id}.faa.stat
    grep "^>" !{id}.faa | sed "s/>//g" > takeid.tsv
    mv !{id}.ffn !{id}_ori.ffn
    seqkit stat !{id}_ori.ffn > !{id}_ori.ffn.stat
    seqkit grep -n -f takeid.tsv !{id}_ori.ffn > !{id}.ffn
    seqkit stat !{id}.ffn > !{id}.ffn.stat

    # extract excluded pattern and remove from gff
    seqkit grep -vsrp "^M[ACDEFGHIKLMNPQRSTVWY]+\\*$" !{id}_ori.faa > vsrp.faa
    grep "^>" vsrp.faa | sed "s/>//g" > excludeid.tsv
    grep -Fv -f excludeid.tsv trna_fixed.gff > !{id}.gff

    # version summary
    seqkit version > version_seqkit.txt
    cat version_seqkit.txt !{trnascan_version} !{infernal_version} !{braker_version} !{diamond_version} > version.txt
    '''
}

workflow {
    fna = Channel.fromFilePairs("input/*.fna", size: 1)    
    rename(fna)
    trnascan(rename.out.combine(params.trnascan_modes))
    trna_gff_O = trnascan.out.gff.filter { id, mode, gff -> mode=='O' }.map { id, mode, gff -> gff }
    trna_ver_O = trnascan.out.version.filter { id, mode, ver -> mode=='O' }.map { id, mode, ver -> ver }
    infernal(rename.out)
    // 乱数文字列チャネルを定義
    //rndCh = Channel.value(
    //    UUID.randomUUID().toString().substring(0,8).replaceAll('-', '')
    //)
    braker(rename.out)
    diamond(braker.out.faa)
    summary(rename.out, trna_gff_O, trna_ver_O, infernal.out.gff, infernal.out.version,braker.out.gff, braker.out.faa, braker.out.ffn, braker.out.version, diamond.out.tsv, diamond.out.version)
}

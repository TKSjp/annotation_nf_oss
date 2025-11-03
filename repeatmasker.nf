nextflow.enable.dsl=2

process repeatmasking {
    publishDir "result/repeatmasker/${id}", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fna)

    output:
    path "RM"
    path "${id}_repeatmasker"
    path "${id}.fna.out"
    path "${id}.fna.out.gff"
    path "${id}.hmask.fna"
    path "${id}.smask.fna"

    shell:
    '''
    BuildDatabase -name !{id} !{fna}
    RepeatModeler -database !{id} -threads !{task.cpus} -LTRStruct
    modeler_dir=$(ls -d RM_* | tail -n 1) 
    mv ${modeler_dir} RM
    cat RM/consensi.fa !{params.customRepeats} > RM/consensus_plus.fa
    RepeatMasker -pa 2 -dir !{id}_repeatmasker -lib RM/consensus_plus.fa !{fna}
    ProcessRepeats -maskSource !{fna} -xsmall -gff !{id}_repeatmasker/!{id}.fna.cat.gz
    mv !{id}.fna.masked !{id}.smask.fna
    mv !{id}_repeatmasker/!{id}.fna.masked  !{id}.hmask.fna
    mv !{id}_repeatmasker/!{id}.fna.out     !{id}.fna.out
    mv !{id}_repeatmasker/!{id}.fna.out.gff !{id}.fna.out.gff
    '''
}

workflow {
    repeatmasking(channel.fromFilePairs('input/*.fna', size: 1))
}

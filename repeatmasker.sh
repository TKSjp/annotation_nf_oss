#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -pe def_slot 4
#$ -l s_vmem=12G,mem_req=12G
#$ -q '!mjobs_rerun.q'

# usage: qsub -N repeatmasker script/repeatmasker.sh

source ~/.bashrc.intr
module use /usr/local/package/modulefiles
module load java/11
module load apptainer

unset JAVA_TOOL_OPTIONS
export JAVA_TOOL_OPTIONS="-XX:+UseG1GC -Xmx10g -Xms10g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=1"
export NXF_VER=23.10.1
export NXF_EXECUTOR=sge
export NXF_OPTS="-Xms10g -Xmx10g"

appname="repeatmasker"

nextflow -log log/${appname}/$(date +%Y%m%d_%H%M%S).txt run script/${appname}.nf -c script/${appname}.config -profile standard -with-report report/${appname}/$(date +%Y%m%d_%H%M%S).html -with-trace trace/${appname}/$(date +%Y%m%d_%H%M%S).txt -work-dir work/${appname}/$(date +%Y%m%d_%H%M%S) -resume

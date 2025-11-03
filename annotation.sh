#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -o log
#$ -e log
#$ -pe def_slot 4
#$ -l s_vmem=8G,mem_req=8G
#$ -q '!mjobs_rerun.q'

# usage: qsub -N annotation script/annotation.sh

source ~/.bashrc.intr
module load java/11
unset JAVA_TOOL_OPTIONS
export JAVA_TOOL_OPTIONS="-XX:+UseG1GC -Xmx8g -Xms8g -XX:ParallelGCThreads=3 -XX:ConcGCThreads=1"
export NXF_VER=23.10.1
export NXF_EXECUTOR=sge
export NXF_OPTS="-Xms8g -Xmx8g"

nextflow -log log/annotation/$(date +%Y%m%d_%H%M%S).txt run script/annotation.nf -c script/annotation.config -profile standard -with-report report/annotation/$(date +%Y%m%d_%H%M%S).html -with-trace trace/annotation/$(date +%Y%m%d_%H%M%S).txt -work-dir work/annotation/$(date +%Y%m%d_%H%M%S) -resume

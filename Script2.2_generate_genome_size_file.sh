#!/bin/bash
#SBATCH --job-name=Script0_merge_fastqs


echo -e "\n######## Script 0 to find for ######## \n"

# rsync -auh --info=progress2 /rds/prj/id_hill_sims_wellcda/genomes/mm10/Sequence/WholeGenomeFasta/genome.fa.fai /scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.genome.fa.fai
cut -f1,2 mm10.genome.fa.fai > mm10.genome.size

echo "##### Finished merging fastqs #####"

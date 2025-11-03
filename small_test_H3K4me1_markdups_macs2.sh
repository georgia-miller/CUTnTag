#!/bin/bash
#SBATCH --job-name=Script1_H3K4me1
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/smalltest_H3K4me1_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/smalltest_H3K4me1_%j.log

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# load the Anaconda module
module load anaconda3/2022.10-gcc-13.2.0

# source conda so it works in non-interactive shells
#source $(dirname $(which conda))/../etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

# set java options so the picard arguments run with these
export JAVA_OPTS="-Xmx30g -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8"

cd /scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling/H3K4me1/H3K4me1_WT_r1/

conda activate CUTnTag_alignment_env
# remove duplicates is set to false so they are only marked


picard MarkDuplicates \
			-I H3K4me1_WT_r1.picardchrsorted.rg.picard.bam \
			-O H3K4me1_WT_r1.picardchrsorted.rg.picard.marked.bam \
			--TMP_DIR ./picard_temp \
			--VALIDATION_STRINGENCY "LENIENT" \
			-METRICS_FILE H3K4me1_WT_r1_PicardMarkDuplicates.txt \
			--REMOVE_DUPLICATES false \
			-Xmx8g

echo -e "\n ######## [`timestamp`] Finished mark duplication for H3K4me1_WT_r1 ######## \n"

intersectBed -v -a H3K4me1_WT_r1.picardchrsorted.rg.picard.marked.bam \
			/scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.blacklisted_and_chrM.sorted.bed > \
			H3K4me1_WT_r1.marked.cleaned.bam

samtools sort -o H3K4me1_WT_r1.marked.cleaned.chrsorted.bam \
			-T H3K4me1_WT_r1.marked.cleaned.chrsorted \
			-@ 16 \
			H3K4me1_WT_r1.marked.cleaned.bam

samtools index H3K4me1_WT_r1.marked.cleaned.chrsorted.bam

conda deactivate

conda activate CUTnTag_macs2_env_2

echo -e "\n ######## [`timestamp`] Starting peak calling for H3K4me1_WT_r1 ######## \n"

macs2 callpeak -t H3K4me1_WT_r1.marked.cleaned.chrsorted.bam \
			-f BAMPE -n H3K4mr1e1_WT_r1 \
			-g mm --keep-dup all \
			--outdir ./ \
			--nolambda --bdg --SPMR \
			--broad --broad-cutoff 0.1

sort -k1,1 -k2,2n H3K4me1_WT_r1_peaks.broadPeak > H3K4me1_WT_r1.sorted.broadPeak

echo -e "\n [`timestamp`] Finished peak calling for H3K4me1_WT_r1 \n"

conda deactivate




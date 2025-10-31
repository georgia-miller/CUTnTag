#!/bin/bash
#SBATCH --job-name=Script1_H3K4me1
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-users=k2477939@kcl.ac.uk


#######################################################
############ set basename and sample names ############
#######################################################

base_name=STM12023_WT_pBMDM

base_name_1=${base_name}_Rep1
base_name_2=${base_name}_Rep2
base_name_3=${base_name}_Rep3

#######################################################
###### check programme & other paths are correct ######
#######################################################

module load anaconda3/personal

PICARD=$HOME/anaconda3/envs/picard/share/picard-2.22.3-0/picard.jar

ATAC_dir_input=$RDS_PROJECT/hill_infection_regulatory_genomics/live/live/processed_data/Thurston_SteE_Stat3_paper_final_figures/processed_datasets_final/pBMDMs_ATACSeq/merged_fastqs
ATAC_dir_output=$RDS_PROJECT/hill_infection_regulatory_genomics/live/live/processed_data/Thurston_SteE_Stat3_paper_final_figures/processed_datasets_final/pBMDMs_ATACSeq/alignment_peakcalling
ATAC_index_genome=$HOME/genomes/mouse/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome

gene_coordinate=$HOME/genomes/mouse/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.normal_chr_only.bed
genome_file=$HOME/genomes/mouse/Mus_musculus/UCSC/mm10/genome.info
genome_fasta=$HOME/genomes/mouse/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa
blacklisted_mitochondrial_regions=$HOME/genomes/mouse/Mus_musculus/UCSC/mm10/mm10.blacklisted_and_chrM.sorted.bed
blacklisted_regions_only=$HOME/genomes/mouse/Mus_musculus/UCSC/mm10/mm10-blacklist.v2.bed
tss=$HOME/genomes/mouse/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/refTSS_v3.3_mouse_coordinate.mm10.bed

######################################################
########## align ATAC-Seq reads with bowtie ##########
######################################################


#############
### Rep 1 ###
#############

mkdir ${ATAC_dir_output}/${base_name_1}
mkdir ${ATAC_dir_output}/${base_name_1}/picard_temp
cd ${ATAC_dir_output}/${base_name_1}

#### trimming ####

source activate Trim_galore

trim_galore --paired --cores 4 --nextera ${ATAC_dir_input}/${base_name_1}_R1.fastq.gz ${ATAC_dir_input}/${base_name_1}_R2.fastq.gz

conda deactivate

#### alignment ####

source activate ATAC_pipeline_1

bowtie2 --threads 8 --very-sensitive -X 1000 -k 10 -x ${ATAC_index_genome} \
-1 ${base_name_1}_R1_val_1.fq.gz -2 ${base_name_1}_R2_val_2.fq.gz \
| samtools view -@ 8 -b -o ${base_name_1}.bam - #align ATAC seq with bowtie

conda deactivate

#### deduplicate ####

source activate picard

java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD SortSam INPUT=${base_name_1}.bam OUTPUT=${base_name_1}.picardchrsorted.bam \
SORT_ORDER=coordinate TMP_DIR=${ATAC_dir_output}/${base_name_1}/picard_temp VALIDATION_STRINGENCY=LENIENT

java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD MarkDuplicates INPUT=${base_name_1}.picardchrsorted.bam OUTPUT=${base_name_1}.deduplicated.bam \
TMP_DIR=${ATAC_dir_output}/${base_name_1}/picard_temp VALIDATION_STRINGENCY=LENIENT METRICS_FILE=${base_name_1}_PicardMarkDuplicates.txt REMOVE_DUPLICATES=true

conda deactivate

#### remove chrM and blacklist reads ####

source activate ATAC_pipeline_1

intersectBed -v -a ${base_name_1}.deduplicated.bam -b ${blacklisted_mitochondrial_regions} > \
${base_name_1}.deduplicated.cleaned.bam

samtools sort -o ${base_name_1}.deduplicated.cleaned.chrsorted.bam -T ${base_name_1}.deduplicated.cleaned.chrsorted -@ 16 \
${base_name_1}.deduplicated.cleaned.bam

samtools index ${base_name_1}.deduplicated.cleaned.chrsorted.bam #index sorted deduplicated bam

conda deactivate

#### cleanup unneeded files ####

rm ${base_name_1}_R1_val_1.fq.gz
rm ${base_name_1}_R2_val_2.fq.gz 
rm ${base_name_1}.bam
rm ${base_name_1}.deduplicated.bam
rm ${base_name_1}.deduplicated.cleaned.bam
rm -r ${ATAC_dir_output}/${base_name_1}/picard_temp

#############
### Rep 2 ###
#############

mkdir ${ATAC_dir_output}/${base_name_2}
mkdir ${ATAC_dir_output}/${base_name_2}/picard_temp
cd ${ATAC_dir_output}/${base_name_2}

#### trimming ####

source activate Trim_galore

trim_galore --paired --cores 4 --nextera ${ATAC_dir_input}/${base_name_2}_R1.fastq.gz ${ATAC_dir_input}/${base_name_2}_R2.fastq.gz

conda deactivate

#### alignment ####

source activate ATAC_pipeline_1

bowtie2 --threads 8 --very-sensitive -X 1000 -k 10 -x ${ATAC_index_genome} \
-1 ${base_name_2}_R1_val_1.fq.gz -2 ${base_name_2}_R2_val_2.fq.gz \
| samtools view -@ 8 -b -o ${base_name_2}.bam - #align ATAC seq with bowtie

conda deactivate

#### deduplicate ####

source activate picard

java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD SortSam INPUT=${base_name_2}.bam OUTPUT=${base_name_2}.picardchrsorted.bam \
SORT_ORDER=coordinate TMP_DIR=${ATAC_dir_output}/${base_name_2}/picard_temp VALIDATION_STRINGENCY=LENIENT

java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD MarkDuplicates INPUT=${base_name_2}.picardchrsorted.bam OUTPUT=${base_name_2}.deduplicated.bam \
TMP_DIR=${ATAC_dir_output}/${base_name_2}/picard_temp VALIDATION_STRINGENCY=LENIENT METRICS_FILE=${base_name_2}_PicardMarkDuplicates.txt REMOVE_DUPLICATES=true

conda deactivate

#### remove chrM and blacklist reads ####

source activate ATAC_pipeline_1

intersectBed -v -a ${base_name_2}.deduplicated.bam -b ${blacklisted_mitochondrial_regions} > \
${base_name_2}.deduplicated.cleaned.bam

samtools sort -o ${base_name_2}.deduplicated.cleaned.chrsorted.bam -T ${base_name_2}.deduplicated.cleaned.chrsorted -@ 16 \
${base_name_2}.deduplicated.cleaned.bam

samtools index ${base_name_2}.deduplicated.cleaned.chrsorted.bam #index sorted deduplicated bam

conda deactivate

#### cleanup unneeded files ####

rm ${base_name_2}_R1_val_1.fq.gz
rm ${base_name_2}_R2_val_2.fq.gz 
rm ${base_name_2}.bam
rm ${base_name_2}.deduplicated.bam
rm ${base_name_2}.deduplicated.cleaned.bam
rm -r ${ATAC_dir_output}/${base_name_2}/picard_temp

#############
### Rep 3 ###
#############

mkdir ${ATAC_dir_output}/${base_name_3}
mkdir ${ATAC_dir_output}/${base_name_3}/picard_temp
cd ${ATAC_dir_output}/${base_name_3}

#### trimming ####

source activate Trim_galore

trim_galore --paired --cores 4 --nextera ${ATAC_dir_input}/${base_name_3}_R1.fastq.gz ${ATAC_dir_input}/${base_name_3}_R2.fastq.gz

conda deactivate

#### alignment ####

source activate ATAC_pipeline_1

bowtie2 --threads 8 --very-sensitive -X 1000 -k 10 -x ${ATAC_index_genome} \
-1 ${base_name_3}_R1_val_1.fq.gz -2 ${base_name_3}_R2_val_2.fq.gz \
| samtools view -@ 8 -b -o ${base_name_3}.bam - #align ATAC-seq with bowtie

conda deactivate

#### deduplicate ####

source activate picard

java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD SortSam INPUT=${base_name_3}.bam OUTPUT=${base_name_3}.picardchrsorted.bam \
SORT_ORDER=coordinate TMP_DIR=${ATAC_dir_output}/${base_name_3}/picard_temp VALIDATION_STRINGENCY=LENIENT

java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD MarkDuplicates INPUT=${base_name_3}.picardchrsorted.bam OUTPUT=${base_name_3}.deduplicated.bam \
TMP_DIR=${ATAC_dir_output}/${base_name_3}/picard_temp VALIDATION_STRINGENCY=LENIENT METRICS_FILE=${base_name_3}_PicardMarkDuplicates.txt REMOVE_DUPLICATES=true

conda deactivate

#### remove chrM and blacklist reads ####

source activate ATAC_pipeline_1

intersectBed -v -a ${base_name_3}.deduplicated.bam -b ${blacklisted_mitochondrial_regions} > \
${base_name_3}.deduplicated.cleaned.bam

samtools sort -o ${base_name_3}.deduplicated.cleaned.chrsorted.bam -T ${base_name_3}.deduplicated.cleaned.chrsorted -@ 16 \
${base_name_3}.deduplicated.cleaned.bam

samtools index ${base_name_3}.deduplicated.cleaned.chrsorted.bam

conda deactivate

#### cleanup unneeded files ####

rm ${base_name_3}_R1_val_1.fq.gz
rm ${base_name_3}_R2_val_2.fq.gz 
rm ${base_name_3}.bam
rm ${base_name_3}.deduplicated.bam
rm ${base_name_3}.deduplicated.cleaned.bam
rm -r ${ATAC_dir_output}/${base_name_3}/picard_temp


#####################################################
############## peak calling with MACS2 ##############
#####################################################


source activate ATAC_pipeline_2

#####################
####### Rep1 ########
#####################

cd ${ATAC_dir_output}/${base_name_1}

macs2 callpeak -t ${base_name_1}.deduplicated.cleaned.chrsorted.bam -f BAMPE -n ${base_name_1} \
-g mm --keep-dup all -p 0.01 \
--outdir ${ATAC_dir_output}/${base_name_1}/ --nolambda --bdg --SPMR

sort -k1,1 -k2,2n ${base_name_1}_peaks.narrowPeak > ${base_name_1}.sorted.narrowPeak

#####################
####### Rep2 ########
#####################

cd ${ATAC_dir_output}/${base_name_2}

macs2 callpeak -t ${base_name_2}.deduplicated.cleaned.chrsorted.bam -f BAMPE -n ${base_name_2} \
-g mm --cutoff-analysis --keep-dup all -p 0.01 \
--outdir ${ATAC_dir_output}/${base_name_2}/ --nolambda --bdg --SPMR

sort -k1,1 -k2,2n ${base_name_2}_peaks.narrowPeak > ${base_name_2}.sorted.narrowPeak

#####################
####### Rep3 ########
#####################

cd ${ATAC_dir_output}/${base_name_3}

macs2 callpeak -t ${base_name_3}.deduplicated.cleaned.chrsorted.bam -f BAMPE -n ${base_name_3} \
-g mm --cutoff-analysis --keep-dup all -p 0.01 \
--outdir ${ATAC_dir_output}/${base_name_3}/ --nolambda --bdg --SPMR

sort -k1,1 -k2,2n ${base_name_3}_peaks.narrowPeak > ${base_name_3}.sorted.narrowPeak


#######################################################
############## merge ATAC-Seq replicates ##############
#######################################################


mkdir ${ATAC_dir_output}/${base_name}
cd ${ATAC_dir_output}/${base_name}

source activate ATAC_pipeline_1

samtools merge -@ 16 ${base_name}.merged.bam \
${ATAC_dir_output}/${base_name_1}/${base_name_1}.deduplicated.cleaned.chrsorted.bam \
${ATAC_dir_output}/${base_name_2}/${base_name_2}.deduplicated.cleaned.chrsorted.bam \
${ATAC_dir_output}/${base_name_3}/${base_name_3}.deduplicated.cleaned.chrsorted.bam

samtools sort -o ${base_name}.merged.chrsorted.bam -T ${base_name}.merged.chrsorted \
-@ 16 ${base_name}.merged.bam

samtools index ${base_name}.merged.chrsorted.bam 

rm ${base_name}.merged.bam

conda deactivate

######################################################
###### call peaks on merged ATAC-Seq replicates ######
######################################################

source activate ATAC_pipeline_2

macs2 callpeak -t ${base_name}.merged.chrsorted.bam -f BAMPE -n ${base_name} \
-g mm --keep-dup all -p 0.01 \
--outdir ${ATAC_dir_output}/${base_name}/ --nolambda --bdg --SPMR

sort -k1,1 -k2,2n ${base_name}_peaks.narrowPeak > ${base_name}.sorted.narrowPeak

conda deactivate

######################################################
############ identify reproducible peaks #############
######################################################

source activate bedtools

bedtools intersect -wa -a ${base_name}.sorted.narrowPeak \
-b ${ATAC_dir_output}/${base_name_1}/${base_name_1}.sorted.narrowPeak > \
${base_name}_merged.intersect_${base_name_1}.sorted.bed

bedtools intersect -wa -a ${base_name}_merged.intersect_${base_name_1}.sorted.bed \
-b ${ATAC_dir_output}/${base_name_2}/${base_name_2}.sorted.narrowPeak > \
${base_name}_merged.intersect_${base_name_1}_${base_name_2}.sorted.bed

bedtools intersect -wa -a ${base_name}_merged.intersect_${base_name_1}_${base_name_2}.sorted.bed \
-b ${ATAC_dir_output}/${base_name_3}/${base_name_3}.sorted.narrowPeak > \
${base_name}.consensus_peaks.stringent.intermediary.bed

bedtools merge -i ${base_name}.consensus_peaks.stringent.intermediary.bed > \
${base_name}.consensus_peaks.stringent.bed

rm ${base_name}_merged.intersect_${base_name_1}.sorted.bed
rm ${base_name}_merged.intersect_${base_name_1}_${base_name_2}.sorted.bed
rm ${base_name}.consensus_peaks.stringent.intermediary.bed
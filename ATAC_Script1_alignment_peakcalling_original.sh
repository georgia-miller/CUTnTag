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

base_name=H3K4me1

# define an array of bio replicate names that can be looped over
replicates=("${base_name}_Rep1" "${base_name}_Rep2" "${base_name}_Rep3")


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

# define timestamp to use for logs
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }


#############################################################
## align ATAC-Seq reads with bowtie & peak call with MACS2 ##
#############################################################

source activate CUTnTag_processing_env

for rep in "${replicates[@]}"; do
	mkdir -p ${ATAC_dir_output}/${rep}/picard_temp
	cd ${ATAC_dir_output}/${rep}

	echo "########[`timestamp`] Starting processing for CUT&Tag ${rep}########"

	#### trimming ####

	trim_galore --paired --cores 4 --nextera ${ATAC_dir_input}/${rep}_R1.fastq.gz ${ATAC_dir_input}/${rep}_R2.fastq.gz

	echo "[`timestamp`] Finished trimming for ${rep}"

	#### alignment ####

	bowtie2 --threads 8 --very-sensitive -X 1000 -k 10 -x ${ATAC_index_genome} \
	-1 ${rep}_R1_val_1.fq.gz -2 ${rep}_R2_val_2.fq.gz \
	| samtools view -@ 8 -b -o ${rep}.bam - #align ATAC seq with bowtie

	echo "[`timestamp`] Finished alignment for ${rep}"

	#### deduplicate ####

	java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD SortSam INPUT=${rep}.bam OUTPUT=${rep}.picardchrsorted.bam \
	SORT_ORDER=coordinate TMP_DIR=${ATAC_dir_output}/${rep}/picard_temp VALIDATION_STRINGENCY=LENIENT

	java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar $PICARD MarkDuplicates INPUT=${rep}.picardchrsorted.bam OUTPUT=${rep}.deduplicated.bam \
	TMP_DIR=${ATAC_dir_output}/${rep}/picard_temp VALIDATION_STRINGENCY=LENIENT METRICS_FILE=${rep}_PicardMarkDuplicates.txt REMOVE_DUPLICATES=true

	echo "[`timestamp`] Finished marking duplicates for ${rep}"

	#### remove chrM and blacklist reads ####

	intersectBed -v -a ${rep}.deduplicated.bam -b ${blacklisted_mitochondrial_regions} > \
	${rep}.deduplicated.cleaned.bam

	samtools sort -o ${rep}.deduplicated.cleaned.chrsorted.bam -T ${rep}.deduplicated.cleaned.chrsorted -@ 16 \
	${rep}.deduplicated.cleaned.bam

	samtools index ${rep}.deduplicated.cleaned.chrsorted.bam #index sorted deduplicated bam

	echo "[`timestamp`] Finished filtering for ${rep}"

	#### cleanup unneeded files ####

	rm  ${rep}_R1_val_1.fq.gz \
		${rep}_R2_val_2.fq.gz \
		${rep}.bam \
		${rep}.deduplicated.bam \
		${rep}.deduplicated.cleaned.bam
	rm -r ${ATAC_dir_output}/${rep}/picard_temp

	echo "[`timestamp`] Finished processing ${rep}"

	echo "[`timestamp`] Starting peak calling for ${rep}"
	
	macs2 callpeak -t ${rep}.deduplicated.cleaned.chrsorted.bam -f BAMPE -n ${rep} \
	-g mm --keep-dup all -p 0.01 \
	--outdir ${ATAC_dir_output}/${rep}/ --nolambda --bdg --SPMR

	sort -k1,1 -k2,2n ${rep}_peaks.narrowPeak > ${rep}.sorted.narrowPeak

	echo "[`timestamp`] Finished peak calling for ${rep}"


	echo "########[`timestamp`] Finished for ${rep} output found in ${ATAC_dir_output}/${rep}########"


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
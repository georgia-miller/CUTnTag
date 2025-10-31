#!/bin/bash
#SBATCH --job-name=Script1_H3K4me1_UI #cleaned_cut&tag_example_script_before.sep.conda.envs
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=Script1_H3K4me1_UI_%j.log
#SBATCH --error=Script1_H3K4me1_UI_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-users=k2477939@kcl.ac.uk


#### MUST CHANGE: SCRIPT AND BASE NAME, DIRECTORY PATHS, MACS2 ARGS DEPENDENT ON TAG ####

#######################################################
############ set basename and sample names ############
#######################################################

base_name=H3K4me1_UI

# define an array of bio replicate names that can be looped over
replicates=("${base_name}_r1" "${base_name}_r2" "${base_name}_r3")


#######################################################
###### check programme & other paths are correct ######
#######################################################

PICARD="${CONDA_PREFIX}/share/picard/picard.jar"

dir_input=/scratch/prj/id_hill_sims_wellcda/CUTnTag/merged_fastqs
dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling
index_genome=/scratch/prj/id_hill_sims_wellcda/CUTnTag/ref_genome/genome # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/Sequence/Bowtie2Index
blacklisted_mitochondrial_regions=/scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.blacklisted_and_chrM.sorted.bed # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/mm10.blacklisted_and_chrM.sorted.bed 

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# load conda
module load anaconda3/2022.10-gcc-13.2.0

#######################################################
####### per replicate: align reads & call peaks #######
#######################################################

conda activate CUTnTag_processing_env
conda list --name CUTnTag_processing_env # list installed packages and versions

for rep in "${replicates[@]}"; do

	#### make a directory per rep and change into it ####
	rep_dir=${dir_output}/${rep}
	mkdir -p ${rep_dir}/picard_temp
	cd ${rep_dir}

	echo "######## [`timestamp`] Process and call peaks for CUT&Tag ${rep} ########"

	#### trimming ####
	trim_galore --paired --cores 4 --nextera ${dir_input}/${rep}*_R1_merged.fastq.gz ${dir_input}/${rep}*_R2_merged.fastq.gz

	echo "[`timestamp`] Finished trimming for ${rep}"

	#### alignment ####

	bowtie2 --threads 8 --very-sensitive -X 1000 -k 10 -x ${index_genome} \
		-1 ${rep}_R1_val_1.fq.gz -2 ${rep}_R2_val_2.fq.gz \
		| samtools view -@ 8 -b -o ${rep}.bam - 

	echo "[`timestamp`] Finished alignment for ${rep}"

	#### mark duplicates ####

	java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar ${PICARD} SortSam INPUT=${rep}.bam OUTPUT=${rep}.picardchrsorted.bam \
		SORT_ORDER=coordinate TMP_DIR=${rep_dir}/picard_temp VALIDATION_STRINGENCY=LENIENT

	# remove duplicates is set to false so they are only marked
	java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar ${PICARD} MarkDuplicates INPUT=${rep}.picardchrsorted.bam OUTPUT=${rep}.marked.bam \
		TMP_DIR=${rep_dir}/picard_temp VALIDATION_STRINGENCY=LENIENT METRICS_FILE=${rep}_PicardMarkDuplicates.txt REMOVE_DUPLICATES=false

	echo "[`timestamp`] Finished marking duplicates for ${rep} found at ${rep}_PicardMarkDuplicates.txt"

	#### remove chrM and blacklist reads ####

	intersectBed -v -a ${rep}.marked.bam -b ${blacklisted_mitochondrial_regions} > \
		${rep}.marked.cleaned.bam

	samtools sort -o ${rep}.marked.cleaned.chrsorted.bam -T ${rep}.marked.cleaned.chrsorted -@ 16 \
		${rep}.marked.cleaned.bam

	samtools index ${rep}.marked.cleaned.chrsorted.bam #index sorted marked bam

	echo "[`timestamp`] Finished filtering for ${rep}"

	#### cleanup unneeded files ####

	rm  ${rep}_R1_val_1.fq.gz \
		${rep}_R2_val_2.fq.gz \
		${rep}.bam \
		${rep}.marked.bam \
		${rep}.marked.cleaned.bam
	rm -r ${rep_dir}/picard_temp


	echo "######## [`timestamp`] Starting peak calling for ${rep}} ########"
	
	macs2 callpeak -t ${rep}.marked.cleaned.chrsorted.bam -f BAMPE -n ${rep} \
		-g mm --keep-dup all -p 0.01 \
		--outdir ${rep_dir}/ --nolambda --bdg --SPMR #can add if needed: --cutoff-analysis  

	sort -k1,1 -k2,2n ${rep}_peaks.narrowPeak > ${rep}.sorted.narrowPeak

	echo "[`timestamp`] Finished peak calling for ${rep}"


	echo "[`timestamp`] Finished for ${rep} output found in ${rep_dir}"
done

######################################################
############## merge CUT&Tag replicates ##############
######################################################

echo "######## [`timestamp`] Merge replicates and call consensus peaks ########"

mkdir -p ${dir_output}/${base_name}
cd ${dir_output}/${base_name}

# here e.g. replicates[0] refers to 1st item of the replicates array = ${base_name}_r1
samtools merge -@ 16 ${base_name}.merged.bam \
	${dir_output}/${replicates[0]}/${replicates[0]}.marked.cleaned.chrsorted.bam \
	${dir_output}/${replicates[1]}/${replicates[1]}.marked.cleaned.chrsorted.bam \
	${dir_output}/${replicates[2]}/${replicates[2]}.marked.cleaned.chrsorted.bam

samtools sort -o ${base_name}.merged.chrsorted.bam -T ${base_name}.merged.chrsorted \
	-@ 16 ${base_name}.merged.bam

samtools index ${base_name}.merged.chrsorted.bam 

rm ${base_name}.merged.bam

echo "[`timestamp`] Finished merging replicates for ${base_name}"

#####################################################
###### call peaks on merged CUT&Tag replicates ######
#####################################################

macs2 callpeak -t ${base_name}.merged.chrsorted.bam -f BAMPE -n ${base_name} \
	-g mm --keep-dup all -p 0.01 \
	--outdir ${dir_output}/${base_name}/ --nolambda --bdg --SPMR

sort -k1,1 -k2,2n ${base_name}_peaks.narrowPeak > ${base_name}.sorted.narrowPeak

echo "[`timestamp`] Finished calling peaks for ${base_name}"

######################################################
############ identify reproducible peaks #############
######################################################

bedtools intersect -wa -a ${base_name}.sorted.narrowPeak \
	-b ${dir_output}/${replicates[0]}/${replicates[0]}.sorted.narrowPeak > \
	${base_name}_merged.intersect_${replicates[0]}.sorted.bed

bedtools intersect -wa -a ${base_name}_merged.intersect_${replicates[0]}.sorted.bed \
	-b ${dir_output}/${replicates[1]}/${replicates[1]}.sorted.narrowPeak > \
	${base_name}_merged.intersect_${replicates[0]}_${replicates[1]}.sorted.bed

bedtools intersect -wa -a ${base_name}_merged.intersect_${replicates[0]}_${replicates[1]}.sorted.bed \
	-b ${dir_output}/${replicates[2]}/${replicates[2]}.sorted.narrowPeak > \
	${base_name}.consensus_peaks.stringent.intermediary.bed

bedtools merge -i ${base_name}.consensus_peaks.stringent.intermediary.bed > \
	${base_name}.consensus_peaks.stringent.bed

rm 	${base_name}_merged.intersect_${replicates[0]}.sorted.bed \
	${base_name}_merged.intersect_${replicates[0]}_${replicates[1]}.sorted.bed \
	${base_name}.consensus_peaks.stringent.intermediary.bed


echo "[`timestamp`] Finished identifying reproducible peaks for ${base_name}"

conda deactivate

echo "######## [`timestamp`] Script completed ########"



#!/bin/bash
#SBATCH --job-name=Script1_H3K4me1_test
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script1_H3K4me1_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script1_H3K4me1_%j.log


#### MUST CHANGE: SCRIPT AND BASE NAME, DIRECTORY PATHS, MACS2 ARGS DEPENDENT ON TAG ####

#######################################################
################ set names & parameters ###############
#######################################################

modification=H3K4me1

if [[ ${modification} == "H3K4me1" ]]; then
	macs2_args="--broad --broad-cutoff 0.1"
    peak_type="broadPeak"
elif [[ ${modification} == "H3K4me3" ]]; then
	macs2_args="-q 0.01"
    peak_type="narrowPeak"
elif [[ ${modification} == "H3K27ac" ]]; then
	macs2_args="-q 0.05"
    peak_type="narrowPeak"
elif [[ ${modification} == "BRG1" ]]; then
	macs2_args="-q 0.05"
    peak_type="narrowPeak"
else
    echo -e "\n######## [`timestamp`] Error: Modification ${modification} is not accepted ######## \n"
    exit 1
fi


dir_input=/scratch/prj/id_hill_sims_wellcda/CUTnTag/merged_fastqs
dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling
index_genome=/scratch/prj/id_hill_sims_wellcda/CUTnTag/ref_genome/genome # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/Sequence/Bowtie2Index
blacklisted_mitochondrial_regions=/scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.blacklisted_and_chrM.sorted.bed # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/mm10.blacklisted_and_chrM.sorted.bed 

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# load the Anaconda module
module load anaconda3/2022.10-gcc-13.2.0

# source conda so it works in non-interactive shells
source $(dirname $(which conda))/../etc/profile.d/conda.sh

#######################################################
############ start to loop over conditions ############
#######################################################

# define an array of base names that can be looped over
conditions=("${modification}_IL10" "${modification}_SteE" "${modification}_UI" "${modification}_WT")


for base_name in "${conditions[@]}"; do

	echo -e "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ [`timestamp`] Starting ${base_name} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n"

	# define an array of bio replicate names that can be looped over
	replicates=("${base_name}_r1" "${base_name}_r2" "${base_name}_r3")


	#######################################################
	############ per replicate: read alignment ############
	#######################################################

	conda activate CUTnTag_alignment_env
	echo -e "\n ######## [`timestamp`] Active environment: $(basename $CONDA_PREFIX) ######## \n"
	#conda list --name CUTnTag_alignment_env # list installed packages and versions

	PICARD="${CONDA_PREFIX}/share/picard/picard.jar"

	for rep in "${replicates[@]}"; do

		#### make a directory per rep and change into it ####
		rep_dir=${dir_output}/${rep}
		mkdir -p ${rep_dir}/picard_temp
		cd ${rep_dir}

		#echo -e "\n ######## [`timestamp`] Process and call peaks for CUT&Tag ${rep} ######## \n"

		#### trimming ####
		#trim_galore --paired --cores 4 --nextera ${dir_input}/${rep}*_R1_merged.fastq.gz ${dir_input}/${rep}*_R2_merged.fastq.gz

		#echo -e "\n [`timestamp`] Finished trimming for ${rep} \n"

		#### alignment ####

		#bowtie2 --threads 8 --very-sensitive -X 1000 -k 10 -x ${index_genome} \
		#	-1 ${rep}_R1_merged_val_1.fq.gz -2 ${rep}_R2_merged_val_2.fq.gz \
		#	| samtools view -@ 8 -b -o ${rep}.bam - 

		echo -e "\n [`timestamp`] Finished alignment for ${rep} \n"

		#### mark duplicates ####

		#java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar ${PICARD} SortSam INPUT=${rep}.bam OUTPUT=${rep}.picardchrsorted.bam \
		#	SORT_ORDER=coordinate TMP_DIR=${rep_dir}/picard_temp VALIDATION_STRINGENCY=LENIENT

		# remove duplicates is set to false so they are only marked
		#java -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8 -Xmx30g -jar ${PICARD} MarkDuplicates INPUT=${rep}.picardchrsorted.bam OUTPUT=${rep}.marked.bam \
		#	TMP_DIR=${rep_dir}/picard_temp VALIDATION_STRINGENCY=LENIENT METRICS_FILE=${rep}_PicardMarkDuplicates.txt REMOVE_DUPLICATES=false

		echo -e "\n [`timestamp`] Finished marking duplicates for ${rep} found at ${rep}_PicardMarkDuplicates.txt \n"

		#### remove chrM and blacklist reads ####

		#intersectBed -v -a ${rep}.marked.bam -b ${blacklisted_mitochondrial_regions} > \
		#	${rep}.marked.cleaned.bam

		#samtools sort -o ${rep}.marked.cleaned.chrsorted.bam -T ${rep}.marked.cleaned.chrsorted -@ 16 \
		#	${rep}.marked.cleaned.bam

		#samtools index ${rep}.marked.cleaned.chrsorted.bam #index sorted marked bam

		echo -e "\n [`timestamp`] Finished filtering for ${rep} \n"

		#### cleanup unneeded files ####

		#rm  ${rep}_R1_val_1.fq.gz \
		#	${rep}_R2_val_2.fq.gz \
		#	${rep}.bam \
		#	${rep}.marked.bam \
		#	${rep}.marked.cleaned.bam
		#rm -r ${rep_dir}/picard_temp
	done

	conda deactivate

	echo -e "\n ######## [`timestamp`] Finished all processing for ${base_name} ######## \n"

	#######################################################
	############# per replicate: peak calling #############
	#######################################################

	conda activate CUTnTag_macs2_env
	echo -e "\n ######## [`timestamp`] Active environment: $(basename $CONDA_PREFIX) ######## \n"
	#conda list --name CUTnTag_alignment_env # list installed packages and versions

	for rep in "${replicates[@]}"; do

		#### change back into the rep specific directory ####
		rep_dir=${dir_output}/${rep}
		cd ${rep_dir}

		echo -e "\n ######## [`timestamp`] Starting peak calling for ${rep}} ######## \n"
		
		#macs2 callpeak -t ${rep}.marked.cleaned.chrsorted.bam -f BAMPE -n ${rep} \
		#	-g mm --keep-dup all --outdir ${rep_dir}/ --nolambda --bdg --SPMR \
		#	${macs2_args} # these will change depending on the modification/tag
		#	# can add if needed: --cutoff-analysis

		#sort -k1,1 -k2,2n ${rep}_peaks.${peak_type} > ${rep}.sorted.${peak_type}

		echo -e "\n [`timestamp`] Finished peak calling for ${rep} \n"


		echo -e "\n [`timestamp`] Finished for ${rep} output found in ${rep_dir} \n"
	done


	######################################################
	############## merge CUT&Tag replicates ##############
	######################################################

	echo -e "\n ######## [`timestamp`] Merge replicates and call consensus peaks ######## \n"

	mkdir -p ${dir_output}/${base_name}
	cd ${dir_output}/${base_name}

	# here e.g. replicates[0] refers to 1st item of the replicates array = ${base_name}_r1
	#samtools merge -@ 16 ${base_name}.merged.bam \
	#	${dir_output}/${replicates[0]}/${replicates[0]}.marked.cleaned.chrsorted.bam \
	#	${dir_output}/${replicates[1]}/${replicates[1]}.marked.cleaned.chrsorted.bam \
	#	${dir_output}/${replicates[2]}/${replicates[2]}.marked.cleaned.chrsorted.bam

	#samtools sort -o ${base_name}.merged.chrsorted.bam -T ${base_name}.merged.chrsorted \
	#	-@ 16 ${base_name}.merged.bam

	#samtools index ${base_name}.merged.chrsorted.bam 

	#rm ${base_name}.merged.bam

	echo -e "\n [`timestamp`] Finished merging replicates for ${base_name} \n"

	#####################################################
	###### call peaks on merged CUT&Tag replicates ######
	#####################################################

	#macs2 callpeak -t ${base_name}.merged.chrsorted.bam -f BAMPE -n ${base_name} \
	#	-g mm --keep-dup all --outdir ${dir_output}/${base_name}/ --nolambda --bdg --SPMR \
	#	${macs2_args} # these will change depending on the modification/tag
	#	# can add if needed: --cutoff-analysis

	#sort -k1,1 -k2,2n ${base_name}_peaks.${peak_type} > ${base_name}.sorted.${peak_type}

	echo -e "\n [`timestamp`] Finished calling peaks for ${base_name} \n"

	######################################################
	############ identify reproducible peaks #############
	######################################################

	#bedtools intersect -wa -a ${base_name}.sorted.${peak_type} \
	#	-b ${dir_output}/${replicates[0]}/${replicates[0]}.sorted.${peak_type} > \
	#	${base_name}_merged.intersect_${replicates[0]}.sorted.bed

	#bedtools intersect -wa -a ${base_name}_merged.intersect_${replicates[0]}.sorted.bed \
	#	-b ${dir_output}/${replicates[1]}/${replicates[1]}.sorted.${peak_type} > \
	#	${base_name}_merged.intersect_${replicates[0]}_${replicates[1]}.sorted.bed

	#bedtools intersect -wa -a ${base_name}_merged.intersect_${replicates[0]}_${replicates[1]}.sorted.bed \
	#	-b ${dir_output}/${replicates[2]}/${replicates[2]}.sorted.${peak_type} > \
	#	${base_name}.consensus_peaks.stringent.intermediary.bed

	#bedtools merge -i ${base_name}.consensus_peaks.stringent.intermediary.bed > \
	#	${base_name}.consensus_peaks.stringent.bed

	#rm 	${base_name}_merged.intersect_${replicates[0]}.sorted.bed \
	#	${base_name}_merged.intersect_${replicates[0]}_${replicates[1]}.sorted.bed \
	#	${base_name}.consensus_peaks.stringent.intermediary.bed


	echo -e "\n [`timestamp`] Finished identifying reproducible peaks for ${base_name} \n"

	conda deactivate

	echo -e "\n ######## [`timestamp`] ${base_name} completed ######## \n"
done


echo -e "\n ######## [`timestamp`] Script completed for ${modification} ######## \n"



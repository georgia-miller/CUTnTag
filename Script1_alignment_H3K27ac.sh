#!/bin/bash
#SBATCH --job-name=Script1_H3K27ac
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script1_H3K27ac_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script1_H3K27ac_%j.log


#### MUST CHANGE: SCRIPT AND MODIFICATION NAME & CHECK DIRECTORIES ####

#######################################################
################ set names & parameters ###############
#######################################################

modification=H3K27ac

echo -e "\n######## Script 1 to do alignments for each replicate for ${modification} ######## \n"

dir_input=/scratch/prj/id_hill_sims_wellcda/CUTnTag/merged_fastqs
dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling/${modification}
index_genome=/scratch/prj/id_hill_sims_wellcda/CUTnTag/ref_genome/genome # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/Sequence/Bowtie2Index
blacklisted_mitochondrial_regions=/scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.blacklisted_and_chrM.sorted.bed # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/mm10.blacklisted_and_chrM.sorted.bed 

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# load the Anaconda module
module load anaconda3/2022.10-gcc-13.2.0

# source conda so it works in non-interactive shells
#source $(dirname $(which conda))/../etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

# activate the alignment conda environment
conda activate CUTnTag_alignment_env
echo -e "\n ######## [$(timestamp)] Active environment: $(basename $CONDA_PREFIX) ######## \n"
#conda list --name CUTnTag_alignment_env # list installed packages and versions

trim_galore --version | head -n 3
bowtie2 --version | head -n 1

# set java options so the picard arguments run with these
export JAVA_OPTS="-Xmx30g -XX:ParallelGCThreads=8 -XX:ConcGCThreads=8"

# set an error trap
set -o errexit # stop if encounter an error
set -o pipefail # stop on pipeline errors (if any command within a piped command fails)
trap 'echo -e "\n [$(timestamp)] Error in ${BASH_COMMAND}. Exiting. \n" >&2; exit 1' ERR # when a command errors, run the command inside the quotes (print timestamp, write the failed command, send the message to the log file and stop the script)

#######################################################
############ start to loop over conditions ############
#######################################################

# define an array of base names that can be looped over
conditions=("${modification}_IL10" "${modification}_SteE" "${modification}_UI" "${modification}_WT")


for base_name in "${conditions[@]}"; do

	echo -e "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ [$(timestamp)] Starting ${base_name} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n"

	# define an array of bio replicate names that can be looped over
	replicates=("${base_name}_r1" "${base_name}_r2" "${base_name}_r3")


	#######################################################
	############ per replicate: read alignment ############
	#######################################################

	for rep in "${replicates[@]}"; do

		#### make a directory per rep and change into it ####
		rep_dir=${dir_output}/${rep}
		mkdir -p ${rep_dir}/picard_temp
		cd ${rep_dir}

		echo -e "\n ######## [$(timestamp)] Starting processing for CUT&Tag ${rep} ######## \n"

		#### trimming ####
		trim_galore --cores 4 --paired --nextera ${dir_input}/${rep}*_R1_merged.fastq.gz ${dir_input}/${rep}*_R2_merged.fastq.gz

		echo -e "\n [$(timestamp)] Finished trimming for ${rep} \n"

		#### alignment ####

		bowtie2 --threads 8 --very-sensitive -X 1000 -k 10 \
			-x ${index_genome} \
			-1 ${rep}_R1_merged_val_1.fq.gz \
			-2 ${rep}_R2_merged_val_2.fq.gz \
			2> ${rep}_bowtie2.log \
			| samtools view -@ 8 -b -o ${rep}.bam - 

		echo -e "\n [$(timestamp)] Finished alignment for ${rep} \n"

		#### mark duplicates ####

		# add dummy read groups so MarkDuplicate works
		picard AddOrReplaceReadGroups \
			-I ${rep}.bam \
			-O ${rep}.RG.bam \
			--RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM ${rep}

		picard SortSam \
			-I ${rep}.RG.bam \
			-O ${rep}.picardchrsorted.bam \
			-SORT_ORDER coordinate \
			--TMP_DIR ${rep_dir}/picard_temp \
			--VALIDATION_STRINGENCY "LENIENT"

		# remove duplicates is set to false so they are only marked
		picard MarkDuplicates \
			-I ${rep}.picardchrsorted.bam \
			-O ${rep}.marked.bam \
			--TMP_DIR ${rep_dir}/picard_temp \
			--VALIDATION_STRINGENCY "LENIENT" \
			-METRICS_FILE ${rep}_PicardMarkDuplicates.txt \
			--REMOVE_DUPLICATES false

		echo -e "\n [$(timestamp)] Finished marking duplicates for ${rep} found at ${rep}_PicardMarkDuplicates.txt \n"

		#### remove chrM and blacklist reads ####

		intersectBed -v \
			-a ${rep}.marked.bam \
			-b ${blacklisted_mitochondrial_regions} > \
			${rep}.marked.cleaned.bam

		samtools sort -o ${rep}.marked.cleaned.chrsorted.bam \
			-T ${rep}.marked.cleaned.chrsorted \
			-@ 16 \
			${rep}.marked.cleaned.bam

		 # index the sorted duplication marked bam
		samtools index ${rep}.marked.cleaned.chrsorted.bam

		echo -e "\n [$(timestamp)] Finished filtering for ${rep} \n"

		#### cleanup unneeded files ####

		rm  ${rep}_R1_merged_val_1.fq.gz \
			${rep}_R2_merged_val_2.fq.gz \
			${rep}.bam \
			${rep}.RG.bam \
			${rep}.picardchrsorted.bam \
			${rep}.marked.bam \
			${rep}.marked.cleaned.bam
		rm -r ${rep_dir}/picard_temp
	done


	echo -e "\n ######## [$(timestamp)] Finished all processing for ${base_name} ######## \n"

done


conda deactivate


echo -e "\n ######## [$(timestamp)] Script completed for ${modification} ######## \n"



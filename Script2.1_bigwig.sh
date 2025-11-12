#!/bin/bash
#SBATCH --job-name=Script2.1_bigwig
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script2.1_bigwig_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script2.1_bigwig_%j.log


#######################################################
################ set names & parameters ###############
#######################################################


echo -e "\n######## Script 2.1 to make normalised bigwig files ######## \n"

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

#define base directory output name
dir=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling

# load the Anaconda module
module load anaconda3/2022.10-gcc-13.2.0

# source conda so it works in non-interactive shells
#source $(dirname $(which conda))/../etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

# activate the SEACR conda environment
conda activate CUTnTag_seacr_env
echo -e "\n #### [$(timestamp)] Active environment: $(basename $CONDA_PREFIX) #### \n"

# set an error trap
set -o errexit # stop if encounter an error
set -o pipefail # stop on pipeline errors (if any command within a piped command fails)
trap 'echo -e "\n [$(timestamp)] Error in ${BASH_COMMAND}. Exiting. \n" >&2; exit 1' ERR # when a command errors, run the command inside the quotes (print timestamp, write the failed command, send the message to the log file and stop the script)


#######################################################
########### start to loop over modifications ##########
#######################################################

# define an array that can be looped over
modifications=("H3K4me1" "H3K4me3" "H3K27ac" "BRG1")

for mod in "${modifications[@]}"; do

echo -e "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ [$(timestamp)] Starting ${mod} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n"

	cd ${dir}/${mod}/

	# define an array of base names that can be looped over
	conditions=("${mod}_IL10" "${mod}_SteE" "${mod}_UI" "${mod}_WT")

	for base_name in "${conditions[@]}"; do

		echo -e "\n ######## Starting ${base_name} ######## \n"

		# define an array of bio replicate names that can be looped over
		replicates=("${base_name}_r1" "${base_name}_r2" "${base_name}_r3")

		rep1=${dir}/${mod}/${replicates[0]}
		rep2=${dir}/${mod}/${replicates[1]}
		rep3=${dir}/${mod}/${replicates[2]}

		# count usable reads in each replicate
		count1=$(samtools view -c -f 0x4 "${rep1}/${replicates[0]}.marked.cleaned.chrsorted.bam")
		count2=$(samtools view -c -f 0x4 "${rep2}/${replicates[1]}.marked.cleaned.chrsorted.bam")
		count3=$(samtools view -c -f 0x4 "${rep3}"/${replicates[2]}.marked.cleaned.chrsorted.bam)

		echo -e "\n Counts for rep 1: ${counts1} \n for rep 2: ${counts2} \n for rep 3: ${counts3} \n"

		# find minimum reads
		minCount=$(printf "%s\n" ${count1} ${count2} ${count3} | sort -n | head -n)
		echo -e "\n The minimum count is: ${minCount} /n"

		for rep in "${replicates[@]}"; do

		# downsample each replicate
		samtools view -@ 16 \
			-b \
			--subsample-seed 42 \
			--subsample-count ${min_reads} \
			-o ${rep}/${rep}.merged.ds.bam \
			${rep}/${rep}.ds.bam

		echo -e "\n Downsampled ${rep} to $(samtools view -c -f 0x2 ${rep}/${rep}.ds.bam) reads \n"

		done

		# merge the 3 replicates
		samtools merge -@ 16 ${base_name}.merged.ds.bam \
			${rep1} \
			${rep2} \
			${rep3}

		samtools sort -o ${base_name}.merged.ds.chrsorted.bam \
		-T ${base_name}.merged.ds.chrsorted \
		-@ 16 ${base_name}.merged.ds.bam

		samtools index ${base_name}.merged.ds.chrsorted.bam 

		# convert to bigwig
		bamCoverage -b ${base_name}.merged.ds.chrsorted.bam \
			-o ${base_name}.merged.ds.bw \
  			--normalizeUsing CPM \
  			--binSize 10

		# for rep in "${replicates[@]}"; do

		# 	# change into the rep specific directory
		# 	cd ${dir}/${mod}/${rep}

		# 	echo -e "\n ######## [$(timestamp)] Converting ${rep}.bam to bigwig ######## \n"
		
		# 	bamCoverage -b ${rep}.marked.cleaned.chrsorted.bam \
  		# 	-o ${rep}.bw \
  		# 	--normalizeUsing CPM \
  		# 	--binSize 10

		# 	echo -e "\n [$(timestamp)] Finished \n"
		# done

		rm  ${base_name}.merged.ds.bam 
			# add downsampled!

	done

	echo -e "\n ######## [$(timestamp)] Finished ${mod} ######## \n"

done


conda deactivate

echo -e "\n ######## [$(timestamp)] Script completed ######## \n"



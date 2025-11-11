#!/bin/bash
#SBATCH --job-name=Script2.1_bigwig
#SBATCH --time=06:00:00
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
dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling

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

	# define an array of base names that can be looped over
	conditions=("${mod}_IL10" "${mod}_SteE" "${mod}_UI" "${mod}_WT")

	for base_name in "${conditions[@]}"; do

		# define an array of bio replicate names that can be looped over
		replicates=("${base_name}_r1" "${base_name}_r2" "${base_name}_r3")
	
		for rep in "${replicates[@]}"; do

			# change into the rep specific directory
			cd ${dir_output}/${mod}/${rep}

			echo -e "\n ######## [$(timestamp)] Converting ${rep}.bam to bigwig ######## \n"
		
			bamCoverage -b ${rep}.marked.cleaned.chrsorted.bam \
  			-o ${rep}.bw \
  			--normalizeUsing CPM \
  			--binSize 10

			echo -e "\n [$(timestamp)] Finished \n"
		done

	done

	echo -e "\n ######## [$(timestamp)] Finished ${mod} ######## \n"

done


conda deactivate

echo -e "\n ######## [$(timestamp)] Script completed ######## \n"



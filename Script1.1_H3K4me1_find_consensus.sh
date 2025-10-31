#!/bin/bash
#SBATCH --job-name=Script1.1_H3K4me1_find_consensus
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script1_H3K4me1_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script1_H3K4me1_%j.log


#### MUST CHANGE: SCRIPT AND MODIFICATION NAME, CHECK DIRECTORIES & MACS2 ARGS ####

echo -e "\n ### Script1.1 Taking processed H3K4me1 peak files per replicate and finding consensus files in =>2 replicates, not all 3 like currently ### \n"

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
dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling/${modification}/
index_genome=/scratch/prj/id_hill_sims_wellcda/CUTnTag/ref_genome/genome # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/Sequence/Bowtie2Index
blacklisted_mitochondrial_regions=/scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.blacklisted_and_chrM.sorted.bed # from: /rds/prj/id_hill_sims_wellcda/genomes/mm10/mm10.blacklisted_and_chrM.sorted.bed 

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# load the Anaconda module
module load anaconda3/2022.10-gcc-13.2.0

# source conda so it works in non-interactive shells
#source $(dirname $(which conda))/../etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

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
	############# per replicate: peak calling #############
	#######################################################

	conda activate CUTnTag_macs2_env
	echo -e "\n ######## [`timestamp`] Active environment: $(basename $CONDA_PREFIX) ######## \n"


	######################################################
	############## merge CUT&Tag replicates ##############
	######################################################

	echo -e "\n ######## [`timestamp`] Merge replicates and call consensus peaks ######## \n"

	mkdir -p ${dir_output}/${base_name}
	cd ${dir_output}/${base_name}



	# create list of all replicate peak files
	peak_files=("${dir_output}/${replicates[0]}/${replicates[0]}.sorted.${peak_type}" \
				"${dir_output}/${replicates[1]}/${replicates[1]}.sorted.${peak_type}" \
				"${dir_output}/${replicates[2]}/${replicates[2]}.sorted.${peak_type}")

	# count peaks between replicates to find overlapping ones
	bedtools multiinter -i ${peak_files[@]} > ${base_name}.multiinter.bed

	# filter to keep peaks present in at least 2 replicates
	awk '$4>=2 {print $1"\t"$2"\t"$3"\t"$4}' ${base_name}.multiinter.bed > ${base_name}.consensus_peaks_raw.bed 

	# merge any overlapping/bookmarkes peaks
	bedtools merge -i ${base_name}.consensus_peaks_raw.bed > ${base_name}.consensus_peaks.bed 

	echo -e "\n [`timestamp`] Finished identifying reproducible peaks for ${base_name} \n"

	conda deactivate

	echo -e "\n ######## [`timestamp`] ${base_name} completed ######## \n"
done


echo -e "\n ######## [`timestamp`] Script completed for ${modification} ######## \n"



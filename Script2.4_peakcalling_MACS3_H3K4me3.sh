#!/bin/bash
#SBATCH --job-name=Script2_H3K4me3
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script2.4_MACS3_H3K4me3_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script2.4_MACS3_H3K4me3_%j.log


#######################################################
################ set names & parameters ###############
#######################################################

modification=H3K4me3

echo -e "\n######## Script 2.4 to call peaks with MACS3 for each replicate and merged samples for ${modification} ######## \n"

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

if [[ ${modification} == "H3K4me1" ]]; then
	macs3_args="--broad --broad-cutoff 0.1"
    peak_type="broadPeak"
elif [[ ${modification} == "H3K4me3" ]]; then
	macs3_args="-q 0.01"
    peak_type="narrowPeak"
elif [[ ${modification} == "H3K27ac" ]]; then
	macs3_args="-q 0.05"
    peak_type="narrowPeak"
elif [[ ${modification} == "BRG1" ]]; then
	macs3_args="-q 0.05"
    peak_type="narrowPeak"
else
    echo -e "\n######## [${timestamp}] Error: Modification ${modification} is not accepted ######## \n"
    exit 1
fi

dir_input=/scratch/prj/id_hill_sims_wellcda/CUTnTag/merged_fastqs
dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling/${modification}
mm10_genome=/scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.genome.size

# source conda so it works in non-interactive shells
#source $(dirname $(which conda))/../etc/profile.d/conda.sh
source /users/k2477939/conda/etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

# set an error trap
set -o errexit # stop if encounter an error
set -o pipefail # stop on pipeline errors (if any command within a piped command fails)
trap 'echo -e "\n [$(timestamp)] Error in ${BASH_COMMAND}. Exiting. \n" >&2; exit 1' ERR # when a command errors, run the command inside the quotes (print timestamp, write the failed command, send the message to the log file and stop the script)


#######################################################
############ start to loop over conditions ############
#######################################################

# define an array of base names that can be looped over
conditions=("${modification}_UI" "${modification}_WT" "${modification}_SteE" "${modification}_IL10")

for base_name in "${conditions[@]}"; do

	echo -e "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ [${timestamp}] Starting ${base_name} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n"

	# define an array of bio replicate names that can be looped over
	replicates=("${base_name}_r1" "${base_name}_r2" "${base_name}_r3")
	
	#######################################################
	############# per replicate: peak calling #############
	#######################################################

	# activate the macs3 conda environment
	conda activate CUTnTag_macs3_env

	echo -e "\n #### [${timestamp}] Active environment: $(basename $CONDA_PREFIX) #### \n"
	echo "MACS3 version: $(macs3 --version)"

	for rep in "${replicates[@]}"; do

		#### change into the rep specific directory ####
		rep_dir=${dir_output}/${rep}
		cd ${rep_dir}

		echo -e "\n ######## [${timestamp}] Starting peak calling for ${rep} ######## \n"
		
		macs3 callpeak -t ${rep}.marked.cleaned.chrsorted.bam \
			-f BAMPE \
			-n ${rep}_macs3 \
			-g mm --keep-dup all \
			--outdir ${rep_dir}/ \
			--nolambda --bdg --SPMR \
			${macs3_args} # these will change depending on the modification/tag
			# can add if needed: --cutoff-analysis (takes much longer)

		sort -k1,1 -k2,2n ${rep}_macs3_peaks.${peak_type} > ${rep}_macs3.sorted.${peak_type}

		# convert to bigwig for IGV visualisation
		sort -k1,1 -k2,2n ${rep}_macs3_treat_pileup.bdg > ${rep}_macs3_treat_pileup.sorted.bdg

		bedGraphToBigWig ${rep}_macs3_treat_pileup.sorted.bdg \
			${mm10_genome} \
			${rep}_macs3.bw

		echo -e "\n [${timestamp}] Finished peak calling for ${rep} \n"


		echo -e "\n [${timestamp}] Finished for ${rep} output found in ${rep_dir} \n"

		rm  ${rep}_macs3_treat_pileup.bdg \
			${rep}_macs3_treat_pileup.sorted.bdg
	done

	conda deactivate

	######################################################
	############## merge CUT&Tag replicates ##############
	######################################################

	## commented out here as alreayd created merged bam and index in script2.3 SEACR peak calling
	# # activate the samtools/bedtools conda environment
	# conda activate CUTnTag_tools_env
	# echo -e "\n #### [${timestamp}] Active environment: $(basename $CONDA_PREFIX) #### \n"

	# echo -e "\n ######## [${timestamp}] Merge replicates and call consensus peaks ######## \n"

	# mkdir -p ${dir_output}/${base_name}
	# cd ${dir_output}/${base_name}

	# # here e.g. replicates[0] refers to 1st item of the replicates array = ${base_name}_r1
	# samtools merge -@ 16 ${base_name}.merged.bam \
	# 	${dir_output}/${replicates[0]}/${replicates[0]}.marked.cleaned.chrsorted.bam \
	# 	${dir_output}/${replicates[1]}/${replicates[1]}.marked.cleaned.chrsorted.bam \
	# 	${dir_output}/${replicates[2]}/${replicates[2]}.marked.cleaned.chrsorted.bam

	# samtools sort -o ${base_name}.merged.chrsorted.bam \
	# 	-T ${base_name}.merged.chrsorted \
	# 	-@ 16 ${base_name}.merged.bam

	# samtools index ${base_name}.merged.chrsorted.bam

	# rm ${base_name}.merged.bam

	# echo -e "\n [${timestamp}] Finished merging replicates for ${base_name} \n"

	# conda deactivate

	#####################################################
	###### call peaks on merged CUT&Tag replicates ######
	#####################################################

	cd ${dir_output}/${base_name}

	# activate the macs3 conda environment
	conda activate CUTnTag_macs3_env
	echo -e "\n #### [${timestamp}] Active environment: $(basename $CONDA_PREFIX) #### \n"

	macs3 callpeak -t ${base_name}.merged.chrsorted.bam \
		-f BAMPE \
		-n ${base_name}_macs3 \
		-g mm --keep-dup all \
		--outdir ${dir_output}/${base_name} \
		--nolambda --bdg --SPMR \
		${macs3_args} # these will change depending on the modification/tag
		# can add if needed: --cutoff-analysis (takes much longer)

	sort -k1,1 -k2,2n ${base_name}_macs3_peaks.${peak_type} > ${base_name}_macs3.sorted.${peak_type}

	# convert to bigwig for IGV visualisation
	sort -k1,1 -k2,2n ${base_name}_macs3_treat_pileup.bdg > ${base_name}_macs3_treat_pileup.sorted.bdg

	bedGraphToBigWig ${base_name}_macs3_treat_pileup.sorted.bdg \
		${mm10_genome} \
		${base_name}_macs3.bw

	echo -e "\n [${timestamp}] Finished calling peaks for ${base_name} \n"

	rm  ${base_name}_macs3_treat_pileup.bdg \
		${base_name}_macs3_treat_pileup.sorted.bdg

	conda deactivate

	######################################################
	############ identify reproducible peaks #############
	######################################################

	# activate the samtools/bedtools conda environment
	conda activate CUTnTag_tools_env
	echo -e "\n #### [${timestamp}] Active environment: $(basename $CONDA_PREFIX) #### \n"

	# create list of all replicate peak files
	peak_files=("${dir_output}/${replicates[0]}/${replicates[0]}_macs3.sorted.${peak_type}" \
				"${dir_output}/${replicates[1]}/${replicates[1]}_macs3.sorted.${peak_type}" \
				"${dir_output}/${replicates[2]}/${replicates[2]}_macs3.sorted.${peak_type}")

	# count peaks between replicates to find overlapping ones
	bedtools multiinter -i ${peak_files[@]} > ${base_name}_macs3.multiinter.bed

	# filter to keep peaks present in at least 2 replicates
	awk '$4>=2 {print $1"\t"$2"\t"$3"\t"$4}' ${base_name}_macs3.multiinter.bed > ${base_name}_macs3.consensus_peaks_raw.bed 

	# merge any overlapping/bookmarked peaks
	bedtools merge -i ${base_name}_macs3.consensus_peaks_raw.bed > ${base_name}_macs3.consensus_peaks.bed

	echo -e "\n [${timestamp}] Finished identifying reproducible peaks (peaks in => 2/3 replicates) for ${base_name} \n"

	echo -e "\n ######## [${timestamp}] ${base_name} completed ######## \n"

	rm  ${base_name}_macs3.multiinter.bed \
		${base_name}_macs3.consensus_peaks_raw.bed

	conda deactivate

done

echo -e "\n ######## [$(timestamp)] Script completed for ${modification} ######## \n"



#!/bin/bash
#SBATCH --job-name=Script2_test_H3K4me1
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script2_test_H3K4me1_WT_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script2_test_H3K4me1_WT_%j.log


#### MUST CHANGE: SCRIPT AND MODIFICATION NAME, CHECK DIRECTORIES & MACS2 ARGS ####


#######################################################
################ set names & parameters ###############
#######################################################

modification=H3K4me1

echo -e "\n######## Script 2 to call peaks for each replicate and merged samples for ${modification} ######## \n"

if [[ ${modification} == "H3K4me1" ]]; then
	AUC="0.02"
elif [[ ${modification} == "H3K4me3" ]]; then
	AUC="0.01"
elif [[ ${modification} == "H3K27ac" ]]; then
	AUC="0.02"
elif [[ ${modification} == "BRG1" ]]; then
	AUC="0.01"
else
    echo -e "\n######## [$(timestamp)] Error: Modification ${modification} is not accepted ######## \n"
    exit 1
fi


dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/alignment_peakcalling/old_${modification}/
mm10_genome=/scratch/prj/id_hill_sims_wellcda/CUTnTag/mm10.genome.size

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# load the Anaconda module
module load anaconda3/2022.10-gcc-13.2.0

# source conda so it works in non-interactive shells
#source $(dirname $(which conda))/../etc/profile.d/conda.sh
eval "$(conda shell.bash hook)"

# activate the SEACR conda environment
conda activate CUTnTag_seacr_env
echo -e "\n #### [$(timestamp)] Active environment: $(basename $CONDA_PREFIX) #### \n"

# print versions of key tools
bedtools --version
samtools --version | head -n 1

# set an error trap
set -o errexit # stop if encounter an error
set -o pipefail # stop on pipeline errors (if any command within a piped command fails)
trap 'echo -e "\n [$(timestamp)] Error in ${BASH_COMMAND}. Exiting. \n" >&2; exit 1' ERR # when a command errors, run the command inside the quotes (print timestamp, write the failed command, send the message to the log file and stop the script)


#######################################################
############ start to loop over conditions ############
#######################################################

# define an array of base names that can be looped over
conditions=("${modification}_WT")


for base_name in "${conditions[@]}"; do

	echo -e "\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ [$(timestamp)] Starting ${base_name} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n"

	# define an array of bio replicate names that can be looped over
	replicates=("${base_name}_r1" "${base_name}_r2")
	
	#######################################################
	############# per replicate: peak calling #############
	#######################################################

	for rep in "${replicates[@]}"; do

		#### change back into the rep specific directory ####
		rep_dir=${dir_output}/${rep}
		cd ${rep_dir}

		echo -e "\n ######## [$(timestamp)] Starting peak calling for ${rep} ######## \n"
		
		# sort be read name so both mates of a paired-end read are next to each other which is needed for bamtobed
		samtools sort -n -o ${rep}.marked.cleaned.qnsorted.bam \
		-T ${rep}.marked.cleaned.qnsorted \
		-@ 16 \
		${rep}.marked.cleaned.chrsorted.bam

		# fix orphaned pairs
		samtools view -b -f 1 -F 0xC ${rep}.marked.cleaned.qnsorted.bam > ${rep}.marked.cleaned.qnsorted.fixed.bam

		samtools sort -n -o ${rep}.marked.cleaned.qnsorted.fixed.sorted.bam \
		-T ${rep}.marked.cleaned.qnsorted.fixed.sorted \
		-@ 16 \
		${rep}.marked.cleaned.qnsorted.fixed.bam

		# for SEACR must convert paired-end BAM files to BED file
		bedtools bamtobed -bedpe -i ${rep}.marked.cleaned.qnsorted.fixed.sorted.bam > ${rep}.bed

		# remove spurious or discordant alignments (following vignette)
		awk '$1==$4 && $6-$2 < 1000 {print $0}' ${rep}.bed > ${rep}.clean.bed

		# extract chr, start and end
		cut -f 1,2,6 ${rep}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${rep}.fragments.bed

		# create fragment coverage bedgraph
		bedtools genomecov -bg -i ${rep}.fragments.bed -g ${mm10_genome} > ${rep}.fragments.bedgraph

		# call SEACR
		bash SEACR_1.3.sh ${rep}.fragments.bedgraph ${AUC} non stringent ${rep}_SEACR


		sort -k1,1 -k2,2n ${rep}_SEACR.stringent.bed > ${rep}_SEACR.stringent.sorted.bed

		echo -e "\n [$(timestamp)] Finished peak calling for ${rep} \n"


		echo -e "\n [$(timestamp)] Finished for ${rep} output found in ${rep_dir} \n"

		#rm  ${rep}.bed \
		#	${rep}.clean.bed \
		#	${rep}.fragments.bed
	done


	######################################################
	############## merge CUT&Tag replicates ##############
	######################################################

	echo -e "\n ######## [$(timestamp)] Merge replicates and call consensus peaks ######## \n"

	mkdir -p ${dir_output}/${base_name}
	cd ${dir_output}/${base_name}

	# here e.g. replicates[0] refers to 1st item of the replicates array = ${base_name}_r1
	samtools merge -@ 16 ${base_name}.merged.bam \
		${dir_output}/${replicates[0]}/${replicates[0]}.marked.cleaned.chrsorted.bam \
		${dir_output}/${replicates[1]}/${replicates[1]}.marked.cleaned.chrsorted.bam 
		#${dir_output}/${replicates[2]}/${replicates[2]}.marked.cleaned.chrsorted.bam

	samtools sort -o ${base_name}.merged.chrsorted.bam \
		-T ${base_name}.merged.chrsorted \
		-@ 16 ${base_name}.merged.bam

	samtools index ${base_name}.merged.chrsorted.bam 

	#rm ${base_name}.merged.bam

	echo -e "\n [$(timestamp)] Finished merging replicates for ${base_name} \n"


	#####################################################
	###### call peaks on merged CUT&Tag replicates ######
	#####################################################

	# sort be read name so both mates of a paired-end read are next to each other which is needed for bamtobed
	samtools sort -n -o ${base_name}.merged.qnsorted.bam \
	-T ${base_name}.merged.qnsorted \
	-@ 16 \
	${base_name}.merged.chrsorted.bam

	# for SEACR must convert paired-end BAM files to BED file
	bedtools bamtobed -bedpe -i ${base_name}.merged.qnsorted.bam  > ${base_name}.bed

	# remove spurious or discordant alignments (following vignette)
	awk '$1==$4 && $6-$2 < 1000 {print $0}' ${base_name}.bed > ${base_name}.clean.bed

	# extract chr, start and end
	cut -f 1,2,6 ${base_name}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${base_name}.fragments.bed

	# create fragment coverage bedgraph
	bedtools genomecov -bg -i ${base_name}.fragments.bed -g ${mm10_genome} > ${base_name}.fragments.bedgraph

	# call SEACR
	bash SEACR_1.3.sh ${base_name}.fragments.bedgraph ${AUC} non stringent ${base_name}_SEACR_merged


	sort -k1,1 -k2,2n ${base_name}_SEACR_merged.stringent.bed > ${base_name}_SEACR_merged.stringent.sorted.bed

	echo -e "\n [$(timestamp)] Finished calling peaks for merged ${base_name} \n"

	#rm ${base_name}.bed \
	#	${base_name}.clean.bed \
	#	${base_name}.fragments.bed


	######################################################
	############ identify reproducible peaks #############
	######################################################

	# create list of all replicate peak files
	peak_files=("${dir_output}/${replicates[0]}/${replicates[0]}_SEACR.stringent.sorted.bed" \
				"${dir_output}/${replicates[1]}/${replicates[1]}_SEACR.stringent.sorted.bed") 
				#"${dir_output}/${replicates[2]}/${replicates[2]}_SEACR.stringent.sorted.bed")

	# count peaks between replicates to find overlapping ones
	bedtools multiinter -i ${peak_files[@]} > ${base_name}.multiinter.bed

	# filter to keep peaks present in at least 2 replicates
	awk '$4>=2 {print $1"\t"$2"\t"$3"\t"$4}' ${base_name}.multiinter.bed > ${base_name}.consensus_peaks_raw.bed 

	# merge any overlapping/bookmarkes peaks
	bedtools merge -i ${base_name}.consensus_peaks_raw.bed > ${base_name}.consensus_peaks.bed 

	echo -e "\n [$(timestamp)] Finished identifying reproducible peaks (peaks in => 2/3 replicates) for ${base_name} \n"


	echo -e "\n ######## [$(timestamp)] ${base_name} completed ######## \n"

done

conda deactivate

echo -e "\n ######## [$(timestamp)] Script completed for ${modification} ######## \n"



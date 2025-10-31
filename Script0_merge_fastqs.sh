#!/bin/bash
#SBATCH --job-name=Script0_merge_fastqs
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script0_merge_fastqs_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script0_merge_fastqs_%j.log

echo "Script0 to merge fastqs for every CUT&Tag sample"

dir_input=/scratch/prj/id_hill_sims_wellcda/CUTnTag/fastqs
dir_output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/merged_fastqs
mkdir -p ${dir_output}

modifications=("H3K4me1" "H3K4me3" "H3K27ac" "BRG1" "RBP1s5p")
conditions=("UI" "IL10" "WT" "SteE")
replicates=("r1" "r2" "r3")
reads=("R1" "R2")

for mod in "${modifications[@]}"; do
	for cond in "${conditions[@]}"; do
		for rep in "${replicates[@]}"; do
			for r in "${reads[@]}"; do

				# define the file names
				file_base=${dir_input}/${mod}_${cond}_${rep}

				files=( ${file_base}_S*_L001_${r}_001.fastq.gz )
				file1=${files[0]}

				files2=( ${file_base}_S*_L002_${r}_001.fastq.gz )
				file2=${files2[0]}
				
				# extract the sample number
				sample=$(basename "$file1" | sed -E 's/.*(S[0-9]+).*/\1/')

				# define the new file name
				new_file=${dir_output}/${mod}_${cond}_${rep}_${sample}_${r}_merged.fastq.gz

				cat "${file1}" "${file2}" > "${new_file}"

				echo "merge ${file1} and ${file2}"
				echo "output is ${new_file}"

			done
		done
	done
done


echo "##### Finished merging fastqs #####"

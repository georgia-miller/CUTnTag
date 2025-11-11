#!/bin/bash
#SBATCH --job-name=Script3_define_promoters_and_enhancers
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script3_%j.log
#SBATCH --error=/scratch/prj/id_hill_sims_wellcda/CUTnTag/logs/Script3_%j.log


# take the differentially accessible peaks between SteE and WT STm infected iBMDMs and using the consensus peaks for CUT&Tag, define them as promoters (H3K4me3) or enhancers (H3K4me1 only)

echo -e "\n######## Script 3 to define promoters and enhancers for SteE diff accessible peaks ######## \n"

base_dir=/scratch/prj/id_hill_sims_wellcda/CUTnTag

ATAC_peaks=${base_dir}/ATAC_seq_top300_steE_vs_WT_peaks_iBMDMs.bed
H3K4me1_WT=${base_dir}/alignment_peakcalling/H3K4me1/H3K4me1_WT/H3K4me1_WT.consensus_peaks.bed
H3K4me3_WT=${base_dir}/alignment_peakcalling/H3K4me3/H3K4me3_WT/H3K4me3_WT.consensus_peaks.bed

# define timestamp to use for logging messages
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

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
########### compare ATAC with CUT&Tag peaks ###########
#######################################################

# look for overlap of the peaks with each modification
bedtools multiinter -i ${ATAC_peaks} ${H3K4me1_WT} ${H3K4me3_WT} \
	-header -names ATAC H3K4me1 H3K4me3 > SteE_vs_WT_peaks_with_mods.bed

# file outputted has columns: chrom, start, end, num, list, ATAC, H3K4me1, H3K4me3

# promoters have H3K4me3, filter for those that have ATAC=1 and H3K4me3=1
# this line: keep header, filter to rows where cols 6 and 7 = 1 and output it
awk 'NR==1 || ($6==1 && $7==1)' SteE_vs_WT_peaks_with_mods.bed > SteE_vs_WT_peaks_promoters.bed

# enhancers have H3K4me1 peaks but not H3K4me3, filter for those that have ATAC=1 and H3K4me1=1, H3K4me3=0
awk 'NR==1 || ($6==1 && $7==0 && $8==1)' SteE_vs_WT_peaks_with_mods.bed > SteE_vs_WT_peaks_enhancers.bed


conda deactivate

echo -e "\n ######## [$(timestamp)] Script completed ######## \n"





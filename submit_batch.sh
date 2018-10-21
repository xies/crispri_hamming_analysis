#!/bin/bash
set -e

#Manage the batch SLURM submission for pairwise Hamming distance analysis

begin=1
total_seqs=16226
Nseqperiter=1000

JOB=pairwise_hamming

MAIL_USER=xies@stanford.edu

INPUT_DIR='/home/xies/data/crispri_hamming'

function submit() {

	GUIDE_NAME=${1}
	INPUT_FILE=${INPUT_DIR}/${GUIDE_NAME}
	ITER_ID=${2}

	OUT_DIR=${INPUT_DIR}/${INPUT_PULSES}.iter_${ITER_ID}
	OUT_NAME=${OUT_DIR}.csv
	STD_OUT=${OUT_DIR}.stdout
	STD_ERR=${OUT_DIR}.stderr
	SCRIPT=${OUT_DIR}.qsub

	CMD="/home/xies/Code/crispri_hamming_distance/pairwise_hamming_distance.py ${INPUT_FILE} ${ITER_ID} ${Nseqperiter} ${OUT_NAME}"

	if [ -e ${SCRIPT} ]
	then
		return
	fi

	echo "Writing script to ${SCRIPT}"
	cat << EOF > ${SCRIPT}
#!/bin/bash
#SBATCH --account=skotheim
#SBATCH --job-name=${JOB}
#SBATCH --output=${STD_OUT}
#SBATCH --error=${STD_ERR}
#SBATCH --time=04:00:00
#SBATCH -c 12
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=${MAIL_USER}
#$ -
${CMD}
EOF

	echo "Submitting: ${INPUT_FILE} #$ITER_ID to be analyzed.."
	sbatch ${SCRIPT}
}

for (( i=$begin; i<=$total_seqs; i=i+$Nseqperiter ))
do
	echo $i
	submit sgRNA_unique.fasta ${i}
done

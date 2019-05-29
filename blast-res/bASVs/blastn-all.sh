#!/bin/bash
#SBATCH -J blastn-01      # job name for array
#SBATCH -n 1                  # Number of cores
#SBATCH -N 1                   # Ensure that all cores are on one machine
#SBATCH -t 0-48:00             # Runtime in D-HH:MM
#SBATCH -p serial_requeue             # Partition to submit to
#SBATCH --mem=64000            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o blastn-01.out      # File to which STDOUT will be written
#SBATCH -e blastn-01.err      # File to which STDERR will be written
#SBATCH --mail-type=END              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=phumphrey@g.harvard.edu  # Email to which notifications will be sent

module load blast
export BLASTDB='/n/regal/informatics_public/ref/ncbi/nt'
cd /n/desai_lab/users/phumphrey/phyllosphere/phy/blast-res
#blastn -query "${1}" -db nt -out "${1}".csv -num_alignments 10 -outfmt '10 qseqid sseqid sgi sseq pident staxid sblastname sblastnames'
blastn -query "${1}" -db nt -out "${1}".csv -num_alignments 10 -outfmt '10 qseqid sseqid sgi sseq pident staxid'

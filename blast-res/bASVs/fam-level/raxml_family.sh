#!/bin/bash
#SBATCH -J raxml-fam      # job name for array
#SBATCH -n 4                  # Number of cores
#SBATCH -N 1                   # Ensure that all cores are on one machine
#SBATCH -t 0-48:00             # Runtime in D-HH:MM
#SBATCH -p serial_requeue             # Partition to submit to
#SBATCH --mem=64000            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o raxml-fam.out      # File to which STDOUT will be written
#SBATCH -e raxml-fam.err      # File to which STDERR will be written
#SBATCH --mail-type=END              # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=phumphrey@g.harvard.edu  # Email to which notifications will be sent

module load gcc/7.1.0-fasrc01 openmpi/2.1.0-fasrc02 RAxML/8.2.11-fasrc02
cd /n/desai_lab/users/phumphrey/phyllosphere/phy/blast-res/fam-level
raxmlHPC-HYBRID-SSE3 -s "${1}" -f a -m GTRCAT -w /n/desai_lab/users/phumphrey/phyllosphere/phy/blast-res/fam-level/ -N 100 -n "${1}_run" -p 12345 -T 4 -x 12345

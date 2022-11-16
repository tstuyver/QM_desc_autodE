#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --time=99:59:00
#SBATCH --constraint=centos7
#SBATCH --partition=sched_mit_ccoley
#SBATCH --nodelist node1238
#SBATCH --mem 192000
#SBATCH --output=dipoles

source /home/tstuyver/.bashrc
conda activate descriptors
module load gaussian/16.c01
GAUSS_SCRDIR=/nobackup1/tstuyver//$SLURM_JOB_NAME-$SLURM_JOB_ID
mkdir -p $GAUSS_SCRDIR
chmod 750 $GAUSS_SCRDIR

python main.py --ismiles dipoles_test.csv --DFT-folder dipoles_DFT --xtb-folder dipoles_XTB_opt

rm -rf $GAUSS_SCRDIR

#!/bin/bash 
#SBATCH -N 1 
#SBATCH -n 48
#SBATCH --time=50:59:00 
#SBATCH --mem 192000 
#SBATCH --output=6416

source /home/gridsan/tstuyver/.bashrc 
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate QM_descriptors
export g16root=/home/gridsan/groups/RMG/Software/gaussian/
export PATH=$g16root/g16/:$g16root/gv:$PATH
GAUSS_SCRDIR=~/scratch/gaussian//$SLURM_JOB_NAME-$SLURM_JOB_ID 
. $g16root/g16/bsd/g16.profile
mkdir -p $GAUSS_SCRDIR   
chmod 750 $GAUSS_SCRDIR  
 
python main.py --ismiles scope6416smiles.csv --xtb_folder xtb_6416 --DFT_folder DFT_6416

echo "jobs done, removing scratch folder: $(date)"
rm -rf $GAUSS_SCRDIR


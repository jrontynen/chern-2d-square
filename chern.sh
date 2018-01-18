#!/bin/sh
#SBATCH --time=0-04:00:00 --mem-per-cpu=5000
#SBATCH -o out/jobc1-%a.out
#SBATCH --array=1-300

cp run_chern.m chern.m cheval.m Struve* $TMPDIR
rm ./out/outc1-$SLURM_ARRAY_TASK_ID.mat
cd $TMPDIR
mkdir $TMPDIR/out
trap "mv -f $TMPDIR/out/* $WRKDIR/matlab/pr2/out/; exit" TERM EXIT

module load matlab
matlab -nojvm -r "run_chern($SLURM_ARRAY_TASK_ID,300)"
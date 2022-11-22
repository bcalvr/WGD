jobname1=`echo $1 | cut -f 1 -d '.'`
split_num=$(wc -l < dateindex_splits.csv)

sbatch  << EOF
#!/bin/bash
#SBATCH --job-name=${jobname1}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200gb
#SBATCH --time=240:00:00
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=4
#SBATCH --output=${jobname1}_%a_err.log
#SBATCH --array=1-${split_num}
cd \$SLURM_SUBMIT_DIR
module load python/3.8.3
row=\$SLURM_ARRAY_TASK_ID
read sidx eidx <<< \$(awk -F',' 'NR=='\$row'{print \$1} NR=='\$row'{print \$2}' dateindex_splits.csv )
echo \$sidx \$eidx 
python $1 \$sidx \$eidx

EOF

jobname=`echo $1 | cut -f 1 -d '.'`

echo "#!/bin/csh"
echo "#SBATCH --job-name=${jobname}"
echo '#SBATCH --nodes=1'
echo '#SBATCH --ntasks=3'
echo '#SBATCH --partition=highmem'
echo '#SBATCH --mem=500GB'
echo '#SBATCH --time=360:00:00'
#echo '#SBATCH --cpus-per-task=4'
echo "#SBATCH --output=${jobname}_err.log"
echo
echo 'cd $SLURM_SUBMIT_DIR'
echo 'module load python/3.8.3'
echo "python jh_cncb_update_pre_processing.py"
echo "python aggregated_mutations0_01_16_22.py & sequence_country_counts.py"
echo "python index_splitter.py"
echo "split_num=$(wc -l < dateindex_splits.csv)"
echo "read sid eid <<< \$(awk -F',' 'NR==1{print \$1} NR=='\$split_num'{print \$2}' dateindex_splits.csv )"
echo "python jh_cncb_update_processing_I.py $sid $eid"
echo "srun -n 1 ./looper.sh IR_calc.py & srun -n 1 ./looper.sh FR_calc.py & srun -n 1 ./looper.sh PP_calc.py"
echo "python jh_cncb_update_processing_II.py"

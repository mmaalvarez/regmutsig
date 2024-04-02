for genome_half in 1st_half 2nd_half null
do
	sbatch run_s.sh $genome_half
done

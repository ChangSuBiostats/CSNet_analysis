suffix=$1
setting=$2
cd ..
fn=results/tables/table_real_data_"$suffix"_"$setting".txt
echo '' > $fn
for n in 150 600
do
	for p in 100 200
	do
	Rscript print_tables.R \
	--result_prefix results/real_data_"$suffix"/K_2/log_var_8.0_equal_strength_TRUE/n_"$n"_p_"$p" \
	--setting "$setting" >> $fn
done
done

echo $fn

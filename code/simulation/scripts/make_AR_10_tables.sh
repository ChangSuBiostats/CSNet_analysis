suffix=$1
setting=$2

cd ..
fn=results/tables/table_AR_10_"$suffix"_"$setting".txt
echo '' > $fn
for n in 150 600
do
	for p in 100 200
	do
	Rscript print_tables.R \
	--result_prefix results/AR_10_rho_0.8_0.8/K_2/log_var_8.0_equal_strength_TRUE/n_"$n"_p_"$p" \
	--setting "$setting" >> $fn
done
done

echo $fn

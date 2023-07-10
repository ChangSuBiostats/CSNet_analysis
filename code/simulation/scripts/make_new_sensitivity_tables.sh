setting=$1
cor_model=$2
cd ..
fn=results/tables/table_new_sensitivity_"$cor_model"_"$setting".txt
echo '' > $fn
for kappa in 0.00 0.60 0.90
do
	for b in -0.20 -0.10 0.00 0.10 0.20 
	do
	Rscript print_tables.R \
	--result_prefix results/sensitivity/"$cor_model"/K_2/log_var_8.0_equal_strength_TRUE/n_600_p_100/kappa_"$kappa"_b_"$b" \
	--setting sensitivity_"$setting" >> $fn
done
done

echo '\\\hline' >> $fn

echo $fn

setting=$1
equal_strength=$2
log_var=$3
cd ..
fn=results/tables/table_K_4_"$setting"_log_var_"$log_var"_equal_strength_"$equal_strength".txt
echo '' > $fn


n=600
p=100
Rscript print_tables.R \
--result_prefix results/AR_10/K_4/log_var_"$log_var"_equal_strength_"$equal_strength"/n_"$n"_p_"$p" \
--setting "$setting" \
--metrics error \
>> $fn

echo '\\\hline' >> $fn


Rscript print_tables.R \
--result_prefix results/AR_10/K_4/log_var_"$log_var"_equal_strength_"$equal_strength"/n_"$n"_p_"$p" \
--setting "$setting" \
--metrics selection \
>> $fn

echo '\\\hline' >> $fn

echo $fn

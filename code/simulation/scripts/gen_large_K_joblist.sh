prefix=K_$1
suffix=equal_strength_$2
log_var=$3

fn="$prefix"_"$suffix"_log_var_"$log_var"_joblist.txt
echo '' > $fn
for n in 600
do
  for p in 100
  do
    for i_rep in {1..200}
    do
      echo "cd ..; ml miniconda; conda activate r_csnet_3; Rscript experiment.R \
--n_rep 200 \
--i_rep "$i_rep" \
--n "$n" \
--p "$p" \
--K "$1" \
--log_var "$log_var" \
--cor_model AR_10 \
--beta1 2 \
--beta2 1 \
--rho1 0.8 \
--rho2 0.8 \
--equal_strength "$2" \
--save_est T" >> $fn
    done
  done
done

echo $fn

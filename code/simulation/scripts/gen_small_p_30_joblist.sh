prefix=small_p
suffix=''

fn="$prefix"_"$suffix"_joblist.txt
echo '' > $fn
for n in 150 600
do
  for p in 30
  do
    for i_rep in {1..200}
    do
      echo "cd ..; ml miniconda; conda activate r_csnet_3; Rscript experiment.R \
--n_rep 200 \
--i_rep "$i_rep" \
--n "$n" \
--p "$p" \
--K 2 \
--log_var 8.0 \
--cor_model AR_10 \
--beta1 2 \
--beta2 1 \
--rho1 0.8 \
--rho2 0.8 \
--save_est T" >> $fn
    done
  done
done

echo $fn

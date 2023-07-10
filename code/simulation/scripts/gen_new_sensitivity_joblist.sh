prefix=$1
suffix=new_sensitivity
fn="$prefix"_"$suffix"_joblist.txt
echo '' > $fn
for kappa in 0.0 0.6 0.9
do
  for b in -0.2 -0.1 0.0 0.1 0.2
  do
    for i_rep in {1..200}
    do
      echo "cd ..; ml miniconda; conda activate r_csnet_3; Rscript experiment.R \
--n_rep 200 \
--i_rep "$i_rep" \
--n 600 \
--p 100 \
--K 2 \
--log_var 8.0 \
--cor_model "$prefix" \
--beta1 2 \
--beta2 1 \
--rho1 0.8 \
--rho2 0.8 \
--sensitivity T \
--kappa "$kappa" \
--b "$b" \
--save_est T" >> $fn
    done
  done
done

echo $fn

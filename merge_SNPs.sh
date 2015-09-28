#! /bin/bash

##INPUT

#folder with file for each pop of SNPs
folder="/home/marco/pool_first/SNP_database/p53"
#names of pops
fold_pops=(07C_snape 07D_snape 07E_snape 07F_snape 07G_snape 07J_snape 07K_snape 07L_snape 07M_snape 07N_snape 07O_snape 07P_snape 07Q_snape 07R_snape 11AA_snape 11AB_snape 11AC_snape 11AE_snape 11AG_snape 11AH_snape 11AJ_snape 11B_snape 11G_snape 11H_snape 11J_snape 11K_snape 11L_F0_snape 11M_snape 11N_snape 11O_snape 11P_snape 11Q_snape 11S_snape 11T_snape 11U_snape 11W_snape 11X_snape 11Z_snape 14A_snape 14B_snape 14C_snape 14D_snape 14E_snape 11D_snape 11V_snape 11R_snape 11A_snape 11C_snape 11E_snape 11F_snape 11Y_snape 07H_snape 13ARE_snape)
#rownames of SNPs
righe="/home/marco/pool_first/SNP_database/p53/SNP_tot_snape"
#output name
nome_out="/home/marco/pool_first/SNP_database/p53_tab_snape"

min_pop=2
min_freq=0.05

################

#ulimit -v 7500000 -m 7500000 

#lascia cosi' di meno da problemi x trial_min.R
ulimit -v 3500000 -m 3500000 

rm $nome_out"_pops"

i=0
len=${#fold_pops[*]}
while [ $i -lt $len ]; do

  echo ${fold_pops[$i]} >> $nome_out"_pops"
  let i++
  
done

#ulimit -a > $nome_out"_ulimit"


echo $folder $righe $nome_out
Rscript merge_SNPs_pops.R $folder $righe $nome_out

# righe é il rownames nome_out_col é colnames

rm $nome_out"_pops"


#Selection of the SNP, remove triallelic.

Rscript trial_min.R $nome_out $min_pop $min_freq $righe $nome_out"_clean"




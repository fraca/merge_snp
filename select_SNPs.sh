#! /bin/bash

##INPUT


fold_pops=(07C 07D 07E 07F 07G 07J 07K 07L 07M 07N 07O 07P 07Q 07R 11C 11AA 11AB 11AC 11AE 11AG 11AH 11AJ 11B 11G 11H 11J 11K 11L_F0 11M 11N 11O 11P 11Q 11S 11T 11U 11W 11X 11Z 11D 11V 11R 11A 11E 11F 11Y 14A 14B 14C 14D 14E 07H 13ARE)

folder="/home/marco/server_vb/bam_ALL/pop"


dir_out=/home/marco/pool_first/SNP_database/p53

int=/home/marco/pool_first/SNP_database/p53_multi_all.BED

################

rm -r $dir_out
mkdir $dir_out




i=0
len=${#fold_pops[*]}
while [ $i -lt $len ]; do
  echo "$i: ${fold_pops[$i]}"
  Rscript SNP2BED.R $folder${fold_pops[$i]}"/p"${fold_pops[$i]}"_SNP_scaf/" $dir_out"/"${fold_pops[$i]}"_med_"
  
  bedtools intersect -u -a $dir_out"/"${fold_pops[$i]}"_med_varscan.BED" -b $int -f 1 | awk ' $4 != "N" {print $1"-"$2"-"$4"\t"$5"\t"$6}' > $dir_out"/"${fold_pops[$i]}"_varscan" 
  rm $dir_out"/"${fold_pops[$i]}"_med_varscan.BED"
  
  cut -f 1 $dir_out"/"${fold_pops[$i]}"_varscan" >> $dir_out"/SNP_tot_varscan_t"

  bedtools intersect -u -a $dir_out"/"${fold_pops[$i]}"_med_snape.BED" -b $int -f 1 | awk ' $4 != "N" {print $1"-"$2"-"$4"\t"$5"\t"$6}'  > $dir_out"/"${fold_pops[$i]}"_snape" 
  #rm $dir_out"/"${fold_pops[$i]}"_med_snape.BED"
  
  cut -f 1 $dir_out"/"${fold_pops[$i]}"_snape" >> $dir_out"/SNP_tot_snape_t"
  
let i++
done


sort -u $dir_out"/SNP_tot_varscan_t" > $dir_out"/SNP_tot_varscan"
rm $dir_out"/SNP_tot_varscan_t"

sort -u $dir_out"/SNP_tot_snape_t" > $dir_out"/SNP_tot_snape"
rm $dir_out"/SNP_tot_snape_t"


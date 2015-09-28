#! /bin/bash

##merge different positions
#used multiIntersectBed (from bedtools 2.22 copy /home/marco/program/bedtools2/bin to /usr/bin)
##INPUT

fold_pops=(07C 07D 07E 07F 07G 07J 07K 07L 07M 07N 07O 07P 07Q 07R 11AA 11AB 11AC 11AE 11AG 11AH 11AJ 11B 11G 11H 11J 11K 11L_F0 11M 11N 11O 11P 11Q 11S 11T 11U 11W 11X 11Z 14A 14B 14C 14D 14E 11D 11V 11R 11A 11C 11E 11F 11Y 07H 13ARE)

folder="/home/marco/server_vb/bam_ALL/pop"

nome_out="p53"




################


scaf_tot=(scaffold_1 scaffold_2 scaffold_3 scaffold_4 scaffold_5 scaffold_6 scaffold_7 scaffold_8) 


unicum_m=""
i=0
len=${#fold_pops[*]}

len2=${#scaf_tot[*]}

while [ $i -lt $len ]; do
  #echo "$i: ${fold_pops[$i]}"
  j=0
  while [ $j -lt $len2 ]; do
    #echo "$i $j"
    unicum_m=$unicum_m$folder${fold_pops[$i]}"/p"${fold_pops[$i]}"_SNP_scaf/"${scaf_tot[$j]}"_multi.BED "
    let j++
  done
  let i++
done
#unicum=${unicum%?} X togliere ultimo char

#echo $unicum_m

multiIntersectBed -i $unicum_m > $nome_out"_multi.BED"
awk '$4=='$len' {print $1"\t"$2"\t"$3}' $nome_out"_multi.BED" > $nome_out"_multi_all.BED"

rm $nome_out"_multi.BED" ##big file to remove.


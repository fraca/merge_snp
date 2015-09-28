#### crea file freq con intestazione like BED file
###INPUT
#folder_in contain each chromosome
#scafs name of each chromosome
#folder_out pathway and output name
###OUTPUT
#folder_out_snape.BED SNP frequencies calculated with snape
#folder_out_varscan.BED SNP frequencies calculated with varscan
#chromosomes position position ref_allele alt_allele alt_allele_freq
# scaffold_1  141  141  C	T	1
# scaffold_1	351	351	C	A	0.424242424242424
# scaffold_1	845	845	T	G	0.0416666666666667
# scaffold_1	849	849	T	C	1
# scaffold_1	1303	1303	C	T	0.024
# scaffold_1	1398	1398	A	G	0.227513227513228
# scaffold_1	1593	1593	C	T	0.0202020202020202
# scaffold_1	1602	1602	C	T	0.370786516853933
# scaffold_1	1623	1623	C	T	0.0236686390532544
# scaffold_1	1825	1825	A	T	0.402173913043478


SNP2BED=function(folder_in,scafs,folder_out) {
  ##SNAPE 
  #folder_in="/home/marco/server_vb/bam_ALL/pop07C/p07C_SNP_scaf/"
  #to remove ref in snape alts
  rem_ref=function(due) {
    return(sub(due[1],"",due[2]))
  }
  tab_tot=matrix(ncol=6)
  colnames(tab_tot)=c("chr","start","end","ref","alts","f_alt")
  for(i in 1:length(scafs)) {
    #cat("snape",scafs[i],"\n")
    tab=read.table(paste(folder_in,scafs[i],".snape",sep=""),stringsAsFactors=FALSE)
    tab[,5]=tab[,5]/(tab[,4]+tab[,5])
    tab[,4]=tab[,8]
    tab=tab[,c(1,2,2,3,4,5)]
    colnames(tab)=c("chr","start","end","ref","alts","f_alt")
    tab_tot=rbind(tab_tot,tab)
  }
  tab_tot=tab_tot[-1,]
  aa=apply(tab_tot[,c("ref","alts")],1,rem_ref)
  tab_tot[,"alts"]=aa
  ##change for problem with Snape if fixed. take only the first nucleotide in alt
  sel=which(lapply(aa,nchar)!=1)
  rem=NULL

  for(i in 1:length(sel)) {
    #cat(tab_tot[sel[i],"alts"],"\n")
    if(tab_tot[sel[i],"f_alt"]==1) {      
      tab_tot[sel[i],"alts"]=substr(tab_tot[sel[i],"alts"],1,1)
    } else {
      rem=c(rem,sel[i])
    }
  }
  if(!is.null(rem)) {
    tab_tot=tab_tot[-rem,]  
  }
  
  write.table(tab_tot,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste(folder_out,"snape.BED",sep=""))        
  
  
  
  ##VARSCAN
  
#   folder_in="/home/marco/server_vb/bam_ALL/pop07C/p07C_SNP_scaf/" 
#   folder_out="prova/bla1"
#   scafs="scaffold_1"
  
  tab_tot=matrix(ncol=6)
  colnames(tab_tot)=c("chr","start","end","ref","alts","f_alt")
  
  for(i in 1:length(scafs)) {
    #cat("varscan",scafs[i],"\n")
    tab=read.table(paste(folder_in,scafs[i],".varscan",sep=""),stringsAsFactors=FALSE,header = TRUE)  
    tab[,"VarFreq"]=gsub("%","",tab[,"VarFreq"])
    mode(tab[,"VarFreq"])="numeric"
    tab[,"VarFreq"]=tab[,"VarFreq"]/100
    tab=tab[,c("Chrom", "Position", "Position","Ref","VarAllele","VarFreq")]
    colnames(tab)=c("chr","start","end","ref","alts","f_alt")
    tab_tot=rbind(tab_tot,tab)
    
  }
  #remove triallelic SNP called by VarScan
  tab_tot=tab_tot[-1,]
  dup=!duplicated(tab_tot$start)
  tab_tot=tab_tot[dup,]
  
  write.table(tab_tot,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste(folder_out,"varscan.BED",sep=""))  
}

#######################################################################################


args <- commandArgs(trailingOnly = TRUE)

folder_in=args[1]
#scafs=args[2]
folder_out=args[2]
scafs=c("scaffold_1", "scaffold_2", "scaffold_3", "scaffold_4", "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8")
SNP2BED(folder_in,scafs,folder_out)
#cat(folder_in,scafs,folder_out,"\n")


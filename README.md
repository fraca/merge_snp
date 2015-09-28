merge_snp
==============

Several script to merge the SNP data from different populations called with the pipeline PWGS.

Software used:

- bedtools (BEDTools, Quinlan et al. 2010)
- R (R Core Team 2015)
- bwa mem (BWA, Li et al. 2013)
- samtools (SAMtools, Li et al. 2009)
- SortSam.jar (Picard tools, http://picard.sourceforge.net)
- MarkDuplicates.jar (Picard tools, http://picard.sourceforge.net)
- identify-genomic-indel-regions.pl (PoPoolation, Kofler et al. 2012)
- filter-pileup-by-gtf.pl (PoPoolation, Kofler et al. 2012)
- snape-pooled (Snape, Raineri et al. 2012)
- VarScan (VarScan, Koboldt et al. 2012)
- NPStat (NPStat, Ferretti et al. 2013)

##merge_positions.sh  

Merge all positions sequenced in all populations. it use the _multi.BED files.


##select_SNPs.sh  

create folder with the SNPs tab for each population analysed with Snape and VarScan.


##merge_SNPs.sh  

merge all the SNPs files.

1. row file (p53_tab_snape)  
file with all the SNPs that are in regions sequenced in all populations. First the ALT allele after the frequencies for each population. Value with NA means that there is the REF allele and is fixed.

2. clean file (p53_tab_snape_clean.BED)  
file in bedtools format scaffold position position REF allele ALT allele frequencies for all the populations. Criteria for selection:
- only biallelic SNPs 
- polymorphic in at least two populations
- at least one population has the frequency more than 0.05.

3. rfixed (p53_tab_snape_rfixed)  
file with SNPs different from reference but fixed in all populations.

4. rmin_freq (p53_tab_snape_rmin_freq)  
file with the discharged SNPs that have frequency less than 0.05 in all populations.

5. rmin_pop (p53_tab_snape_rmin_pop)  
file with the discharged SNPs that are polymorphic in one population.

6. rtri (p53_tab_snape_rtri)  
file with the discharged triallelic SNPs.





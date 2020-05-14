setwd("~/Desktop/Sarah/Salmon/TransAtlanticPaper/Recombination/")
library(genepopedit)
library(windowscanr)
library(ggplot2)
library(cowplot)

#This R script includes
  # 1 - Plots for recombination rate across chromosomes (note any misplaced SNPs were removed to calculate recombination rate)
      #outlier regions were determined as values that were +/- 2 SD from the mean
  # 2 - Gene density (or count) per 1MBP were calculated using windowscanr
        #outlier regions were determined as values that were +/- 2 SD from the mean
  # 3 - Plots for gene density were made and combined with those for recombination using cowplot
  
#Outlier script (taken from Brenna Forester's RDA vignette)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#Read in recombination rate for Canada and Europe
Canada_data=read.csv("RecombinationRate_allSNPs_Chr_Canada.csv", header=T, fill=T, stringsAsFactors = T)
EU_data=read.csv("RecombinationRate_allSNPs_Chr_Europe.csv", header=T, fill=T, stringsAsFactors = T)

#Note - plotting for each chromosome was done seperately using recombination information
#Kept in potential misplaced snps for the CM x BP plot
#but misplaced snps were removed to show recombination rate across the chromosome (smoothed line plot)
#Note that recombination rate is in cm/Mbp

############ SSA 06
#Europe data
EU_data_Chr6=EU_data[which(EU_data$Chromosome==6),]

#remove 2 misplaced SNPs
EU_data_Chr6_2removed=EU_data_Chr6[which(EU_data_Chr6$SNP_ID!="ESTNV_34441_695" & EU_data_Chr6$SNP_ID!="ESTNV_24999_794" ),]

#Canada
ssa06_canada_bp=Canada_data[which(Canada_data$Chromosome==6),]
#remove extra snps (misplaced)
ssa06_canada_bp_removed=ssa06_canada_bp[which(ssa06_canada_bp$SNP_ID!= "ESTNV_24999_794" & 
                                                ssa06_canada_bp$SNP_ID!= "ESTNV_25601_615" &  
                                                ssa06_canada_bp$SNP_ID!= "ESTNV_30516_1089"  ),]

Merged_eu_candaa_chr06_Recom=ggplot(data=ssa06_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.0e7, ymax=0, ymin=125, xmin=7.3e7, fill="gray")+
  geom_point(data=ssa06_canada_bp,aes(x=start, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr6,aes(x=start, y=Female.map, col="Europe"))+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+ theme_bw()+
  ylab("CM")+xlab("Ssa06 Position")+  theme(panel.grid = element_blank(), legend.position = "top")

#Check number of outliers for Ssa06 Europe and North AM
outliers(ssa06_canada_bp_removed$CM_BP,2)
#North AM - just 0 cM/mbp is outlier

outliers(EU_data_Chr6_2removed$CM_BP,2)
#Europe - just 51.42029 cM/mbp is outlier


#Merged plot - canada europe cm/Mbp
merged_ssa06_CMbp=ggplot(data=ssa06_canada_bp_removed)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.0e7, ymax=5, ymin=0, xmin=7.3e7, fill="gray")+
  geom_smooth(data=ssa06_canada_bp_removed,
              aes(x=start, y=CM_BP, col="North Am"), stat="identity")+
  geom_smooth(data=EU_data_Chr6_2removed,
              aes(x=start, y=CM_BP, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa06 Position")+theme_bw()+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,5), expand = F)+
  NULL

########################################################
################ Ssa13
########################################################

EU_data_Chr13=EU_data[which(EU_data$Chromosome==13),]
ssa13_canada_bp=Canada_data[which(Canada_data$Chromosome==13),]

#remove misplaced SNPs
ssa13_canada_bp_removed=ssa13_canada_bp[which(ssa13_canada_bp$SNP_ID!= "ESTNV_32655_334" ),]


Merged_eu_candaa_chr13_Recom=ggplot(data=ssa13_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.5e7, ymax=0, ymin=125, xmin=7.8e7, fill="gray")+
  geom_point(data=ssa13_canada_bp,aes(x=start, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr13,aes(x=start, y=Female.map, col="Europe"))+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+theme_bw()+
  ylab("CM")+xlab("Ssa13 Position")+  theme(panel.grid = element_blank(), legend.position = "top")

 
#Merged plot _ canada europe -cm/Mbp

merged_ssa13_CMbp=ggplot(data=ssa13_canada_bp_removed)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.5e7, ymax=5, ymin=0, xmin=7.8e7, fill="gray")+
  geom_smooth(data=ssa13_canada_bp_removed,
              aes(x=start, y=CM_BP, col="North Am"), stat="identity")+
  geom_smooth(data=EU_data_Chr13,
              aes(x=start, y=CM_BP, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa13 Position")+theme_bw()+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,3), expand = F)+
  NULL


#Check number of outliers for Ssa13 Europe and North AM
outliers(ssa13_canada_bp_removed$CM_BP,2)

#North AM - just 0 cM/mbp is outlier (two values)
outliers(EU_data_Chr13$CM_BP,2)

#Europe - just 0 cM/mbp is outlier (three values)
EU_data_Chr13[which(EU_data_Chr13$CM_BP<=0),]  #no recombination until 2.6 Mbp?
ssa13_canada_bp_removed[which(ssa13_canada_bp_removed$CM_BP<=0),]  #no recombination until 5.4Mbp?


########################################################
################ Ssa16
########################################################

EU_data_Chr16=EU_data[which(EU_data$Chromosome==16),]
ssa16_canada_bp=Canada_data[which(Canada_data$Chromosome==16),]



Merged_eu_candaa_chr16_Recom=ggplot(data=ssa16_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=3.6e7, ymax=0, ymin=125, xmin=3.7e7, fill="gray")+
  geom_point(data=ssa16_canada_bp,aes(x=start, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr16,aes(x=start, y=Female.map, col="Europe"))+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+theme_bw()+
  ylab("CM")+xlab("Ssa16 Position")+  theme(panel.grid = element_blank(), legend.position = "top")


#Merged plot _ canada europe

merged_ssa16_CMbp=ggplot(data=ssa16_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=3.6e7, ymax=5, ymin=0, xmin=3.7e7, fill="gray")+
  geom_smooth(data=ssa16_canada_bp,
              aes(x=start, y=CM_BP, col="North Am"), stat="identity")+
  geom_smooth(data=EU_data_Chr16,
              aes(x=start, y=CM_BP, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa16 Position")+theme_bw()+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,3), expand = F)+
  NULL


#Check number of outliers for Ssa16 Europe and North AM
outliers(ssa16_canada_bp$CM_BP,2)
ssa16_canada_bp[which(ssa16_canada_bp$CM_BP<=0),]  #7 SNPS
ssa16_canada_bp[which(ssa16_canada_bp$CM_BP>=2.36670),]  #2 SNPS

#North AM - just 0 cM/mbp is outlier (two values)
outliers(EU_data_Chr16$CM_BP,2)
#Europe - just 0 cM/mbp is outlier (three values)
EU_data_Chr16[which(EU_data_Chr16$CM_BP>=2.423),]  #no recombinatio

########################################################
################ Ssa19
########################################################

EU_data_Chr19=EU_data[which(EU_data$Chromosome==19),]

##Remove 1 misplaced SNP
EU_data_Chr19_remove1=EU_data_Chr19[which(EU_data_Chr19$SNP_ID!="ESTV_16749_913"),]

ssa19_canada_bp=Canada_data[which(Canada_data$Chromosome==19),]
#remove 1 snps misplaced
ssa19_canada_bp_remove1=ssa19_canada_bp[which(ssa19_canada_bp$SNP_ID!="ESTNV_32145_342"),]


#Merge Position / CM
merged_ssa19_CM_Position_noremove=ggplot(data=ssa19_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=5.4e7, ymax=0, ymin=125, xmin=5.5e7, fill="gray")+
  geom_point(data=ssa19_canada_bp,aes(x=start, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr19,aes(x=start, y=Female.map,  col="Europe"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+theme_bw()+
  ylab("CM")+xlab("Ssa19 Position")+  theme(panel.grid = element_blank(), legend.position = "top")


merged_ssa19_CMbp_remove=ggplot(data=ssa19_canada_bp_remove1)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=5.4e7, ymax=0, ymin=5, xmin=5.5e7, fill="gray")+
  geom_smooth(data=ssa19_canada_bp_remove1,
              aes(x=start, y=CM_BP, col="North Am"), stat="identity")+
  geom_smooth(data=EU_data_Chr19_remove1,
              aes(x=start, y=CM_BP, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa19 Position")+theme_bw()+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,4), expand = F)+
  NULL

#Check number of outliers for Ssa19 Europe and North AM
outliers(ssa19_canada_bp_remove1$CM_BP,2) #only high recombination - no low
ssa19_canada_bp_remove1[which(ssa19_canada_bp_remove1$CM_BP>=3.687),]  #2 SNPS

outliers(EU_data_Chr19_remove1$CM_BP,2) #1 high value and 2 low values
#Europe - just 0 cM/mbp is outlier (three values)
EU_data_Chr19_remove1[which(EU_data_Chr19_remove1$CM_BP>2.208),]  #no recombinatio

#PLOTS
cowplot::plot_grid(merged_ssa06_CMbp, merged_ssa13_CMbp, merged_ssa16_CMbp, merged_ssa19_CMbp_remove,
                   Merged_eu_candaa_chr06_Recom,     Merged_eu_candaa_chr13_Recom,  Merged_eu_candaa_chr16_Recom,     merged_ssa19_CM_Position_noremove,  ncol=4,
                   labels=c("A", "B", "C", "D", "E", "F", "G", "H"))

#######################

#Read in gene annotations - downloaded .csv of all genes on NCBI for Atlantic Salmon
anno=read.csv("ALL_ANNOTAIONS.csv", header=T)

#add extra chromosome - without "ssa"
anno$chromosome2=as.factor(anno$chromosome)

#Remove any with NA for chromosome number
anno2=anno[!is.na(anno$chromosome2),]

#Add extra column wiht just '1's that will be used for counting
anno2$extra=rep(1)

#Get size of gene regions
anno2$size=anno2$end_position_on_the_genomic_accession-anno2$start_position_on_the_genomic_accession

#Use window scanr to get number of genes (sum of 1s) in each 1mbp window  - using start position 
gene_wind=winScan(x = anno2, groups = "chromosome2", position = "start_position_on_the_genomic_accession", win_size = 1000000, win_step = 1000000, 
                  values = "extra", funs = "sum")

#Plot results
genes_ssa06=ggplot(data=gene_wind[which(gene_wind$chromosome2=="ssa06"),], aes(win_mid, extra_n))+geom_rect(xmax=7.0e7, ymax=70, ymin=0, xmin=7.3e7, fill="gray")+theme_bw()+
  geom_area(stat="identity")+theme(panel.grid = element_blank(), legend.position = "top")+ylab("Gene count")+xlab("Position")+coord_cartesian(expand=F)

genes_ssa13=ggplot(data=gene_wind[which(gene_wind$chromosome2=="ssa13"),], aes(win_mid, extra_n))+geom_rect(xmax=7.5e7, ymax=70, ymin=0, xmin=7.8e7, fill="gray")+theme_bw()+
  geom_area(stat="identity")+theme(panel.grid = element_blank(), legend.position = "top")+ylab("Gene count")+xlab("Position")+coord_cartesian(expand=F)

genes_ssa16=ggplot(data=gene_wind[which(gene_wind$chromosome2=="ssa16"),], aes(win_mid, extra_n))+geom_rect(xmax=3.6e7, ymax=70, ymin=0, xmin=3.7e7, fill="gray")+theme_bw()+
  geom_area(stat="identity")+theme(panel.grid = element_blank(), legend.position = "top")+ylab("Gene count")+xlab("Position")+coord_cartesian(expand=F)

genes_ssa19=ggplot(data=gene_wind[which(gene_wind$chromosome2=="ssa19"),], aes(win_mid, extra_n))+geom_rect(xmax=5.4e7, ymax=70, ymin=0, xmin=5.5e7, fill="gray")+theme_bw()+
  geom_area(stat="identity")+theme(panel.grid = element_blank(), legend.position = "top")+ylab("Gene count")+xlab("Position")+coord_cartesian(expand=F)


#Run outlier script to add arrows/asterisk to figures (added in illustrator)
outliers(gene_wind$extra_n[which(gene_wind$chromosome2=="ssa06")], 2)
outliers(gene_wind$extra_n[which(gene_wind$chromosome2=="ssa13")], 2)
outliers(gene_wind$extra_n[which(gene_wind$chromosome2=="ssa16")], 2)
gene_wind[which(gene_wind$chromosome2=="ssa16" & gene_wind$extra_n>43),]
outliers(gene_wind$extra_n[which(gene_wind$chromosome2=="ssa19")], 2)
gene_wind[which(gene_wind$chromosome2=="ssa19" & gene_wind$extra_n>37),]

#Save PDF 14 x 12 - Figure 5
cowplot::plot_grid(  Merged_eu_candaa_chr06_Recom,     Merged_eu_candaa_chr13_Recom,  Merged_eu_candaa_chr16_Recom,     merged_ssa19_CM_Position_noremove, 
                     merged_ssa06_CMbp, merged_ssa13_CMbp, merged_ssa16_CMbp, merged_ssa19_CMbp_remove,
                     genes_ssa06, genes_ssa13, genes_ssa16, genes_ssa19,ncol=4,
                     labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))


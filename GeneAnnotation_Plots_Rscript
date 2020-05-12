setwd("~/Desktop/Sarah/Salmon/TransAtlanticPaper/")
library(ggplot2)

#This script plots the gene annotation information (downloaded from NCBI) with the q-values from pcadapt
#SNPs that are highly significant (low q-values) were determined in this script
#These SNPs were highlighted in the figure along with the gene that they are found within
#This script could be improved :) 


#All annotations from genome - downloaded all genes for Salmo salar from NCBI in .csv format
annot=data.table::fread("Recombination/ALL_ANNOTAIONS.csv")

#Pcadapt results (PC=1; K=1)
pcadapt=read.table("PCADAPT_185KLoci_Results_pvals_K1_maf001.txt", header=T)
pcadapt$Position=as.numeric(as.character(pcadapt$Position))
pcadapt=pcadapt[order(pcadapt$Qval, decreasing = T),]


#####Ssa13 annoations
#Annotation for Ssa13
ssa13_mbp=annot[which(annot$chromosome=="ssa13" & annot$start_position_on_the_genomic_accession>=7.4e07 & 
                        annot$end_position_on_the_genomic_accession< 7.8e07),]

#Check for top SNP
#topSNP in >7.7e7, position=77351757
ssa13_mbp[which(ssa13_mbp$start_position_on_the_genomic_accession>7.71e07 &
                  ssa13_mbp$start_position_on_the_genomic_accession< 7.74e07),]
#Get minimum q-value for Chr 13
min(pcadapt$qval1[which(pcadapt$Chr==13 & pcadapt$Position >= 7.4e7 & pcadapt$Position <= 7.8e7 )])

#Get SNP that has this minimum q-value
pcadapt[which(pcadapt$Chr==13 & pcadapt$Position >= 7.4e7 & pcadapt$Position <= 7.8e7 & pcadapt$qval1 <1.8e-17 ),]

#First plot q-values within the region of interest
plot_ssa13=ggplot(data=pcadapt[which(pcadapt$Chr==13 & pcadapt$Position >= 7.4e7 & pcadapt$Position <= 7.8e7),])+
  geom_point(aes(x=Position, y=-log10(qval1)), size=2.5, alpha=0.7, col="dodgerblue4")+
  xlim(7.4e07, 7.8e07)+ scale_y_continuous(limits = c(0,22), expand = c(0,0))+
  ylab("-log10(q)")+xlab("Ssa13 position")

#Next add  gene annotations (using geom_rect)
#And highlight genes of interest (based on looking at row in dataframe that belongs to the gene of interest) - this could be improved :) 
#In this case the gene is in row 5 and 77 in the dataframe ssa13_mbp
#Added gene label manual (see "label" in aes) - looked up gene name to find abbreviation
plot2_ssa13=plot_ssa13+ 
  geom_hline(yintercept = 19, col="gray70")+
  geom_rect(data=ssa13_mbp, aes(xmin=ssa13_mbp$start_position_on_the_genomic_accession,
                                xmax=ssa13_mbp$end_position_on_the_genomic_accession, 
                                ymin= rep(20, nrow(ssa13_mbp)), ymax=rep(18, nrow(ssa13_mbp))), fill="gray90", col="gray70")+
  geom_rect(data=ssa13_mbp, aes(xmin=ssa13_mbp$start_position_on_the_genomic_accession[5],
                                 xmax=ssa13_mbp$end_position_on_the_genomic_accession[5], 
                                 ymin= rep(20, 1), ymax=rep(18,1)), fill="pink", col='black')+
  #Added labels for genes that had highly significant SNPs (peaks in q-values) 
  #note these numbers pertain to these genes - adds name of gene (description)
  annotate(geom = "text", x = ssa13_mbp$start_position_on_the_genomic_accession[5], y=21,
           label="Cxadr" , size=4, angle=45)+
  geom_rect(data=ssa13_mbp, aes(xmin=ssa13_mbp$start_position_on_the_genomic_accession[77],
                                 xmax=ssa13_mbp$end_position_on_the_genomic_accession[77], 
                                 ymin= rep(20, 1), ymax=rep(18,1)), fill="goldenrod", col="black")+
  annotate(geom = "text", x = ssa13_mbp$end_position_on_the_genomic_accession[77], y=21,
           label="SLC25A14-l", size=4, angle=45)+
  NULL

#Highlight these highly significant snps
sign_snp=pcadapt[which(pcadapt$Chr==13 & pcadapt$Position>7.48e7 &  pcadapt$Position<7.55e7 & pcadapt$qval1< 1.5e-15),]
sign_snp2=pcadapt[which(pcadapt$Chr==13 & pcadapt$Position>7.7e7 &  pcadapt$Position<7.76e7 & pcadapt$qval1< 6.5e-13),]

#Add red diamonds for highly significant snps 
plot3_ssa13=plot2_ssa13+geom_point(data=sign_snp,  aes(x=Position, y=-log10(qval1)),col="red", pch=18,  size=4)+
  geom_point(data=sign_snp2,  aes(x=Position, y=-log10(qval1)),col="red", pch=18,  size=4)+theme_bw()+theme(panel.grid = element_blank())

#Check plot - merge all plots below
plot3_ssa13


#####Ssa06 annoations
#A similar script to aboe is used to plot gene annotations for Ssa06

#Pcadapt for ssa06
pcadapt_ssa06=pcadapt[which(pcadapt$Chr==6 & pcadapt$Position > 7e07 & pcadapt$Position < 7.3e07),]

ssa06_3MBP=annot[which(annot$chromosome=="ssa06" & annot$start_position_on_the_genomic_accession>7e+07 & annot$end_position_on_the_genomic_accession< 7.3e07),]

#Plot pcadapt q-values for Ssa06 region of interest
plot_ssa06=ggplot(pcadapt[which(pcadapt$Chr==6),])+
  geom_point(aes(x=Position, y=-log10(qval1)), size=2.5, alpha=0.7, col="dodgerblue4")+
  xlim(7e07, 7.3e07)+ scale_y_continuous(limits = c(0,80), expand = c(0,0))+
  ylab("-log10(q)")+xlab("Ssa06 position")


#Top SNP for Ssa06 (two peaks of interest)
snp1=pcadapt[which(pcadapt$Chr==6  & pcadapt$Position > 7e07 & pcadapt$Position <7.3e07 &
                pcadapt$qval1<1e-58 ),]
snp2=pcadapt[which(pcadapt$Chr==6  & pcadapt$Position > 7.2e07 & pcadapt$Position <7.3e07 &
                pcadapt$qval1<1e-35 ),]

#Add gene annotations to the plot and highlight important genes
##Highlight the two snps of interest (snp1 and snp2) - red diamond
#Also highlighted an additional gene (PPP2R5C-l) because was important in another study (see paper)
plot2_ssa06= plot_ssa06 + geom_hline(yintercept = 72.5, col="gray70")+
  geom_point(data=snp2, aes(x=Position, y=-log10(qval1)),col="red", pch=18,  size=4)+
  geom_point(data=snp1, aes(x=Position, y=-log10(qval1)),col="red", pch=18,  size=4)+
  geom_rect(data=ssa06_3MBP, aes(xmin=ssa06_3MBP$start_position_on_the_genomic_accession,
                                xmax=ssa06_3MBP$end_position_on_the_genomic_accession, 
                                ymin= rep(70, 62), ymax=rep(75, 62)), fill="gray90", col="gray70")+
  geom_rect(data=ssa06_3MBP, aes(xmin=ssa06_3MBP$start_position_on_the_genomic_accession[11],
                                 xmax=ssa06_3MBP$end_position_on_the_genomic_accession[11], 
                                 ymin= rep(70, 1), ymax=rep(75,1)), fill="cornflowerblue", col="black")+
  annotate(geom = "text", x = ssa06_3MBP$end_position_on_the_genomic_accession[11], y=78,
           label="EML1-l" , size=4, angle=45)+
  geom_rect(data=ssa06_3MBP, aes(xmin=ssa06_3MBP$start_position_on_the_genomic_accession[51],
                                 xmax=ssa06_3MBP$end_position_on_the_genomic_accession[51], 
                                 ymin= rep(70, 1), ymax=rep(75,1)), fill="palegreen", col="black")+
  annotate(geom = "text", x = ssa06_3MBP$end_position_on_the_genomic_accession[51], y=78,
           label="snap25a" , size=4, angle=45)+theme_bw()+theme(panel.grid = element_blank())+
  geom_rect(data=ssa06_3MBP, aes(xmin=ssa06_3MBP$start_position_on_the_genomic_accession[43],
                               xmax=ssa06_3MBP$end_position_on_the_genomic_accession[43], 
                               ymin= rep(70, 1), ymax=rep(75,1)), fill="yellow", col="black")+
  annotate(geom = "text", x = ssa06_3MBP$end_position_on_the_genomic_accession[43], y=78,
           label="PPP2R5C-l" , size=4, angle=45)+theme_bw()+theme(panel.grid = element_blank())

plot2_ssa06

####
######Ssa16 annoations

ssa16_1MBP=annot[which(annot$chromosome=="ssa16" & annot$start_position_on_the_genomic_accession>=3.597e7 & 
                         annot$end_position_on_the_genomic_accession< 3.7e07),]
#Max SNPs
pcadapt[which(pcadapt$Chr==16 & pcadapt$qval1< 1e-30),]

#Plot q-values for Ssa16 in region of interest
plot_ssa16=ggplot(pcadapt[which(pcadapt$Chr==16),])+
  geom_point(aes(x=Position, y=-log10(qval1)), size=2.5, alpha=0.7, col="dodgerblue4")+
  xlim(3.59e07, 3.7e07)+ scale_y_continuous(limits = c(0,65), expand = c(0,0))+
  ylab("-log10(q)")+xlab("Ssa16 position")

##Add gene of interest (one is Uncharacterized - so also included the next closest gene)
#Highlight top SNP (based on lowest qvalue)
plot2_ssa16= plot_ssa16+ geom_hline(yintercept = 57.5, col="gray70")+
  geom_point(data=pcadapt[which(pcadapt$Chr==16 & pcadapt$qval1< 2.3e-43),],
             aes(x=pcadapt$Position[which(pcadapt$Chr==16 & pcadapt$qval1< 2.3e-43
                                            )],  y=-log10(pcadapt$qval1[which(pcadapt$Chr==16 & pcadapt$qval1< 2.3e-43)])),
             col="red", pch=18,  size=4)+
  geom_rect(data=ssa16_1MBP, aes(xmin=ssa16_1MBP$start_position_on_the_genomic_accession,
                                 xmax=ssa16_1MBP$end_position_on_the_genomic_accession, 
                                 ymin= rep(56, nrow(ssa16_1MBP)), ymax=rep(59, nrow(ssa16_1MBP))), fill="gray90", col="gray70")+
  geom_rect(data=ssa16_1MBP, aes(xmin=ssa16_1MBP$start_position_on_the_genomic_accession[5],
                                 xmax=ssa16_1MBP$end_position_on_the_genomic_accession[5], 
                                 ymin= rep(56, 1), ymax=rep(59,1)), fill="indianred", col="black")+
  annotate(geom = "text", x = ssa16_1MBP$start_position_on_the_genomic_accession[5], y=61,
           label="SV2B-l" , size=4, angle=45)+
  geom_rect(data=ssa16_1MBP, aes(xmin=ssa16_1MBP$start_position_on_the_genomic_accession[6],
                                 xmax=ssa16_1MBP$end_position_on_the_genomic_accession[6], 
                                 ymin= rep(56, 1), ymax=rep(59,1)), fill="purple", col="black")+
  annotate(geom = "text", x = ssa16_1MBP$end_position_on_the_genomic_accession[6], y=61,
           label="Unchr" , size=4, angle=45)+theme(panel.grid = element_blank())+theme_bw()+theme(panel.grid = element_blank())

plot2_ssa16


######Ssa19 annoations

ssa19_1MBP=annot[which(annot$chromosome=="ssa19" & annot$end_position_on_the_genomic_accession >5.3e+07 & annot$end_position_on_the_genomic_accession< 5.5e07),]

#Plot q-values for Ssa19 region of interest
plot_ssa19=ggplot(pcadapt[which(pcadapt$Chr==19),])+
  geom_point(aes(x=Position, y=-log10(qval1)), size=2.5, alpha=0.7, col="dodgerblue4")+
  xlim(5.3e07, 5.5e07)+ scale_y_continuous(limits = c(0,32), expand = c(0,0))+
  ylab("-log10(q)")+xlab("Ssa19 position")


#Add gene annotations
#Highlight gene of interest 
#And highlight most significant SNP in region
plot2_ssa19 = plot_ssa19+ geom_hline(yintercept = 27, col="gray70")+
  geom_point(data=pcadapt[which(pcadapt$Chr==19 & pcadapt$qval1< 2.3e-23),],
             aes(x=pcadapt$Position[which(pcadapt$Chr==19 & pcadapt$qval1< 2.3e-23
             )],  y=-log10(pcadapt$qval1[which(pcadapt$Chr==19 & pcadapt$qval1< 2.3e-23)])),
             col="red", pch=18,  size=4)+
  geom_rect(data=ssa19_1MBP, aes(xmin=ssa19_1MBP$start_position_on_the_genomic_accession,
                                 xmax=ssa19_1MBP$end_position_on_the_genomic_accession, 
                                 ymin= rep(26, nrow(ssa19_1MBP)), ymax=rep(28, nrow(ssa19_1MBP))), fill="gray90", col="gray70")+ 
                              theme(panel.grid = element_blank())+
  geom_rect(data=ssa19_1MBP, aes(xmin=ssa19_1MBP$start_position_on_the_genomic_accession[44],
                                 xmax=ssa19_1MBP$end_position_on_the_genomic_accession[44], 
                                 ymin= 26, ymax=28), fill="purple", col="black")+
 annotate(geom = "text", x = ssa19_1MBP$end_position_on_the_genomic_accession[44], y=30,
         label=ssa19_1MBP$Symbol[44],  size=4, angle=45)+theme_bw()+
  theme(panel.grid = element_blank())



###Merge all plots together with cowplot
plot1=cowplot::plot_grid(plot2_ssa06, plot3_ssa13, nrow=2)
plot2=cowplot::plot_grid(plot2_ssa16, plot2_ssa19, nrow=1)

#Save as 16 x 10 
cowplot::plot_grid(plot1, plot2,rel_heights = c(2,1)  ,nrow=2)
#Edited this a bit in illustrator after exporting

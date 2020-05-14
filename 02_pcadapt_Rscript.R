#R packages
library(pcadapt)
library(qvalue)
library(data.table)
library(qqman)
library(ggplot2)

#Set working directory
setwd("~/Desktop/Sarah/Salmon/TransAtlanticPaper/")

#Run PCAdapt with new file
nor_can_salmon185K <- read.pcadapt("Nor_Can_GOOD_SNPs_maf001.bed",type="bed")

#Use k=2 to get split between Europe an North America on both axes for plotting PCA
pcadapt_nor_can_salmon_2pc <- pcadapt(nor_can_salmon185K, method = "mahalanobis",min.maf = 0.01, K=2)

#Plot screeplot
plot(pcadapt_nor_can_salmon_2pc,option="screeplot")

#Plot individual scores
plot(pcadapt_nor_can_salmon_2pc,"scores")

#Check variance explained by each axis
pcadapt_nor_can_salmon_2pc$singular.values^2

#Read in individual data
file_info <- fread("Nor_Can_GOOD_SNPs_maf001.fam", header=F)
pop_list2 <- as.data.frame(cbind(file_info$V1, file_info$V2)) #individual and Pop info

#List of EU and North American pops
pop_info_plot <- read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/pop_info_CAN_NOR.txt", header=F)

#Quick plot of pcadapt
plot(pcadapt_nor_can_salmon_2pc,option="scores", Pop=pop_info_plot$V3)

#Data for plotting
plot_data <- as.data.frame(cbind(pcadapt_nor_can_salmon_2pc$scores,pop_info_plot$V3))
colnames(plot_data) <- c("PC1", "PC2", "Pop")
plot_data$Pop <- as.factor(plot_data$Pop) #Pop here indicates Europe vs North America (2 is Europe)

#save as 8x8 pdf
ggplot(data=plot_data, aes(x=PC1, y=PC2))+
  geom_point(size=4, aes(fill=Pop), shape=21, color="black", alpha=0.7)+
  scale_fill_manual(values=c("blue", "red"))+theme_bw()+
  theme(panel.grid = element_blank(), axis.title = element_text(size=20), axis.text = element_text(size=15))
#Canada in blue and Europe in red

##################################################################
#Rerun with K=1 to get outlier loci associated with axis 1 (axis that seperates Europe and North America)
pcadapt_nor_can_salmon_1pc <- pcadapt(nor_can_salmon185K,
                                      method = "mahalanobis", K=1, min.maf = 0.01)

#Get outliers based on q-values
alpha <- 0.05
qval1 <- qvalue(pcadapt_nor_can_salmon_1pc$pvalues)$qvalues
outliers1 <- which(qval1<alpha)
length(outliers1) #number of outliers 

#Read table of SNPs in same order as genotype file (used .map file as .txt)
SNPs <- read.table("Nor_Can_GOOD_SNPs_maf001.map")

#Create dataframe with qvalues
qval1 <- as.data.frame(qval1) 

#Combine SNP names and qvalue data
results_pc1 <- cbind(qval1, SNPs) #qvalues are in order of SNPs in map file 
colnames(results_pc1)<-c("qval1", "Chr", "SNP", "x", "Position")

#Save as final dataframe
Results2_PC1 <- results_pc1

#Make sure Chr and Position are numeric values
Results2_PC1$Chr=as.numeric(as.character(Results2_PC1$Chr))
Results2_PC1$Position=as.numeric(as.character(Results2_PC1$Position))

#Save data for other analyses
write.table(Results2_PC1, "PCADAPT_185KLoci_Results_pvals_K1_maf001.txt", col.names = T, row.names = F, quote = F, sep="\t")
write.table(Results2_PC1[which(Results2_PC1$qval1<0.05),], "PCADAPT_outliers_qvals_185K_K1_maf001.txt", col.names = T, row.names = F, quote = F, sep="\t")

#number of outliers
nrow(Results2_PC1[which(Results2_PC1$qval1<0.05),])
#Percentage of outliers
nrow(Results2_PC1[which(Results2_PC1$qval1<0.05),])/nrow(Results2_PC1)

#Remove blanks/NAs as won't work in Manhattan plot
#save as 18 x 6 pdf
manhattan(x=Results2_PC1[!is.na(Results2_PC1$Chr),],
          chr="Chr", bp ="Position", p = "qval1", snp = "SNP",
          genomewideline = F,
          suggestiveline =  F,  #col = c("dodgerblue4", "gray48"),
          ylim=c(0,90),col=c("gray50", "dodgerblue4"),
          highlight=Results2_PC1$SNP[order(Results2_PC1$qval1, decreasing=F)][1:100])
#Added gene names manually in Illustrator


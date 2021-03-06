#Packages
library(genepopedit)
library(reshape2)
library(ggplot2)
library(gplots)
library(cowplot)

#This R script include 
  # 1 - Plotting mean LD across chromosome (outlier regions determined by top 5% cutoff) and smoothed line for FST

#Script for Ssa06 included here - same for other chromosomes

#Set directory
setwd("/Users/Bradbury/Desktop/Sarah/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/AlleleFreq/Whole_Chromosome_LD/")

###########

#Read in fst values - as they were included in the plot too - sliding window values
fst=data.table::fread("/Users/Bradbury/Desktop/Sarah/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/Fst/FST_SlidingWindow_RESULTS_185K_norCan_maf001.txt", header=T)

##Sliding window LD

####Read in LD values (r2)
data_Chr6_nor=read.table("LD_50nearbySNPWindow_Chr6_NORWAY_185K.txt", header=T)
data_Chr6_canada=read.table("LD_50nearbySNPWindow_Chr6_CANADA_maf001.txt", header=T)

#considered SNPs in the top 5% as those that are outliers
max(data_Chr6_nor$Mean_LD, na.rm=T)
nrow(data_Chr6_nor[!is.na(data_Chr6_nor$Mean_LD),])*0.05

Chr6_Norway=ggplot2::ggplot()+
  geom_rect(aes(ymin=0, ymax=0.85, xmin=7.0e7, xmax=7.3e7), col="gray80", alpha=0.4)+
  geom_point(data=data_Chr6_nor, aes(data_Chr6_nor$Position, data_Chr6_nor$Mean_LD, col="Norway"), alpha=0.5, size=2,inherit.aes = F)+
  scale_color_manual(values=c( "dodgerblue3"))+
  #coord_cartesian(ylim=c(0,0.1))+
  geom_line(data=fst[which(fst$CHR==6),], aes(x = win_mid, y=FST_mean),lwd=1.2,
            col="navy", stat = "identity")+
  geom_hline(yintercept = min(data_Chr6_nor$Mean_LD[order(data_Chr6_nor$Mean_LD,decreasing = T)][1: 382]), lty=2)+
  
  ylab("Mean LD (Sliding window)")+xlab("Position on Chromosome 6")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.border = element_blank(), text = element_text(size=20))+geom_smooth()+
  NULL

#considered SNPs in the top 5% as those that are outliers
max(data_Chr6_canada$Mean_LD, na.rm=T)
max(fst$FST_mean, na.rm=T)
nrow(data_Chr6_canada[!is.na(data_Chr6_canada$Mean_LD),])*0.05

Chr6_Canada=ggplot2::ggplot()+
  geom_rect(aes(ymin=0, ymax=0.75, xmin=7.0e7, xmax=7.3e7), col="gray80", alpha=0.4)+
  geom_point(data=data_Chr6_canada, aes(data_Chr6_canada$Position,data_Chr6_canada$Mean_LD, col="Canada"), alpha=0.5, size=2)+
  scale_color_manual(values=c("firebrick2"))+
  #coord_cartesian(ylim=c(0,0.1))+
  geom_line(data=fst[which(fst$CHR==6),], aes(x = win_mid, y=FST_mean),lwd=1.2,
            col="navy", stat = "identity")+
  geom_hline(yintercept = min(data_Chr6_canada$Mean_LD[order(data_Chr6_canada$Mean_LD,decreasing = T)][1: 367]), lty=2)+
  
  ylab("Mean LD (Sliding window)")+xlab("Position on Chromosome 6")+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(),
        axis.line.y = element_line(),
        panel.border = element_blank(), text = element_text(size=20))+geom_smooth()+
  NULL


plot_grid(Chr6_Canada, Chr6_Norway)



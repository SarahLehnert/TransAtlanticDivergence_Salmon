#Libraries
library(ggplot2)
library(cowplot)

#This R script includes smoothed line plot for chromosomes of interest
#This is for Figure 3 in the manuscript
#Plots include result from Fst, pcadapt, twisst, and Tajima's D
#Data were read in and plots were created and combined using cowplot
#Data sets were examined for max/min values to set y coordinates
#Plots were further editted in Illustrator

###############

##Plots for FST
sliding_FST=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/Fst/FST_SlidingWindow_RESULTS_185K_norCan_maf001.txt", header=T)
sliding_FST$Chr=sliding_FST$CHR

#only keep windows with >= 5 SNPs used for calculations
sliding_FST=sliding_FST[which(sliding_FST$FST_n >=5),]

plot_chr6_fst=ggplot() + 
  geom_rect(aes(ymin=0, 
                ymax=max(sliding_FST$FST_mean[which(sliding_FST$Chr==6)])+0.05,
                xmin=7.0e7, xmax=7.3e7  ), fill="gray80")+
  coord_cartesian(ylim=c(0, max(sliding_FST$FST_mean[which(sliding_FST$Chr==6)])+0.05),
                  expand=F)+
  geom_smooth(data = sliding_FST[which(sliding_FST$Chr==6),], 
              aes(x=win_mid, y=FST_mean, colour = "All"), #stat="identity")+
              method="loess",  span = 0.1, se=FALSE) + 
  scale_colour_manual(values="blue")+ylab("Mean FST")+xlab("")+
  theme_classic() + #geom_vline(xintercept = 3.6e+7) + geom_vline(xintercept = 3.7e+7)
  NULL

plot_chr13_fst=ggplot() +
  geom_rect(aes(ymin=0, 
                ymax=max(sliding_FST$FST_mean[which(sliding_FST$Chr==13)])+0.05,
                xmin=7.5e7, xmax=7.8e7  ), fill="gray80")+
  coord_cartesian(ylim=c(0, max(sliding_FST$FST_mean[which(sliding_FST$Chr==13)])+0.05),
                  expand=F)+
  geom_smooth(data = sliding_FST[which(sliding_FST$Chr==13),], 
              aes(x=win_mid, y=FST_mean, colour = "All"),
              #stat="identity")+
              method="loess",  span = 0.1, se=FALSE) + 
  scale_colour_manual(values="darkgreen")+ylab("")+xlab("")+
  theme_classic() + #geom_vline(xintercept = 3.6e+7) + geom_vline(xintercept = 3.7e+7)
  NULL

plot_chr16_fst=ggplot() +
  geom_rect(aes(ymin=0, 
                ymax=max(sliding_FST$FST_mean[which(sliding_FST$Chr==16)])+0.05,
                xmin=3.6e7, xmax=3.7e7 ), fill="gray80")+
  coord_cartesian(ylim=c(0,  max(sliding_FST$FST_mean[which(sliding_FST$Chr==16)])+0.05),
                  expand=F)+
  geom_smooth(data = sliding_FST[which(sliding_FST$Chr==16),], 
              aes(x=win_mid, y=FST_mean, colour = "All"), #stat="identity")+
              method="loess",  span = 0.1, se=FALSE) + ylab("")+xlab("")+
  scale_colour_manual(values="orange")+
  
  theme_classic() + 
  NULL


plot_chr19_fst=ggplot() +
  geom_rect(aes(ymin=0, 
                ymax=max(sliding_FST$FST_mean[which(sliding_FST$Chr==19)])+0.05,
                xmin=5.4e7, xmax=5.5e7 ), fill="gray80")+
  coord_cartesian(ylim=c(0, max(sliding_FST$FST_mean[which(sliding_FST$Chr==19)])+0.05),
                  expand=F)+
  
  geom_smooth(data = sliding_FST[which(sliding_FST$Chr==19),], 
              aes(x=win_mid, y=FST_mean, colour = "All"), #stat="identity")+
              method="loess",  span = 0.1, se=FALSE) + ylab("")+xlab("")+
  scale_colour_manual(values="red")+
  
  theme_classic() + 
  NULL

##Plots for PCAdapt
sliding_pcadapt=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/pcadapt/PCA_SlidingWindow_k1_Maf001_canNor_185K.txt", header=T)
sliding_pcadapt$logp_mean=-log10(sliding_pcadapt$qval1_mean)

#only keep windows with >= 5 SNPs used for calculations
sliding_pcadapt=sliding_pcadapt[which(sliding_pcadapt$qval1_n >= 5),]

plot_chr6_pcadapt=ggplot() + 
  geom_rect(aes(ymin=0, 
                ymax=max(sliding_pcadapt$logp_mean[which(sliding_pcadapt$Chr==6)])+0.1,
                xmin=7.0e7, xmax=7.3e7), fill="gray80")+ylab("Mean -log10(q)")+xlab("")+
  coord_cartesian(ylim=c(0, max(sliding_pcadapt$logp_mean[which(sliding_pcadapt$Chr==6)])+0.1),
                  expand=F)+
  geom_smooth(data = sliding_pcadapt[which(sliding_pcadapt$Chr==6),], 
              aes(x=win_mid, y=logp_mean, colour = "All"),# stat="identity")+
              method="loess",  span = 0.1, se=FALSE) + 
  scale_colour_manual(values="blue")+
  theme_classic() + #geom_vline(xintercept = 3.6e+7) + geom_vline(xintercept = 3.7e+7)
  NULL

plot_chr13_pcadapt=ggplot() +
  geom_rect(aes(ymin=0, 
                ymax=max(sliding_pcadapt$logp_mean[which(sliding_pcadapt$Chr==13)])+0.1,
                xmin=7.5e7, xmax=7.8e7), fill="gray80")+
  coord_cartesian(ylim=c(0, max(sliding_pcadapt$logp_mean[which(sliding_pcadapt$Chr==13)])+0.1),
                  expand=F)+ylab("")+xlab("")+
  geom_smooth(data = sliding_pcadapt[which(sliding_pcadapt$Chr==13),], 
              aes(x=win_mid, y=logp_mean, colour = "All"),
              #  stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  scale_colour_manual(values="darkgreen")+
  theme_classic() + #geom_vline(xintercept = 3.6e+7) + geom_vline(xintercept = 3.7e+7)
  NULL


plot_chr16_pcadapt=ggplot() +
  geom_rect(aes(ymin=0, 
                ymax=0.3,
                xmin=3.6e7, xmax=3.7e7), fill="gray80")+
  coord_cartesian(ylim=c(0, 0.3),
                  expand=F)+ylab("")+xlab("")+
  # geom_point(data = sliding_pcadapt[which(sliding_pcadapt$Chr==16),], 
  #            aes(x=win_mid, y=logp_mean, colour = "All",alpha=0.3))+
  
  geom_smooth(data = sliding_pcadapt[which(sliding_pcadapt$Chr==16),], 
              aes(x=win_mid, y=logp_mean, colour = "All"), #stat="identity")+
              method="loess",  span = 0.1, se=FALSE) + 
  scale_colour_manual(values="orange")+
  
  theme_classic() +
  NULL


plot_chr19_pcadapt=ggplot() +
  geom_rect(aes(ymin=0, 
                ymax=0.4,
                xmin=5.4e7, xmax=5.5e7), fill="gray80")+
  coord_cartesian(ylim=c(0,0.4),
                  expand=F)+ ylab("")+xlab("")+
  geom_smooth(data = sliding_pcadapt[which(sliding_pcadapt$Chr==19),], 
              aes(x=win_mid, y=logp_mean, colour = "All"), #stat="identity")+
              method="loess",  span = 0.1, se=FALSE) + 
  # geom_point(data = sliding_pcadapt[which(sliding_pcadapt$Chr==19),], 
  #             aes(x=win_mid, y=logp_mean, colour = "All",alpha=0.3))+
  scale_colour_manual(values="red")+
  
  theme_classic() +
  NULL




##Plots for Tajimas
#Outliers for Tajima's D were determiend based on window values lower than the 5th percentile. (indicative of postive selection)
#These were indicated by asterisks in the plot (added in ggplot2)
#Tajima's D was calculated with PopGenome R package

sliding_Tajima6=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/Tajima/Updated_184K/TajimaD_each_Chr6.txt", header=T)

plot_chr6_tajima=ggplot() +
  geom_rect(aes(ymin=min(sliding_Tajima6$TajimaD_NA)-0.5, 
                ymax=max(sliding_Tajima6$TajimaD_ALL)+0.5,
                xmin=7.0e7, xmax=7.3e7), fill="gray80")+
  coord_cartesian(ylim=c(min(sliding_Tajima6$TajimaD_NA)-0.5, max(sliding_Tajima6$TajimaD_ALL)+0.5),
                  expand=F)+
  geom_smooth(data = sliding_Tajima6,
              aes(x=Position , y=TajimaD_EU, colour = "EU"),# stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  geom_smooth(data = sliding_Tajima6,
              aes(x=Position , y=TajimaD_NA, colour = "NA"),# stat="identity")+
              method="loess",  span = .1, se=FALSE) +  ylab("Tajima's D")+xlab("")+
  theme_classic()  +
  geom_text(aes(x=7.2e7, y=7), label="*", size=10, col="hotpink")+
  scale_colour_manual(values=c("hotpink", "skyblue2"))
NULL

sliding_Tajima13=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/Tajima/Updated_184K/TajimaD_each_Chr13.txt", header=T)

plot_chr13_tajima=ggplot() + 
  geom_rect(aes(ymin=min(sliding_Tajima13$TajimaD_NA)-0.5, 
                ymax=max(sliding_Tajima13$TajimaD_ALL)+0.5,
                xmin=7.5e7, xmax=7.8e7), fill="gray80")+
  coord_cartesian(ylim=c(min(sliding_Tajima13$TajimaD_NA)-0.5, max(sliding_Tajima13$TajimaD_ALL)+0.5),
                  expand=F)+
  geom_smooth(data = sliding_Tajima13,
              aes(x=Position , y=TajimaD_EU, colour = "EU"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  geom_smooth(data = sliding_Tajima13,
              aes(x=Position , y=TajimaD_NA, colour = "NA"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  theme_classic() +
  ylab("")+xlab("")+
  geom_text(aes(x=7.7e7, y=7), label="*", size=10, col="hotpink")+
  scale_colour_manual(values=c("hotpink", "skyblue2"))

NULL

sliding_Tajima16=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/Tajima/Updated_184K/TajimaD_each_Chr16.txt", header=T)

plot_chr16_tajima=ggplot() +
  geom_rect(aes(ymin=min(sliding_Tajima16$TajimaD_NA)-0.5, 
                ymax=max(sliding_Tajima16$TajimaD_ALL)+0.5,
                xmin=3.6e7, xmax=3.7e7), fill="gray80")+
  coord_cartesian(ylim=c(min(sliding_Tajima16$TajimaD_NA)-0.5, max(sliding_Tajima16$TajimaD_ALL)+0.5),
                  expand=F)+
  geom_smooth(data = sliding_Tajima16,
              aes(x=Position , y=TajimaD_EU, colour = "EU"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) +  ylab("")+xlab("")+
  geom_smooth(data = sliding_Tajima16,
              aes(x=Position , y=TajimaD_NA, colour = "NA"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  geom_text(aes(x=3.65e7, y=6.5), label="*", size=10, col="hotpink")+
  geom_text(aes(x=3.6e7, y=6.5), label="*", size=10, col="skyblue2")+
  scale_colour_manual(values=c("hotpink", "skyblue2"))+
  theme_classic()  
NULL


sliding_Tajima19=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/Tajima/Updated_184K/TajimaD_each_Chr19.txt", header=T)

plot_chr19_tajima=ggplot() +
  geom_rect(aes(ymin=min(sliding_Tajima19$TajimaD_NA)-0.5, 
                ymax=max(sliding_Tajima19$TajimaD_ALL)+0.5,
                xmin=5.4e7, xmax=5.5e7), fill="gray80")+ ylab("")+xlab("")+
  coord_cartesian(ylim=c(min(sliding_Tajima19$TajimaD_NA)-0.5, max(sliding_Tajima19$TajimaD_ALL)+0.5),
                  expand=F)+
  geom_smooth(data = sliding_Tajima19,
              aes(x=Position , y=TajimaD_EU, colour = "EU"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  geom_smooth(data = sliding_Tajima19,
              aes(x=Position , y=TajimaD_NA, colour = "NA"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  scale_colour_manual(values=c("hotpink", "skyblue2"))+
  theme_classic()  
NULL


##Plotting results from twisst

###Tree concordance

#Results for Ssa06
tree_6=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/twisst/chr6_NAEU_twist.txt", header=T)
#The weighting of topology 1 was calculated as a proportion
tree_6$Prop_Topo1=tree_6$topo1/((tree_6$topo1+tree_6$topo2+tree_6$topo3))

#Results for Ssa13, The weighting of topology 1 was calculated as a proportion
tree_13=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/twisst/chr13_NAEU_twist.txt", header=T)
tree_13$Prop_Topo1=tree_13$topo1/((tree_13$topo1+tree_13$topo2+tree_13$topo3))

#Results for Ssa16, The weighting of topology 1 was calculated as a proportion
tree_16=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/twisst/chr16_NAEU_twist.txt", header=T)
tree_16$Prop_Topo1=tree_16$topo1/((tree_16$topo1+tree_16$topo2+tree_16$topo3))

#Results for Ssa19, The weighting of topology 1 was calculated as a proportion
tree_19=read.table("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/twisst/chr19_NAEU_twist.txt", header=T)
tree_19$Prop_Topo1=tree_19$topo1/((tree_19$topo1+tree_19$topo2+tree_19$topo3))

plot_chr6_tree=ggplot() +
  geom_rect(aes(ymin=0.3, 
                ymax=1.05,
                xmin=7.0e7, xmax=7.3e7), fill="gray80")+
  coord_cartesian(ylim=c(0.3, 1.05),
                  expand=F)+ ylab("Tree concordance")+xlab("")+
  geom_hline(yintercept = 0.8, lty=2)+
  geom_smooth(data = tree_6,
              aes(x=start , y=Prop_Topo1,  col="All"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  scale_colour_manual(values="blue")+
  
  theme_classic()  
NULL


plot_chr13_tree=ggplot() +
  geom_rect(aes(ymin=0.3 ,
                ymax=1.05,
                xmin=7.5e7, xmax=7.8e7), fill="gray80")+
  geom_hline(yintercept = 0.8, lty=2)+
  coord_cartesian(ylim=c(0.3, 1.05),
                  expand=F)+ ylab("")+xlab("")+
  geom_smooth(data = tree_13,
              aes(x=start , y=Prop_Topo1, col="All"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  scale_colour_manual(values="darkgreen")+
  
  theme_classic()  
NULL



plot_chr16_tree=ggplot() +
  geom_rect(aes(ymin=0.3, 
                ymax=1.05,
                xmin=3.6e7, xmax=3.7e7), fill="gray80")+
  coord_cartesian(ylim=c(0.3, 1.05),
                  expand=F)+ ylab("")+xlab("")+
  geom_hline(yintercept = 0.8, lty=2)+
  geom_smooth(data = tree_16,
              aes(x=start , y=Prop_Topo1, col="All"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  scale_colour_manual(values="orange")+
  
  theme_classic()  
NULL


plot_chr19_tree=ggplot() +
  geom_rect(aes(ymin=0.3, 
                ymax=1.05,
                xmin=5.4e7, xmax=5.5e7), fill="gray80")+
  coord_cartesian(ylim=c(0.3, 1.05),
                  expand=F)+ ylab("")+xlab("")+
  geom_hline(yintercept = 0.8, lty=2)+
  geom_smooth(data = tree_19,
              aes(x=start , y=Prop_Topo1, col="All"), #stat="identity")+
              method="loess",  span = .1, se=FALSE) + 
  scale_colour_manual(values="red")+
  theme_classic()  
NULL

#NOTE - checked significance for Tajima's D for new results -184 K SNPs  - updated in figure - 
#Outliers for Tajima's D were determiend based on window values lower than the 5th percentile. (indicative of postive selection)
#added asterisks to plots

#Save as PDF -- size? 14 x 14
plot_grid( plot_chr6_fst, plot_chr13_fst,  plot_chr16_fst,plot_chr19_fst,
           plot_chr6_pcadapt, plot_chr13_pcadapt,  plot_chr16_pcadapt,plot_chr19_pcadapt,
           plot_chr6_tree, plot_chr13_tree, plot_chr16_tree, plot_chr19_tree,
           plot_chr6_tajima, plot_chr13_tajima, plot_chr16_tajima, plot_chr19_tajima,
           nrow=4)


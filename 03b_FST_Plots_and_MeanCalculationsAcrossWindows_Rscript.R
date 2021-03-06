#Set directory
setwd("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/Fst/")

#Libraries
library(data.table)
library(qqman)
library(windowscanr)
library(cowplot)
library(ggplot2)

#This R script include 
  # 1 - Manhattan plot for all FST values
  # 2 - Calculations of mean FST per 1 MBP window 
  # 3 - Manhattan plot of mean FST values (with outliers)
  # 4 - Smoothed line plots of mean FST values for chromosomes of interest

#Note that the same script was used to for mean pcadapt results, where mean q-value per 1mbp was
#calculated in the same way using windowscanr, and plotted with qqman for manhattan and ggplot2 for smoothed line plots 

#In both FST and q-values - outliers (based on 1mbp windows) were determined as values that were
#greater than ±3 SD from the mean, using the outlier function below

####################### Manhattan plot of all FST values ####################### 
#Read in FST values
fst_results=fread("Nor_Can_GOOD_SNPs_maf001_fst.fst", header=T)

#Check top 100 snps
all_data=fst_results[order(fst_results$FST, decreasing = T),]
all_data[1:100,]

#A few SNPs have no chromosome postion
all_data[is.na(all_data$CHR)]

#Check max FST value
max(all_data$FST)

#Save manhattan as 18 x 6 inch pdf
#plot for Figure 1
manhattan(x = all_data[!is.na(all_data$CHR),], chr = "CHR", bp = "POS", snp = "SNP",
          p = "FST", logp = F, cex=1, col=c("gray50", "dodgerblue4"), ylab="FST", 
          highlight=as.character(all_data$SNP[order(all_data$FST, decreasing=T)])[1:100],
         ylim=c(0.7,1))
#Note only plotted top values (FST> 0.7) to show

####################### Get mean FST per 1Mbp window  ####################### 

## For WindowScanr - to calculate mean FST per 1Mbp window

#add extra fake position at 1MBP for the max. --- see below

fake_data = list()

#for loop gets max value for each chromosome
#If this is not done, then it does not calculate values for last window on chromosome
#Note that used FST = NA for these 'fake' position - so does not influence results

for(i in 1:30){
  dat_EU=as.data.frame((max(all_data$POS[which(all_data$CHR==i)])))
  dat_EU$i <- i  # maybe you want to keep track of which iteration produced it
  fake_data[[i]] <- dat_EU # add it to your list  }
}
fake_data_chr = do.call(rbind, fake_data)
colnames(fake_data_chr)=c("Max_Pos", "Chr")
fake_data_chr$Fake_Pos=fake_data_chr$Max_Pos+1e6

#Create fake data frame to add to include extra windows at the end of chromosome ... Fake position have extra 1MBP added to the chromosome
#Make dataframe the same as the one you want to add to

head(all_data)
fake_data_to_add_FST=as.data.frame(cbind("CHR"=fake_data_chr$Chr, 
  "SNP"=paste0("AX-fake",1:30), 
  "POS"= fake_data_chr$Fake_Pos,
   "NMISS"=rep(1400), "FST"=rep(NA) ))

#remove Chr that have NA from original dataframe 
all_data=all_data[!is.na(all_data$CHR),]

#rbind both datasets (fake and real)
all_data_plus_extra=rbind(as.data.frame(all_data), fake_data_to_add_FST)

#Change anything to numeric if needed
all_data_plus_extra$POS=as.numeric(as.character(all_data_plus_extra$POS))
all_data_plus_extra$FST=as.numeric(as.character(all_data_plus_extra$FST))

#Run windowscanr - get mean FST for 1Mbp windows across each chromosome
fst_sliding=winScan(x = all_data_plus_extra, groups = "CHR", position = "POS",
        win_size = 1000000, win_step = 1000000,  values = "FST", funs = "mean")

#Add "chr/position' column (called snp) for reference
fst_sliding$snp=interaction(fst_sliding$CHR, fst_sliding$win_start, sep="_")
fst_sliding$FST_mean=as.numeric(as.character(fst_sliding$FST_mean))

#Retain windows with at least 1 SNP used for the calculations
fst_sliding2 <- fst_sliding[which(fst_sliding$FST_n>0),]

#Save results
write.table(fst_sliding2, "FST_SlidingWindow_RESULTS_185K_norCan_maf001.txt", col.names=T, row.names = F, quote = F, sep="\t")


############# Plot manhattan of mean values per windows

#Ensure chromosome is numeric
fst_sliding2$CHR=as.numeric(as.character(fst_sliding2$CHR))

#Check min and max values for FST plot
min(fst_sliding2$FST_mean)
max(fst_sliding2$FST_mean)

#Remove windows with less than 5 SNPs (only keep those with =5 or more)
fst_sliding2_5snpsplus=fst_sliding2[which(fst_sliding2$FST_n >= 5),]

#Plot manhattan of values
manhattan(x = fst_sliding2_5snpsplus, chr = "CHR", bp = "win_start", snp = "snp",
          p = "FST_mean", logp = F, ylim=c(0,0.7), col=c("gray50", "dodgerblue4"))

#Identify outliers (script taken from Brenna Forester's RDA vignette; https://popgen.nescent.org/2018-03-27_RDA_GEA.html)

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#Get Fst values
loadings=as.data.frame(fst_sliding2_5snpsplus$FST_mean)
#Add info chr/position (called snp) as row names
rownames(loadings)=fst_sliding2_5snpsplus$snp

#Get fst vlaues that are 3 standard deviations from the mean
cand_fst=outliers(loadings[,1], 3)

#Get minimum value of the outliers - to use for subsetting dataset
minimum_Fst_outlier=min(cand_fst)

#List of windows with FST outliers
fst_outliers=fst_sliding2_5snpsplus[which(fst_sliding2_5snpsplus$FST_mean>=minimum_Fst_outlier),]

#Save outlier windows to compare with Pcadapt
write.table(fst_outliers, "FST_SlidingWindow_Outlier_Windows_185K.txt", quote = F, row.names = F, col.names = T, sep="\t")

#Save manhattan as 12 x 6 inch pdf
#Highlight outlier values
#Plot for Figure 2
manhattan(x = fst_sliding2_5snpsplus, chr = "CHR", bp = "win_start", snp = "snp",
          p = "FST_mean", logp = F, ylim=c(0,0.7), # suggestiveline = 0.529, 
          highlight =fst_outliers$snp, ylab="Mean FST sliding window", 
          col=c("gray50", "dodgerblue4"))


#Chromosome specific plots - smoothed line - similar to plots in Figure 3 (see other script for combining all plots for Fig. 3)

chr6_plot=ggplot(fst_sliding2_5snpsplus[which(fst_sliding2_5snpsplus$CHR==6),])+#geom_point(aes(x=win_start, y=FST_mean))+
  geom_rect(aes(xmin=7e+07, xmax= 7.3e+07, ymin=0, ymax=0.6), fill="gray80")+
  geom_smooth(aes(x=win_mid, y=FST_mean), stat = "identity")+theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 6 position")+ylab("Mean FST")

chr13_plot=ggplot(fst_sliding2_5snpsplus[which(fst_sliding2_5snpsplus$CHR==13),])+#geom_point(aes(x=win_start, y=FST_mean))+
  geom_rect(aes(xmin=7.5e+07, xmax= 7.8e+07, ymin=0, ymax=0.6), fill="gray80")+
  geom_smooth(aes(x=win_mid, y=FST_mean), stat = "identity")+theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 13 position")+ylab("Mean FST")


chr16_plot=ggplot(fst_sliding2_5snpsplus[which(fst_sliding2_5snpsplus$CHR==16),])+#geom_point(aes(x=win_start, y=FST_mean))+
  geom_rect(aes(xmin=3.6e+07, xmax= 3.7e+07, ymin=0, ymax=0.6), fill="gray80")+
    geom_smooth(aes(x=win_mid, y=FST_mean), stat = "identity")+theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 16 position")+ylab("Mean FST")


chr19_plot=ggplot(fst_sliding2_5snpsplus[which(fst_sliding2_5snpsplus$CHR==19),])+#geom_point(aes(x=win_start, y=FST_mean))+
  geom_rect(aes(xmin=5.4e+07, xmax=5.5e+07, ymin=0, ymax=0.6), fill="gray80")+
  geom_smooth(aes(x=win_mid, y=FST_mean), stat = "identity")+theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 19 position")+ylab("Mean FST")


plot_grid(chr6_plot, chr13_plot, chr16_plot, chr19_plot,nrow = 1)


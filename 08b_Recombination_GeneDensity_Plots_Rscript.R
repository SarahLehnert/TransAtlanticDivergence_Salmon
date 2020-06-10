setwd("~/Desktop/Sarah/Salmon/TransAtlanticPaper/Recombination/")
library(genepopedit)
library(windowscanr)
library(ggplot2)
library(MareyMap)

#This R script includes
# 1 - Plots for recombination rate across chromosomes (note any misplaced SNPs were removed to calculate recombination rate)
#outlier regions were determined as values that were +/- 2 SD from the mean
# 2 - Gene density (or count) per 1MBP were calculated using windowscanr
#outlier regions were determined as values that were +/- 2 SD from the mean
# 3 - Plots for gene density were made and combined with those for recombination using cowplot
# 4 - Also includes subsetting of files for MareyMap to calculate recombination rate, running MareyMap GUI, as well as reading in the results for plotting


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
#but misplaced snps were removed to show recombination rate across the chromosome (smoothed line plot) - manually removed in MareyMap GUI (see notes)
#Note that recombination rate is in cm/Mbp

############ SSA 06
#Europe data
EU_data_Chr6=EU_data[which(EU_data$Chromosome==6),]

#Some SNPs were removed initially as they obviously deviated from patterns -- although additional SNPs were manually removed in MareyMap upon close inspection

#remove 2 misplaced SNPs
EU_data_Chr6_2removed=EU_data_Chr6[which(EU_data_Chr6$SNP_ID!="ESTNV_34441_695" & EU_data_Chr6$SNP_ID!="ESTNV_24999_794" ),]

#Canada
ssa06_canada_bp=Canada_data[which(Canada_data$Chromosome==6),]
#remove extra snps (misplaced)
ssa06_canada_bp_removed=ssa06_canada_bp[which(ssa06_canada_bp$SNP_ID!= "ESTNV_24999_794" & 
                                                ssa06_canada_bp$SNP_ID!= "ESTNV_25601_615" &  
                                                ssa06_canada_bp$SNP_ID!= "ESTNV_30516_1089"  ),]

#save for marey map package...
canada_ssa06_for_mareymaplib<- ssa06_canada_bp_removed[,c(2,3,5,22,23)]
canada_ssa06_for_mareymaplib$species=rep("Salmo_NorthAM")
colnames(canada_ssa06_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
canada_ssa06_for_mareymaplib$phys=(canada_ssa06_for_mareymaplib$phys1+canada_ssa06_for_mareymaplib$phys2)/2
head(canada_ssa06_for_mareymaplib)

#europe too
EU_ssa06_for_mareymaplib<- EU_data_Chr6_2removed[,c(2,4,6,20,21)]
EU_ssa06_for_mareymaplib$species=rep("Salmo_EU")
colnames(EU_ssa06_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
EU_ssa06_for_mareymaplib$phys=(EU_ssa06_for_mareymaplib$phys1+EU_ssa06_for_mareymaplib$phys2)/2
head(EU_ssa06_for_mareymaplib)


#Export for "MareyMap" package -- and read in results for plotting 
#Note some additional markers were removed in Marey Map if they appeared to deviate from expected (this was done by removing "valid" mark)
write.table(canada_ssa06_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa06_Canada_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(EU_ssa06_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa06_Europe_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")

#MareyMap uses a GUI so recombination rate must be calculated outside of R codes
#To start MareyMap use:
startMareyMapGUI()
#Note that some loci were removed manually in MareyMap to avoid negative recombination rate
#Results from MareyMap can be exported and then read into R for plotting below

###

Merged_eu_candaa_chr06_Recom=ggplot(data=ssa06_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.0e7, ymax=0, ymin=125, xmin=7.3e7, fill="gray")+
  geom_point(data=ssa06_canada_bp,aes(x=(start+end)/2, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr6,aes(x=(start+end)/2, y=Female.map, col="Europe"))+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+ theme_bw()+
  ylab("CM")+xlab("Ssa06 Position")+  theme(panel.grid = element_blank(), legend.position = "top")


#Merged plot - canada europe cm/Mbp
#Recombination rate now calcaulted in MareyMap -- corrected from previous version...
marey_map_results_Ssa06<-  read.table("Ssa06_Canada_results_mareymap.txt", header=T)
#negative values were changed to 0 recombination rate 
marey_map_results_Ssa06[which(marey_map_results_Ssa06$loess<0),]
marey_map_results_Ssa06$loess[which(marey_map_results_Ssa06$loess<0)] <- 0
max(marey_map_results_Ssa06$loess, na.rm=T)

marey_map_results_Ssa06_eu<-  read.table("Ssa06_Europe_results_mareymap.txt", header=T)
#negative values were changed to 0 recombination rate 
marey_map_results_Ssa06_eu[which(marey_map_results_Ssa06_eu$loess<0),]
marey_map_results_Ssa06_eu$loess[which(marey_map_results_Ssa06_eu$loess<0)] <- 0
max(marey_map_results_Ssa06_eu$loess, na.rm=T)

#outliers
#Check number of outliers for Ssa06 Europe and North AM
outliers(marey_map_results_Ssa06_eu$loess[which(marey_map_results_Ssa06_eu$loess>=0)],2)
min(outliers(marey_map_results_Ssa06_eu$loess[which(marey_map_results_Ssa06_eu$loess>=0)],2))
oultiers_ssa06_eu<-marey_map_results_Ssa06_eu[which(marey_map_results_Ssa06_eu$loess>=min(outliers(marey_map_results_Ssa06_eu$loess[which(marey_map_results_Ssa06_eu$loess>=0)],2)) ),]


#North AM 
outliers(marey_map_results_Ssa06$loess[which(marey_map_results_Ssa06$loess>=0)],2)
min(outliers(marey_map_results_Ssa06$loess[which(marey_map_results_Ssa06$loess>=0)],2))
oultiers_ssa06_canada<-marey_map_results_Ssa06[which(marey_map_results_Ssa06$loess>=min(outliers(marey_map_results_Ssa06$loess[which(marey_map_results_Ssa06$loess>=0)],2)) ),]

#plot MareyMap results for recombination rate 
#Note all negative values were converted to 0 recombination
merged_ssa06_CMbp=ggplot(data=marey_map_results_Ssa06[which(marey_map_results_Ssa06$loess>=0),])+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.0e7, ymax=21, ymin=0, xmin=7.3e7, fill="gray")+
  geom_smooth(data=marey_map_results_Ssa06[which(marey_map_results_Ssa06$loess>=0),],
              aes(x=phys, y=loess, col="North Am"), stat="identity")+
  geom_smooth(data=marey_map_results_Ssa06_eu[which(marey_map_results_Ssa06_eu$loess>=0),],
            aes(x=phys, y=loess, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa06 Position")+theme_bw()+
  geom_text(data=oultiers_ssa06_canada, aes(x=phys, y=20.5, label="*", col="North Am"), size=8)+
  geom_text(data=oultiers_ssa06_eu, aes(x=phys, y=20.1, label="*", col="Europe"), size=8)+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,21), expand = F)+
  NULL


########################################################
################ Ssa13
########################################################

EU_data_Chr13=EU_data[which(EU_data$Chromosome==13),]
ssa13_canada_bp=Canada_data[which(Canada_data$Chromosome==13),]

#remove misplaced SNPs
ssa13_canada_bp_removed=ssa13_canada_bp[which(ssa13_canada_bp$SNP_ID!= "ESTNV_32655_334" ),]



#save for marey map package...
canada_ssa13_for_mareymaplib<- ssa13_canada_bp_removed[,c(2,3,5,22,23)]
canada_ssa13_for_mareymaplib$species=rep("Salmo_NorthAM")
colnames(canada_ssa13_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
canada_ssa13_for_mareymaplib$phys=(canada_ssa13_for_mareymaplib$phys1+canada_ssa13_for_mareymaplib$phys2)/2
head(canada_ssa13_for_mareymaplib)

#europe too
EU_ssa13_for_mareymaplib<- EU_data_Chr13[,c(2,4,6,20,21)]
EU_ssa13_for_mareymaplib$species=rep("Salmo_EU")
colnames(EU_ssa13_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
EU_ssa13_for_mareymaplib$phys=(EU_ssa13_for_mareymaplib$phys1+EU_ssa13_for_mareymaplib$phys2)/2
head(EU_ssa13_for_mareymaplib)


#Export for "MareyMap" package -- and read in results for plotting 
#Note some additional markers were removed in Marey Map if they appeared to deviate from expected (this was done by removing "valid" mark)
write.table(canada_ssa13_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa13_Canada_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(EU_ssa13_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa13_Europe_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")

#MareyMap uses a GUI so recombination rate must be calculated outside of R codes
#To start MareyMap use:
startMareyMapGUI()
#Note that some loci were removed manually in MareyMap to avoid negative recombination rate
#Results from MareyMap can be exported and then read into R for plotting below


Merged_eu_candaa_chr13_Recom=ggplot(data=ssa13_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.5e7, ymax=0, ymin=125, xmin=7.8e7, fill="gray")+
  geom_point(data=ssa13_canada_bp,aes(x=(start+end)/2, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr13,aes(x=(start+end)/2, y=Female.map, col="Europe"))+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+theme_bw()+
  ylab("CM")+xlab("Ssa13 Position")+  theme(panel.grid = element_blank(), legend.position = "top")


#Merged plot - canada europe cm/Mbp
#Recombination rate now calcaulted in MareyMap -- corrected from previous version...
marey_map_results_Ssa13 <-  read.table("Ssa13_Canada_results_mareymap.txt", header=T)
#Negative values were changed to 0
marey_map_results_Ssa13[which(marey_map_results_Ssa13$loess<0),]
marey_map_results_Ssa13$loess[which(marey_map_results_Ssa13$loess<0)] <- 0
max(marey_map_results_Ssa13$loess, na.rm=T)

marey_map_results_Ssa13_eu<-  read.table("Ssa13_Europe_results_mareymap.txt", header=T)
#Negative values were changed to 0
marey_map_results_Ssa13_eu[which(marey_map_results_Ssa13_eu$loess<0),]
marey_map_results_Ssa13_eu$loess[which(marey_map_results_Ssa13_eu$loess<0)] <- 0
max(marey_map_results_Ssa13_eu$loess, na.rm=T)

#outliers
#Check number of outliers for Ssa13 Europe and North AM
outliers(marey_map_results_Ssa13_eu$loess[which(marey_map_results_Ssa13_eu$loess>=0)],2)
min(outliers(marey_map_results_Ssa13_eu$loess[which(marey_map_results_Ssa13_eu$loess>=0)],2))
oultiers_ssa13_eu<-marey_map_results_Ssa13_eu[which(marey_map_results_Ssa13_eu$loess>=min(outliers(marey_map_results_Ssa13_eu$loess[which(marey_map_results_Ssa13_eu$loess>=0)],2)) ),]



#North AM - just 0 cM/mbp is outlier (two values)
outliers(marey_map_results_Ssa13$loess[which(marey_map_results_Ssa13$loess>=0)],2)
min(outliers(marey_map_results_Ssa13$loess[which(marey_map_results_Ssa13$loess>=0)],2))
oultiers_ssa13_canada<-marey_map_results_Ssa13[which(marey_map_results_Ssa13$loess>=min(outliers(marey_map_results_Ssa13$loess[which(marey_map_results_Ssa13$loess>=0)],2)) ),]


#with negative values set to 0 
merged_ssa13_CMbp=ggplot(data=marey_map_results_Ssa13[which(marey_map_results_Ssa13$loess>=0),])+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=7.5e7, ymax=6, ymin=0, xmin=7.8e7, fill="gray")+
  geom_smooth(data=marey_map_results_Ssa13[which(marey_map_results_Ssa13$loess>=0),],
              aes(x=phys, y=loess, col="North Am"), stat="identity")+
  geom_smooth(data=marey_map_results_Ssa13_eu[which(marey_map_results_Ssa13_eu$loess>=0),],
              aes(x=phys, y=loess, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa13 Position")+theme_bw()+
  geom_text(data=oultiers_ssa13_canada, aes(x=phys, y=4.8, label="*", col="North Am"), size=8)+
  geom_text(data=oultiers_ssa13_eu, aes(x=phys, y=4.5, label="*", col="Europe"), size=8)+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,6), expand = F)+
  NULL


########################################################
################ Ssa16
########################################################

EU_data_Chr16=EU_data[which(EU_data$Chromosome==16),]
ssa16_canada_bp=Canada_data[which(Canada_data$Chromosome==16),]

#save for marey map package...
canada_ssa16_for_mareymaplib<- ssa16_canada_bp[,c(2,3,5,22,23)]
canada_ssa16_for_mareymaplib$species=rep("Salmo_NorthAM")
colnames(canada_ssa16_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
canada_ssa16_for_mareymaplib$phys=(canada_ssa16_for_mareymaplib$phys1+canada_ssa16_for_mareymaplib$phys2)/2
head(canada_ssa16_for_mareymaplib)

#europe too
EU_ssa16_for_mareymaplib<- EU_data_Chr16[,c(2,4,6,20,21)]
EU_ssa16_for_mareymaplib$species=rep("Salmo_EU")
colnames(EU_ssa16_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
EU_ssa16_for_mareymaplib$phys=(EU_ssa16_for_mareymaplib$phys1+EU_ssa16_for_mareymaplib$phys2)/2
head(EU_ssa16_for_mareymaplib)


#Export for "MareyMap" package -- and read in results for plotting 
#Note some additional markers were removed in Marey Map if they appeared to deviate from expected (this was done by removing "valid" mark)
write.table(canada_ssa16_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa16_Canada_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(EU_ssa16_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa16_Europe_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")

#MareyMap uses a GUI so recombination rate must be calculated outside of R codes
#To start MareyMap use:
startMareyMapGUI()
#Note that some loci were removed manually in MareyMap to avoid negative recombination rate
#Results from MareyMap can be exported and then read into R for plotting below

#Plot physical x genetic map position
Merged_eu_candaa_chr16_Recom=ggplot(data=ssa16_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=3.6e7, ymax=0, ymin=125, xmin=3.7e7, fill="gray")+
  geom_point(data=ssa16_canada_bp,aes(x=(start+end)/2, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr16,aes(x=(start+end)/2, y=Female.map, col="Europe"))+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+theme_bw()+
  ylab("CM")+xlab("Ssa16 Position")+  theme(panel.grid = element_blank(), legend.position = "top")



#Merged plot - canada europe cm/Mbp
#Recombination rate now calcaulted in MareyMap -- corrected from previous version...
marey_map_results_Ssa16 <-  read.table("Ssa16_Canada_results_mareymap.txt", header=T)
#Negative values were changed to 0
marey_map_results_Ssa16[which(marey_map_results_Ssa16$loess<0),]
marey_map_results_Ssa16$loess[which(marey_map_results_Ssa16$loess<0)] <- 0
max(marey_map_results_Ssa16$loess, na.rm=T)

marey_map_results_Ssa16_eu<-  read.table("Ssa16_Europe_results_mareymap.txt", header=T)
#Negative values were changed to 0
marey_map_results_Ssa16_eu[which(marey_map_results_Ssa16_eu$loess<0),]
marey_map_results_Ssa16_eu$loess[which(marey_map_results_Ssa16_eu$loess<0)] <- 0
max(marey_map_results_Ssa16_eu$loess, na.rm=T)

#outliers
#Check number of outliers for Ssa16 Europe and North AM
outliers(marey_map_results_Ssa16_eu$loess[which(marey_map_results_Ssa16_eu$loess>=0)],2)
min(outliers(marey_map_results_Ssa16_eu$loess[which(marey_map_results_Ssa16_eu$loess>=0)],2))
oultiers_ssa16_eu<-marey_map_results_Ssa16_eu[which(marey_map_results_Ssa16_eu$loess>=min(outliers(marey_map_results_Ssa16_eu$loess[which(marey_map_results_Ssa16_eu$loess>=0)],2)) ),]

#North AM - just 0 cM/mbp is outlier (two values)
outliers(marey_map_results_Ssa16$loess[which(marey_map_results_Ssa16$loess>=0)],2)
min(outliers(marey_map_results_Ssa16$loess[which(marey_map_results_Ssa16$loess>=0)],2))
oultiers_ssa16_canada<-marey_map_results_Ssa16[which(marey_map_results_Ssa16$loess>=min(outliers(marey_map_results_Ssa16$loess[which(marey_map_results_Ssa16$loess>=0)],2)) ),]


merged_ssa16_CMbp=ggplot(data=marey_map_results_Ssa16[which(marey_map_results_Ssa16$loess>=0),])+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=3.6e7, ymax=7, ymin=0, xmin=3.7e7, fill="gray")+
  geom_smooth(data=marey_map_results_Ssa16[which(marey_map_results_Ssa16$loess>=0),],
              aes(x=phys, y=loess, col="North Am"), stat="identity")+
  geom_smooth(data=marey_map_results_Ssa16_eu[which(marey_map_results_Ssa16_eu$loess>=0),],
              aes(x=phys, y=loess, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa16 Position")+theme_bw()+
  geom_text(data=oultiers_ssa16_canada, aes(x=phys, y=6.3, label="*", col="North Am"), size=8)+
  geom_text(data=oultiers_ssa16_eu, aes(x=phys, y=6, label="*", col="Europe"), size=8)+
  #geom_point(data= map_bp[which(map_bp$V2!="ssa06"  & map_bp$Chromosome==6),], aes(x=start, y=Female.map, col="dodgerblue3"))+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,7), expand = F)+
  NULL



########################################################
################ Ssa19
########################################################

EU_data_Chr19=EU_data[which(EU_data$Chromosome==19),]

##Remove 1 misplaced SNP
EU_data_Chr19_remove1=EU_data_Chr19[which(EU_data_Chr19$SNP_ID!="ESTV_16749_913"),]

ssa19_canada_bp=Canada_data[which(Canada_data$Chromosome==19),]
#remove 1 snps misplaced
ssa19_canada_bp_remove1=ssa19_canada_bp[which(ssa19_canada_bp$SNP_ID!="ESTNV_32145_342"),]

#save for marey map package...
canada_ssa19_for_mareymaplib<- ssa19_canada_bp_remove1[,c(2,3,5,22,23)]
canada_ssa19_for_mareymaplib$species=rep("Salmo_NorthAM")
colnames(canada_ssa19_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
canada_ssa19_for_mareymaplib$phys=(canada_ssa19_for_mareymaplib$phys1+canada_ssa19_for_mareymaplib$phys2)/2
head(canada_ssa19_for_mareymaplib)

#europe too
EU_ssa19_for_mareymaplib<- EU_data_Chr19_remove1[,c(2,4,6,20,21)]
EU_ssa19_for_mareymaplib$species=rep("Salmo_EU")
colnames(EU_ssa19_for_mareymaplib)=c("mkr", "map", "gen", "phys1", "phys2", "set")
#take average of start and end position physical
EU_ssa19_for_mareymaplib$phys=(EU_ssa19_for_mareymaplib$phys1+EU_ssa19_for_mareymaplib$phys2)/2
head(EU_ssa19_for_mareymaplib)


#Export for "MareyMap" package -- and read in results for plotting 
#Note some additional markers were removed in Marey Map if they appeared to deviate from expected (this was done by removing "valid" mark)
write.table(canada_ssa19_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa19_Canada_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")
write.table(EU_ssa19_for_mareymaplib[,c("set", "map", "mkr", "phys", "gen")], "Ssa19_Europe_forMareyMap.txt", quote = F, row.names = F, col.names = T, sep="\t")


#MareyMap uses a GUI so recombination rate must be calculated outside of R codes
#To start MareyMap use:
startMareyMapGUI()
#Note that some loci were removed manually in MareyMap to avoid negative recombination rate
#Results from MareyMap can be exported and then read into R for plotting below

#Merge Position / CM
merged_ssa19_CM_Position_noremove=ggplot(data=ssa19_canada_bp)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=5.4e7, ymax=0, ymin=125, xmin=5.5e7, fill="gray")+
  geom_point(data=ssa19_canada_bp,aes(x=(start+end)/2, y=Female.map, col="North Am"))+
  geom_point(data=EU_data_Chr19,aes(x=(start+end)/2, y=Female.map,  col="Europe"))+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+theme_bw()+
  ylab("CM")+xlab("Ssa19 Position")+  theme(panel.grid = element_blank(), legend.position = "top")





#Merged plot - canada europe cm/Mbp
#Recombination rate now calcaulted in MareyMap -- corrected from previous version...
marey_map_results_Ssa19 <-  read.table("Ssa19_Canada_results_mareymap.txt", header=T)
marey_map_results_Ssa19[which(marey_map_results_Ssa19$loess<0),]
marey_map_results_Ssa19$loess[which(marey_map_results_Ssa19$loess<0)] <- 0
max(marey_map_results_Ssa19$loess, na.rm=T)

marey_map_results_Ssa19_eu<-  read.table("Ssa19_Europe_results_mareymap.txt", header=T)
marey_map_results_Ssa19_eu[which(marey_map_results_Ssa19_eu$loess<0),]
marey_map_results_Ssa19_eu$loess[which(marey_map_results_Ssa19_eu$loess<0)] <- 0
max(marey_map_results_Ssa19_eu$loess, na.rm=T)

#outliers
#Check number of outliers for Ssa19 Europe and North AM
outliers(marey_map_results_Ssa19_eu$loess[which(marey_map_results_Ssa19_eu$loess>=0)],2)
min(outliers(marey_map_results_Ssa19_eu$loess[which(marey_map_results_Ssa19_eu$loess>=0)],2))
oultiers_ssa19_eu<-marey_map_results_Ssa19_eu[which(marey_map_results_Ssa19_eu$loess>=min(outliers(marey_map_results_Ssa19_eu$loess[which(marey_map_results_Ssa19_eu$loess>=0)],2)) ),]


#North AM - just 0 cM/mbp is outlier (two values)
outliers(marey_map_results_Ssa19$loess[which(marey_map_results_Ssa19$loess>=0)],2)
min(outliers(marey_map_results_Ssa19$loess[which(marey_map_results_Ssa19$loess>=0)],2))
oultiers_ssa19_canada<-marey_map_results_Ssa19[which(marey_map_results_Ssa19$loess>=min(outliers(marey_map_results_Ssa19$loess[which(marey_map_results_Ssa19$loess>=0)],2)) ),]


merged_ssa19_CMbp_remove=ggplot(data=ssa19_canada_bp_remove1)+#scale_y_continuous(limits = c(0, 10), expand = c(0,0))+
  geom_rect(xmax=5.4e7, ymax=0, ymin=9, xmin=5.5e7, fill="gray")+
  geom_smooth(data=marey_map_results_Ssa19[which(marey_map_results_Ssa19$loess>=0),],
              aes(x=phys, y=loess, col="North Am"), stat="identity")+
  geom_smooth(data=marey_map_results_Ssa19_eu[which(marey_map_results_Ssa19_eu$loess>=0),],
              aes(x=phys, y=loess, col="Europe"), stat="identity")+
  scale_color_manual(values=c("dodgerblue3", "indianred2"))+
  ylab("cM/Mbp")+xlab("Ssa19 Position")+theme_bw()+
  geom_text(data=oultiers_ssa19_canada, aes(x=phys, y=8.4, label="*", col="North Am"), size=8)+
  geom_text(data=oultiers_ssa19_eu, aes(x=phys, y=8, label="*", col="Europe"), size=8)+
  theme(panel.grid = element_blank(), legend.position = "top")+ coord_cartesian(ylim=c(0,9), expand = F)+
  NULL


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






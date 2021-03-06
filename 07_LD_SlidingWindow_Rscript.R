##Rscript for LD across chromosomes

#set directory
setwd("/Users/Bradbury/Desktop/Sarah/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/AlleleFreq/Whole_Chromosome_LD/")

#Packages
library(genepopedit)
library(reshape2)
library(ggplot2)
library(gplots)

#This R script include 
  # 1 - Reading in r2 matrix from plink for the whole chromosome (continent seperately)
  # 2 - A for loop that calculated mean r2 (LD) between SNPs within a 50 SNP window - to calculate mean LD across the chromosome
#This script is for chromosome 6 but same was applied to other chromosomes of interest

###########
##Sliding window LD

#LD matrix was created in plink using '--r2 square' to create the square matrix of LD between all snps on the chromosome
#LD was calculated seperately for Europe and Canada

#Read LD matrix - one matrix for full chromosome
LD_matrix=data.table::fread("Chr6_norway_LD.ld", header=F)
dim(LD_matrix)

#read Map positions for SNPs - map file from plink
map=data.table::fread("Chr6_nor_can.map")
colnames(LD_matrix)=map$V2 #change column names of matrix to SNP names

#Set as matrix
LD_matrix=as.matrix(LD_matrix) 

#Change diagnoal to NAs
diag(LD_matrix)<-NA

#Start and end of each matrix to calculate mean for 50 SNP at each snp position
#Run with full matrix  (make sure only diagonal is NA; ie. not missing upper or lower matrix values)
#Start at column 25 and get mean pairwise LD between the previous 25 and next 25 SNPs  (50 SNPs)

result<-NULL

#Run loop
for(i in 1:(nrow(LD_matrix)-50)){
  row_end=i+50
  column_start=i+25
  result[i]<-mean(LD_matrix[i:row_end,column_start], na.rm=T)
}

#Check results
print(result)
str(LD_matrix)
str(result) #Results will be shorter than LD matrix by 50 rows


y=nrow(LD_matrix)-25 #Get number of rows in LD matrix without last 25 (remove last 25 SNPs)
loci=map$V2[26:y] #Get SNP IDs starting at SNP 25 and going to end with last 25 removed (remove first 25 SNPs + last 25 SNPs)
all_LD_15=cbind(result, loci) #cbind resutls from matrix calculations with Loci info

all_LD_15 #Has Mean LD (for window) with SNP name

#Merge loci LD/name with Locus position for plot
data_for_plot=merge(all_LD_15, map, by.x=2, by.y=2)
head(data_for_plot)

#Data for plot
data_for_plot=data_for_plot[,c(1,2,3,5)]
colnames(data_for_plot)=c("Locus", "Mean_LD", "Chr",  "Position")

#Change to numeric
data_for_plot$Position=as.numeric(as.character(data_for_plot$Position))
data_for_plot$Mean_LD=as.numeric(as.character(data_for_plot$Mean_LD))

#Quick plot
plot(data_for_plot$Position,data_for_plot$Mean_LD) #NOTE** see plot for nicer ggplot script BELOW

#Save as file and replicate same analyses for Canadian populations
write.table(data_for_plot, file="LD_50nearbySNPWindow_Chr6_NORWAY_185K.txt",quote = F, col.names = T, row.names = F)

#Same script for Canada populations:
#Read plink results Pop2 (LD matrix)
LD_matrix_can=data.table::fread("Chr6_canada_LD.ld", header=F)

#Same map (same SNPs were used for each here)
map[1:10,]
colnames(LD_matrix_can)=map$V2

LD_matrix_can=as.matrix(LD_matrix_can)
diag(LD_matrix_can)<-NA

#Start and end of each matrix to calculate mean

#Run with full matrix 
#Start at column 25 and get mean pairwise LD between the previous 25 and next 25 SNPs 
#Not NAs are ignored for calculations (note for fixed loci, LD isn't calculated)
result_can<-NULL
for(i in 1:(nrow(LD_matrix_can)-50)){
  row_end=i+50
  column_start=i+25
  result[i]<-mean(LD_matrix_can[i:row_end,column_start], na.rm=T)
}
print(result_can)

#can use same map as europe
y=nrow(LD_matrix_can)-25
loci=map$V2[26:y]

all_LD_15_can=cbind(result_can, loci)

#Merge loci LD/name with Locus position for plot

#Merge loci LD/name with Locus position for plot
data_for_plot_can=merge(all_LD_15_can, map, by.x=2, by.y=2)
head(data_for_plot_can)
data_for_plot_can=data_for_plot_can[,c(1,2,3,5)]
colnames(data_for_plot_can)=c("Locus", "Mean_LD", "Chr",  "Position")
data_for_plot_can$Position=as.numeric(as.character(data_for_plot_can$Position))
data_for_plot_can$Mean_LD=as.numeric(as.character(data_for_plot_can$Mean_LD))

plot(data_for_plot_can$Position,data_for_plot_can$Mean_LD)

#save canadian data
write.table(data_for_plot_can, file="LD_50nearbySNPWindow_Chr6_CANADA_maf001.txt",quote = F, col.names = T, row.names = F)

#See other script for plotting results

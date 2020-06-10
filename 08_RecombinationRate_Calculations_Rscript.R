#Open libraries
library(seqRFLP)
library(dplyr)

setwd("~/Desktop/Sarah/Salmon/TransAtlanticPaper/Recombination/")

#This R script includes
# 1 - Create a fasta file of available SNPs from the atlantic salmon linkage map 
#    -This fasta file was then blasted against the genome (not included here) and a .csv of hits was saved
# 2 - Results of BLAST were then filtered here based on e-value and identity to select the best match
# 3 - Physical genomic position and linkage map were compared for available SNPs (for Europe and North America seperately)
# 4 - Recombination rates were previously calculated as cm/MBP using linkage map and physical position
#  --> but this has now been calculated using MareyMap package (see additional R script for details)
#  --> this was done after manuscript appeared online but before final publication (corrected in proofs)

#First use published sequences from linkage map to get genomic positions (physical) 

setwd("~/Desktop/Sarah/Salmon/TransAtlanticPaper/Recombination/")

#First use published sequences from linkage map to get genomic positions (physical) 

#Read linkage map info - get sequences for SNPs
df=read.csv("European_LinkageMap.csv", header=T)
df2=df[,c(10,8)]
#Create fasta file
df.fasta = dataframe2fas(df2, file="df.fasta")
#Save as fasta file
writeLines(df.fasta, "EU_map_FASTA.fasta") 
#Next used BLAST online to get position of all sequences in genome of Atlantic salmon

#Read in North American linkage map (same SNPs - but less) -
df_can=read.csv("NorthAmericaLinkageMap2.csv", header=T)

#######After running sequences through BLAST - read in results
hits=read.csv("EuropeanMapAlignment-HitTable.csv", header=F)
colnames(hits)=c("SNP", "CHR", "Iden", "size", "x", "y", "QueryStart", "QueryEnd", "start", "end", "Eval", "MaxScore")

#Keep only hits for each snp with the maximum vscore - and then lowest evalue
hit_results=hits %>% group_by(SNP) %>% top_n(1, MaxScore)
hit_results2=hit_results %>% group_by(SNP) %>% top_n(1, 1/Eval)

#Note some SNPs are duplicated (Remove these for analyses below)
#means that evalues and max scores are the same (just keep one - no reason to select one of another)
sum(duplicated(hit_results$SNP))

#Get chromosome names to merge with data from NCBI
Chrom_names=read.table("Chr_report.txt")

#Set as data frame for Hits
hit_results2=as.data.frame(hit_results)

#Merge hit info with chromosome numbers
hits_all=merge(x=Chrom_names, y=hit_results2, by=4, by.y=2)

#Create reduced data frame (less columns) - snp info (original dataframe)
df5=df[,c(10,1:7)]

#Make sure strings are characters if necessary
df5=as.data.frame(df5)
hits_all=as.data.frame(hits_all)
hits_all$SNP=as.character((hits_all$SNP))
df5$dbSNP.access.=as.character((df5$dbSNP.access.))

#Check number of duplicate SNPs included (remove these below)
nrow(hits_all[duplicated(hits_all$SNP),])
nrow((hits_all[!duplicated(hits_all$SNP),]))

#Merge SNP hits wiht map infor
map_bp=merge(x=df5, y=hits_all, by.x=1, by.y=6)
head(map_bp)

#Merge canadian map with SNP info hit
map_canada_Bp=merge(x=df_can, y=hits_all, by.x=7, by.y=6)

#########
######### Note --- It has come to my attention that this is not the appropraite way to calculate recombination rate
### Therefore - recombination was calculated using the MareyMap R package - see other script for more details... 

#############################################################################################################################################
######### These estimates of cM/MBP were not used in the paper - but note other parts of the script are relevant ############################
#############################################################################################################################################

#Create cM/MBP value -- note it is actually cm/MBP (but R heading was CM_BP) 
#THIS WAS NOT USED: - 
#This is not the appropriate way to calculate recombination rate (now use MareyMap package in R -see other R script)
#map_bp$CM_BP=map_bp$Female.map/((map_bp$start/1e6))
#map_canada_Bp$CM_BP=map_canada_Bp$Female.map/(map_canada_Bp$start/1e6)

#Remove duplicated SNPs for Europe
test=map_bp[!duplicated(map_bp$SNP), ]
nrow(test)-sum(is.na(test$Female.map)) #Really for female map there are 5422 SNPs total (no duplicated)

#Remove duplicated SNPs for north America
test2=map_canada_Bp[!duplicated(map_canada_Bp$SNP), ]
nrow(test2)-sum(is.na(test2$Female.map)) #Really for female canadian map there are 3080 SNPs total (no duplicated)

#Calculated Male values - but did not use in paper - this is cm/MBP
#This is not the appropriate way to calculate recombination rate (now use MareyMap package in R -see other R script)
#map_bp$CM_BP_male=map_bp$Male.map/((map_bp$start/1e6))
#map_canada_Bp$CM_BP_male=map_canada_Bp$Male.map/(map_canada_Bp$start/1e6)

#Add column for chromosome blast - to match numeric values of Ssa chromosome
map_canada_Bp$Chr_blast=map_canada_Bp$V2
map_canada_Bp$Chr_blast=as.numeric(map_canada_Bp$Chr_blast)

#Based on linkage map information -- keep only SNPs where the Chr from blast matches the chromosome from linkage map
map_correct_Chr=map_canada_Bp[which(map_canada_Bp$Chromosome==map_canada_Bp$Chr_blast),]
#remove duplicated SNPs
map_correct_Chr=map_correct_Chr[!duplicated(map_correct_Chr$SNP_ID),]

#Save results for North America
write.csv(map_correct_Chr, "RecombinationRate_allSNPs_Chr_Canada.csv", quote = F, row.names = F, col.names = T, sep="\t")

#Do same for European Linkage Map/blast info
map_bp$Chr_blast=map_bp$V2
map_bp$Chr_blast=as.numeric(map_bp$Chr_blast)

#Based on linkage map information -- keep only SNPs where the Chr from blast matches the chromosome from linkage map
map_correct_Chr_EU=map_bp[which(map_bp$Chromosome==map_bp$Chr_blast),]
map_correct_Chr_EU=map_correct_Chr_EU[!duplicated(map_correct_Chr_EU$SNP_ID),]

#Save results for Europe
write.csv(map_correct_Chr_EU[!is.na(map_correct_Chr_EU$Female.map),], "RecombinationRate_allSNPs_Chr_Europe.csv", quote = F, row.names = F, col.names = T, sep="\t")


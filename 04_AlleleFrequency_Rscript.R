setwd("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/")

#Packages
library(genepopedit)
library(reshape2)
library(gplots)

#This R script include 
 # 1- Generating allele frequencies for populations using genepopedit
 # 2- Adjusting frequencies so that the allele that is considered the 'major' allele is found in Canada (results are the same but this is used for plotting purposes)
 # 3- Plotting allele frequencies in a heatmap with populations aligned by latitude within continents
#This script was used for Ssa06, and the same script was applied to other chromosomes of interest


#Get MAF across Chr 6 region - from subsetted genepop file
#This should be done in plink for larger datasets
freq=genepop_allelefreq("AlleleFreq/AlleleFreq_region/Chr6/LargerRegion/genepop_chr6_67-76.txt", wide = F, fullpanel = F)

#Transpose data into wide format
freq_test=reshape(freq, timevar="Loci", idvar=c("Population"), direction="wide")

#Make dataframe
freq_test=as.data.frame(freq_test)

#Add Pop info to dataframe
pop_info=read.table("AlleleFreq/Poplist.txt", header=F)
freq_test_new=merge(x=pop_info, by.x=1, y=freq_test, by.y=1)

#Order sites by country
freq_test2=freq_test_new[order(freq_test_new$V3),]

#Add latitude data to allele frequency table
loc=read.csv("AlleleFreq/All_80_Locations.csv", header=T)
new_merged_freq=merge(y=freq_test2, by.y=1, x=loc, by.x=2)
freq_by_lat=new_merged_freq[order(new_merged_freq$Lat, decreasing=T),]

#Check dim and change last column to final column of data frame here
dim(freq_by_lat)
max_col=ncol(freq_by_lat)
freq_latitude_3=freq_by_lat[,c(1,3, 11:max_col)] #Inlcude Pop, Lat, and MAF for each locus

#change values so alleles are changed based on frequency in Canada (if maf < 0.5 change all values to 1-maf, if maf >0.5 leave it as is..)
#Allele frequencies in Canada should be closer to 1, so some SNPs had to be changed if it was closer to 1 in Europe
#Just used for plotting purposes to show differences between continents in regions of interest.
for(i in 2:ncol(freq_latitude_3)){
  
  if(freq_latitude_3[80,i]< 0.5){
    freq_latitude_3[,i]=1-freq_latitude_3[,i]
  }
  if(freq_latitude_3[80,i]> 0.5){
    freq_latitude_3[,i]=freq_latitude_3[,i]
  }
}

#Remove Pop and Latitude info and just use MAF for plotting allele frequency heatmap
data.table::fwrite(freq_latitude_3, "Allelefreqs_chr6_region_by_LAT_larger67-76MbpREgion.txt", col.names = T, row.names = F, quote = F, sep="\t")

#To use for plotting
ncol_use=ncol(freq_latitude_3)
freq_test_new_3=as.matrix(freq_latitude_3[,3:ncol_use])

#Read in map positions for the SNPs
map=read.table("AlleleFreq/AlleleFreq_region/Chr6/LargerRegion/chr6_goodLoc_67-76.map")

#Set colour of all snp positions to "gray"
map$col=rep("gray")
#Then change region of interest to 'white' so this area is highlighted in the heatmap
map$col[which(map$V4>70000000 & map$V4<73000000)]<-"white"

#Save heatmap
mypalette4<-colorRampPalette(c("yellow","purple"))

pdf(file="AlleleFreq_across_pops_Chr6_1MBP_region_67-76_outside_NEW.pdf", width = 13, height=12, bg = "white")
gplots::heatmap.2(freq_test_new_3,
                  Rowv=FALSE, #rows should be reordered as required
                  Colv = "Rowv", #columns should be treated as rows
                  dendrogram="none", #no trees
                  scale="none",
                  breaks=100, #number of break points used to bin into colours
                  col=mypalette4,
                  trace="none", #whether lines should be drawn between cols, rows,
                  margins=c(5,5),#margins for column names and row names
                  labRow= freq_latitude_3$Genepop_Name,
                  labCol= " ",
                  ColSideColors = as.character(map$col),
                  key=TRUE,
                  keysize = 1,
                  density.info="none")
dev.off()



#Packages
library(genepopedit)
library(reshape2)
library(gplots)

#Set directory
setwd("/Users/Bradbury/Desktop/Sarah/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/AlleleFreq/Region_LD/LargerRegion_Chr6_13/")

#This R script includes
  # 1 - Reading in r2 matrix (LD) generated in plink and plotting LD as a heatmap
#This is the script for Ssa06 - the same script was used for each chromosome


#LD matrix is created in plink using --r2 square to create a square matrix of LD values (r2)
#This is done seperately for Europe and Canada

##LD heatmap - read in LD matrix for Canada
mat_CAN=data.table::fread("chr6_canada_LD_70-73.ld", header=F)
mat_CAN2=as.matrix(mat_CAN)
dim(mat_CAN2)

#Read in LD matrix for Norway
mat_NOR=data.table::fread("chr6_norway_LD_70-73.ld", header=F)
mat_NOR2=as.matrix(mat_NOR)
dim(mat_NOR2)

#Combine matrix for Canada and Norway (make sure same number of SNPs are used for each - do not filter)

#Create matrix with upper tri as Pop1 and lower tri as Pop2
new<-mat_NOR2
diag(new) <- 1
new[upper.tri(new)] <- mat_CAN2[upper.tri(mat_CAN2)]
diag(new) <- NA

#colour palette
mypalette4<-colorRampPalette(c("aliceblue", "blue", "red"))

#Save ad PDF
pdf(file="LDmatrix_Ssa06_can_vs_norway.pdf", width = 13, height=12, bg = "white")
gplots::heatmap.2(new,
                  Rowv=FALSE, #rows should be reordered as required
                  Colv = "Rowv", #columns should be treated as rows
                  dendrogram="none", #no trees
                  scale="none",
                  breaks=100, #number of break points used to bin into colours
                  col=mypalette4,
                  trace="none", #whether lines should be drawn between cols, rows,
                  margins=c(5,5),#margins for column names and row names
                  labRow= " ",
                  labCol= " ",
                 #ColSideColors = as.character(map$col),
                  key=TRUE,
                  xlab="Norway",
                  ylab="Canada",
                  keysize = 1,
                  density.info="none"
                  #lmat = rbind(c(3,1),c(4,2)), #order of elements (1 is heatmap, 4 is key) #makes the plot a 2x2 grid, each "cell" has 1 element - row tree, key, col tree, and heatmap
                  #lhei = c(20,5), #row height for plot elements
                  #lwid = c(8,30)  #column width for plot elements)
)
dev.off()

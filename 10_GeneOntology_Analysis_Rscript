#Packages
library(topGO)

#Set working directory
setwd("~/Desktop/SalmoAnalysis/Transatlanticpaper/Revised_184K_SNPS/Top_1_FST/")

#This Rscript contains analysis for topGO
#This R script was accessed/modified from Wellband et al. 2019 Mol Ecol (see https://github.com/kylewellband/ssa_cast2016/blob/master/00_scripts/07_outlier_enrichment_analysis.R)
#Prior to running this script, BEDTOOLS was used to generate '.anno' files for the reference (all snps) and outlier datasets
#BEDTOOLS was used with Ssal_ICSASG_v2-18092015.gff3 from NCBI (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/233/375/GCA_000233375.4_ICSASG_v2/)

#all gene ontology for salmo
#Read in downloaded file with GO annotation for the Salmo Salar Genome (ICSASG_v2) - downloaded from salmobase (accessed: Sept 2018)
go.anno <- read.csv("~/Desktop/SalmoAnalysis/Ssal_ICSASG_v2_GOAccession.csv", header = T, stringsAsFactors = F, sep="\t")
go.anno[,1] <- gsub("\\.t[0-9]*", "", go.anno[,1])
go.anno=go.anno[,c(1,4)]

#Outliers (intersected) - top 1% of FST snps (1843 SNPs)
outliers <- read.table("TOP_1%SNPs_FST.accnos", stringsAsFactors = F)[,1]
outliers <- gsub("\\.t[0-9]*", "", outliers)

#Get outliers in annotated gene ontology
outliers.w.go <- outliers[outliers %in% go.anno[,1]]

#all genes -all SNPs used in fst analysis (184 K SNPs)
#For the reference dataset genes within 10Kbp of the 184,295 polymorphic SNPs were used
all.220K.genes <- read.table("ALLSNPS_FST.accnos", stringsAsFactors = F)[,1]
all.220K.genes.w.go <- all.220K.genes[all.220K.genes %in% go.anno[,1]]

#Create function for topGo
CreateGene2GOMapping <- function(x) {
  nr <- nrow(x)
  map <- list()
  for(r in 1:nr) {
    map[[as.character(x[r, 1])]] <- append(map[[as.character(x[r, 1])]], x[r,2])
  }
  return(map)
}

# create GO mapping, takes about 20 minutes
gene2GO.map <- CreateGene2GOMapping(go.anno)
gene2GO.map <- gene2GO.map[names(gene2GO.map) %in% all.220K.genes.w.go]
gene.list <- factor(as.integer(all.220K.genes.w.go %in% outliers.w.go))
names(gene.list) <- all.220K.genes.w.go

transAtl.GOdata <- new("topGOdata",
                      description = "trans atlantic outliers", ontology = "BP",
                      allGenes = gene.list,
                      nodeSize = 5,
                      annot = annFUN.gene2GO, gene2GO = gene2GO.map)

transAtl.result.Fisher <- runTest(transAtl.GOdata, algorithm = "weight01", statistic = "fisher")

geneData(transAtl.result.Fisher)
hist(score(transAtl.result.Fisher), 50, xlab = "p-values")

#Save results
write.table(GenTable(transAtl.GOdata, transAtl.result.Fisher, topNodes = length(transAtl.result.Fisher@score)), 
            "TransAtlantic_Top1%_SNPs_GO_184k.txt", quote = F, row.names = F , sep="\t", col.names = T)

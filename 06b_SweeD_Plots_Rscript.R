#Libraries 
library(ggplot2)

#Rscript for plotting CLR/Sweed results (Figure in Supplement)

#This R script include 
  # 1 - Reading in SweeD results and plotting them with outliers indicated by the top 5% of SNP values

#Set directory
setwd("/Volumes/Sarah Lehnert DFO/FROM_COMPUTERS/Windows_Computer/TransAtlantic_Genomewide_213K/185K_SNPS_Nor_Can_NEW/SweeD/")

##########  Ssa06

#Read in values for Europe and North America
sweeps6=read.table("SweeD_Report.chr6_sweed_europe", header=T, skip=2, fill=T)
sweeps_NA6=read.table("SweeD_Report.chr6_sweed_northam", header=T, skip=2, fill=T)

#Get number of SNPs in top 5% 
keep_eu=nrow(sweeps6)*0.05
keep_na=nrow(sweeps_NA6)*0.05

#Get top SNPs - considered signifcant as represent the top 5% of SNP values
top_na=sweeps_NA6[order(sweeps_NA6$Likelihood, decreasing = T),]
sweeps_sign_NA=top_na[1:keep_na,]
top_EU=sweeps6[order(sweeps6$Likelihood, decreasing = T),]
sweeps_sign_EU=top_EU[1:keep_eu,]

#Max values in top 5% for Ssa06
max(sweeps_sign_NA$Likelihood)
max(sweeps_sign_EU$Likelihood)

#Plot results wiht outliers as asterisks
sweep_6_plot=ggplot(sweeps6)+
  geom_rect(aes(xmin=7e+07, xmax= 7.3e+07, ymin=0, ymax=6), fill="gray80")+
  geom_point(aes(x=Position, y=Likelihood), col="dodgerblue3",)+
  geom_point(data=sweeps_NA6, aes(x=Position, y=Likelihood), col="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 6 position")+ylab("Likelihood")+
  geom_text(data=sweeps_sign_NA,
          nudge_y = -0.1,
          aes(x=Position, y=Likelihood),
          label = "*",col="red", size=12)+
  geom_text(data=sweeps_sign_EU, nudge_y = -0.1,
          aes(x=Position, y=Likelihood),  label = "*",col="dodgerblue3", size=12)
NULL


################################

#Read in values for Ssa13 -- the same as above was done for this chromosome
sweeps13=read.table("SweeD_Report.chr13_sweed_europe", header=T, skip=2, fill=T)
sweeps_NA13=read.table("SweeD_Report.chr13_sweed_northam", header=T, skip=2, fill=T)

#Top 5% of SNPs
keep_eu_13=nrow(sweeps13)*0.05
keep_na_13=nrow(sweeps_NA13)*0.05

#top North Am. snps
top_na_13=sweeps_NA13[order(sweeps_NA13$Likelihood, decreasing = T),]
sweeps_sign_NA_13=top_na_13[1:keep_na_13,]

#top European snps
top_EU_13=sweeps13[order(sweeps13$Likelihood, decreasing = T),]
sweeps_sign_EU_13=top_EU_13[1:keep_eu_13,]

#Max value of likelihoods for plotting
max(sweeps_sign_EU_13$Likelihood)
max(sweeps_sign_NA_13$Likelihood)

sweep_13_plot=ggplot(sweeps13)+  
  geom_rect(aes(xmin=7.5e+07, xmax= 7.8e+07, ymin=0, ymax=8), fill="gray80")+
  geom_point(aes(x=Position, y=Likelihood), col="dodgerblue3")+
  geom_point(data=sweeps_NA13, aes(x=Position, y=Likelihood), col="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 13 position")+ylab("Likelihood")+
  geom_text(data=sweeps_sign_NA_13, 
            nudge_y = -0.2,
            aes(x=Position, y=Likelihood),
            label = "*",col="red", size=12)+
  geom_text(data=sweeps_sign_EU_13, nudge_y = -0.2,
            aes(x=Position, y=Likelihood),  label = "*",col="dodgerblue3", size=12)
NULL


##########Ssa 16
#Read in values for Ssa16 -- the same as above was done for this chromosome
sweeps16=read.table("SweeD_Report.chr16_sweed_europe", header=T, skip=2, fill=T)
sweeps_NA16=read.table("SweeD_Report.chr16_sweed_northam", header=T, skip=2, fill=T)

#Top 5 percent of SNPs
keep_eu_16=nrow(sweeps16)*0.05
keep_na_16=nrow(sweeps_NA16)*0.05

top_na_16=sweeps_NA16[order(sweeps_NA16$Likelihood, decreasing = T),]
sweeps_sign_NA_16=top_na_16[1:keep_na_16,]

top_EU_16=sweeps16[order(sweeps16$Likelihood, decreasing = T),]
sweeps_sign_EU_16=top_EU_16[1:keep_eu_16,]

#Get max values for likelihood (plot)
max(sweeps_sign_EU_16$Likelihood)
max(sweeps_sign_NA_16$Likelihood)

#Plot
sweep_16_plot=ggplot(sweeps16)+  
  geom_rect(aes(xmin=3.6e+07, xmax= 3.7e+07, ymin=0, ymax=7), fill="gray80")+
  geom_point(aes(x=Position, y=Likelihood), col="dodgerblue3")+
  geom_point(data=sweeps_NA16, aes(x=Position, y=Likelihood), col="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 16 position")+ylab("Likelihood")+
    geom_text(data=sweeps_sign_NA_16,
         nudge_y = -0.125,
            aes(x=Position, y=Likelihood),
           label = "*",col="red", size=12)+
  geom_text(data=sweeps_sign_EU_16, nudge_y = -0.125,
            aes(x=Position, y=Likelihood),  label = "*",col="dodgerblue3", size=12)+
  #geom_smooth(data=sweeps_NA16, aes(x=Position, y=Likelihood), col="red", span=0.1, se =F)+
  #geom_smooth(data=sweeps16, aes(x=Position, y=Likelihood), col="dodgerblue3", span=0.1, se = F)
NULL

###########

#ssa19
#Read in Ssa19 results - same script as above
sweeps19=read.table("SweeD_Report.chr19_sweed_europe", header=T, skip=2, fill=T)
sweeps_NA19=read.table("SweeD_Report.chr19_sweed_northam", header=T, skip=2, fill=T)

#Get top 5% of SNPs
keep_eu_19=nrow(sweeps19)*0.05
keep_na_19=nrow(sweeps_NA19)*0.05

top_na_19=sweeps_NA19[order(sweeps_NA19$Likelihood, decreasing = T),]
sweeps_sign_NA_19=top_na_19[1:keep_na_19,]

top_EU_19=sweeps19[order(sweeps19$Likelihood, decreasing = T),]
sweeps_sign_EU_19=top_EU_19[1:keep_eu_19,]

#Check max values for likelihoods
max(sweeps_sign_NA_19$Likelihood)
max(sweeps_sign_EU_19$Likelihood)

#plotting
sweep_19_plot=ggplot(sweeps19)+  
geom_rect(aes(xmin=5.4e+07, xmax= 5.5e+07, ymin=0, ymax=6), fill="gray80")+
    geom_point(aes(x=Position, y=Likelihood), col="dodgerblue3")+
     geom_point(data=sweeps_NA19, aes(x=Position, y=Likelihood), col="red")+
  theme_bw()+theme(panel.grid = element_blank())+
  xlab("Chromosome 19  position")+ylab("Likelihood")+ylim(0,6)+
      geom_text(data=sweeps_sign_NA_19,
        nudge_y = -0.125,
           aes(x=Position, y=Likelihood),
      label = "*",col="red", size=12)+
  geom_text(data=sweeps_sign_EU_19, nudge_y = -0.125,
            aes(x=Position, y=Likelihood),  label = "*",col="dodgerblue3", size=12)+
  #geom_smooth(data=sweeps_NA16, aes(x=Position, y=Likelihood), col="red", span=0.1, se =F)+
  #geom_smooth(data=sweeps16, aes(x=Position, y=Likelihood), col="dodgerblue3", span=0.1, se = F)
  NULL

######### Put all plots together
#save as 14 x 10 - if appropriate
plot_grid(sweep_6_plot, sweep_13_plot, sweep_16_plot,sweep_19_plot, nrow=4)

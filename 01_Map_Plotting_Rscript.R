#Set directory
setwd("~/Desktop/Sarah/Salmon/TransAtlanticPaper/map/")

#Example of map for Canadian populations
#The same script was used for Norway with different lat/long range

## Load Libraries --------------
library(marmap) # for bathymetry if needed
library(ggplot2) 

#read coordinates for all sites - to add points for locations
coords=read.csv("All_80_Locations.csv", header=T)
#make sure all coordinates are within the limits of the map range (otherwise will add points anyway)

#Now get a bathy map for Canada area - Approximately Maine to Northern Labrador/QC
getNOAA.bathy(-49,-70,42,60, res=1,keep=T) -> bathydata

#Can save as PDF with this - optional (use with dev.off() below if saving)
#pdf(file="Map_autoplot_Canada.pdf", width = 13, height=12, bg = "white")

#Plot map
map=autoplot(bathydata, geom=c("r", "c"), colour="grey", size=0.1) +
  scale_fill_etopo() +
  geom_contour(aes(z=z), breaks=c(-100, -200, -500, -1000, -2000, -4000), colour="grey90", size=0.01) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(legend.position = "none")

#make sure there are no coordinates (populations) outside of map region (otherwise it will add those too without map area )
#Add population points
map + geom_point(data=coords[which(coords$Location=="CAN"),],
                 aes(x=Long, y=Lat), col="red", size=3, inherit.aes = F)

#dev.off() #finish saving to PDF if using

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

###### Hakalau Metacommunity 
###### Coding by Chris Wall, Cameron Egan
###### Nicole Hynson Lab, Pacific Biosciences Research Center, UH Mānoa
###### October 2019

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

# FYI all code running from R project, no "setwd()" needed; reference folders of interest in project wd


###### ###### ###### ###### 
########## packages
###### ###### ###### ###### 

if (!require("pacman")) install.packages("pacman"); library(pacman) # for rapid install if not in library
if (!require("BiocManager")) install.packages("BiocManager") # for buoinformatic packages
if (!require("devtools")) install.packages("devtools") # for developement tools

devtools::install_github('oswaldosantos/ggsn')

pacman::p_load("ade4", "multtest","phyloseq","rhdf5","ggplot2","colorspace","stringi", "geosphere", "ggplot2", "ggmap", "dplyr", "gridExtra", "geosphere", "sf", "raster", "spData", "tmap", "leaflet", "mapview", "shiny", "fossil", "RgoogleMaps", "devtools", "ggsn", "vegan", "multcomp", "dplyr", "grid", "scales", "gridExtra", "emmeans", "multcompView", "ggpubr", "Rmisc", "purrr", "RVAideMemoire", "RColorBrewer", "vegan")


##########################################################################
### Genearate MAPS of Hakalau showing habitat types and sampling plots ###
##########################################################################


##########################
####   Hawaii Map     ####
##########################

# Load API key (confidential)
API.key<-read.csv("data/API_Egan.csv")
API.key<-API.key[1,1]

register_google(key=API.key)
hi_map<-get_map(location=c(-157,20.5),zoom=7,maptype="satellite",color="color")
hi_map_for_man <- ggmap(hi_map) +
  geom_point(aes(x = -155.320, y = 19.83), pch=23,colour="black",fill="red", size = 2, stroke=0.5) +
  xlab("Longitude") + ylab("Latitude") +
  scale_y_continuous(limits=c(18.8, 22.2))+
  scale_x_continuous(limits=c(-160.2, -154)) +
  theme(axis.text=element_text(colour="black",size=8)) +
  theme(axis.title=element_text(colour="black",size=8)) +
  ggsn::scalebar(x.min=-160.2, x.max=-154, y.min=18.8, y.max=22.2, dist=50, dist_unit="km", transform=TRUE,
                 st.bottom=FALSE, st.size=2, box.fill=c("gray50", "white"), model="WGS84",st.color="white", border.size=0.5)

plot(hi_map_for_man)
ggsave("figures/hi_map.pdf", width= 5,height=3, plot=hi_map_for_man)


##########################
####   Forest Plots   ####
##########################
hakalau_map_zoom <-get_map(location=c(-155.320,19.83), zoom=14, maptype="satellite", color="color")

haka_metadata <-read.csv("data/haka_soil_metadata 2.csv", header=TRUE, row.names=1)

sampling_plots <- ggmap(hakalau_map_zoom) + 
  geom_point(data=haka_metadata,aes(x=lon,y=lat,fill=HabitatType),pch=21, stroke=0.3,colour="white",size=3) +
  scale_fill_manual(values=c("#88A550","#336B87")) +
  theme(legend.text=element_text(size=6),legend.title = element_text(size=8), legend.key=element_blank()) +
  xlab("Longitude") + ylab("Latitude") +
  scale_y_continuous(limits=c(19.81, 19.84)) +
  scale_x_continuous(limits=c(-155.328, -155.295)) +
  theme(axis.text=element_text(colour="black",size=8)) +
  theme(axis.title=element_text(colour="black",size=8)) +
  ggsn::scalebar(x.min=-155.345, x.max=-155.295, y.min=19.81, y.max=19.84, dist=0.5, dist_unit="km",
                 transform=TRUE, st.bottom=FALSE, st.size=2, box.fill=c("gray50", "white"), 
                 model="WGS84",st.color="white", border.size=0.5) 

plot(sampling_plots)
ggsave("figures/Sampling_plots.pdf", width= 6, height=5, plot=sampling_plots)
  


         ##############################################################
         #### Calculate lat and lon of each sample using geosphere ####
         ##############################################################

## GPS coordinates were only taken for focal trees (denoted by bearing = 0, and distance = 0 in metadata). 
## Need to calculate lat and lon for all samples to do spatial/metacommunity analysis
p <- as.matrix(cbind(haka_metadata$lon,haka_metadata$lat))
new_lat_lon <- as.data.frame(destPoint(p,haka_metadata$Bearing,haka_metadata$Distance))
colnames(new_lat_lon) <- c("Longitude","Latitude")
rownames(new_lat_lon) <- rownames(haka_metadata)
haka_meta <- cbind(haka_metadata,new_lat_lon)

write.csv(haka_meta, "output/haka_meta_dist.csv")

## PLOT MEAN GPS lat-long
lats<-aggregate(lat~HabitatType+Plot, data=haka_meta, FUN=mean)
longs<-aggregate(lon~HabitatType+Plot, data=haka_meta, FUN=mean)
meta_dist<-cbind(lats, longs[3])
write.csv(meta_dist, "output/meta_dist.csv")


## df to calculate geographic distances between SAMPLES
library(fossil)
haka_gps <- cbind(new_lat_lon)
geog_dists_km<-earth.dist(haka_gps[,1:2],dist=TRUE)
geog_dists_m <- (geog_dists_km * 1000)
geog_dists_m <- as.matrix(geog_dists_m)
rownames(geog_dists_m) <- rownames(haka_meta)
colnames(geog_dists_m)<- t(rownames(haka_meta))

write.csv(geog_dists_m,file="output/haka_dists.csv")

haka_dists <- lower.tri(geog_dists_m)


### df to calculate geographic distances between PLOTS
plot.dist<-as.data.frame(haka_gps)
plot.dist$sampleID<-as.factor(rownames(plot.dist))
plot.dist$HabitatType <- substr(plot.dist$sampleID, 0, 2) # extract first 2 letters of ID
plot.dist$Plot <- substr(plot.dist$sampleID, 0, 3)
plot.dist$Host <- substr(plot.dist$sampleID, 4, 5)
plot.dist<-na.omit(plot.dist)
plot.dist$HabitatType<-revalue(plot.dist$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest"))

# make a mean lat-long for each plot (pooling all samples in a plot)
plot.m.lat<-aggregate(Latitude~HabitatType+Plot, data=plot.dist, FUN=mean)
plot.m.long<-aggregate(Longitude~HabitatType+Plot, data=plot.dist, FUN=mean)
plot.m.dist<-cbind(plot.m.lat,plot.m.long[3])
write.csv(plot.m.dist,file="output/plot.m.dist.csv")

# Distance to neighboring plot in m
g_dists_km<-earth.dist(plot.m.dist[,3:4],dist=TRUE)
g_dists_m <- as.matrix(g_dists_km*1000)
rownames(g_dists_m) <- c("AK1", "AK2", "AK3", "AK4", "AK5", "AK6", "RO1", "RO2", "RO3", "RO4", "RO5", "RO6")
colnames(g_dists_m)<- t(plot.m.dist[2])

write.csv(g_dists_m,file="output/haka_plot_dists.csv")

#####################
### SEQUENCE DATA ###
#####################

   ###################################
   ####Import files and pre-process###
   ###################################

#Import files generated in QIIME
otu <- as.matrix(read.csv("data/haka_soil_ESV_table.csv", header = TRUE,row.names = 1))
haka_otu <- t(otu)
taxmat <- as.matrix(read.csv("data/haka_soil_taxonomy.csv", header = TRUE,row.names = 1))
sample <- haka_meta
all.equal(rownames(haka_otu), rownames(sample))

#Subset sample file to contain only the same rows as the otu file. Samples have been filtered out in
#QIIME for various reasons (too few reads, no Glom hits etc.)

#First create a vector of otu file rownames and call it "keep"
rows_to_keep<-rownames(haka_otu)

#Then subset sample file using the keep vector
sample<-sample[rows_to_keep,]
all.equal(rownames(haka_otu), rownames(sample))

#Order host names how you want them to appear in figures
sample$Host<- factor(sample$Host,levels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                                          "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                     labels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                              "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Write subsetted sample data as a dataframe
sampledata = sample_data(data.frame(Transect=sample$Transect, HabitatType=sample$HabitatType,
                                    Plot = sample$Plot, lat = sample$lat, lon = sample$lon, 
                                    TreeID= sample$TreeID, Host= sample$Host, 
                                    Bearing= sample$Bearing, Distance= sample$Distance,
                                    xCoords= sample$xCoords, yCoords= sample$yCoords,
                                    SampleType=sample$SampleType,
                                    OM= sample$OM, N= sample$N, P= sample$P.ppm, K= sample$K.ppm,
                                    Mg= sample$Mg.ppm, Ca= sample$Ca.ppm, Na= sample$Na.ppm, S= sample$S.ppm,
                                    pH= sample$pH,
                                    H= sample$H, CEC= sample$CationExchangeCapacity,
                                    K.cation= sample$K.cation, Mg.cation= sample$Mg.cation,
                                    Ca.cation= sample$Ca.cation, H.cation= sample$H.cation,
                                    Na.cation= sample$Na.cation, Longitude= sample$Longitude, 
                                    Latitude= sample$Latitude, stringsAsFactors = FALSE))

row.names(sampledata) <- row.names(sample)


#Change each file to the phyloseq format
OTU = otu_table(haka_otu, taxa_are_rows = FALSE)
physeq = phyloseq(OTU)
TAX = tax_table(taxmat)
haka_soil_physeq = merge_phyloseq(OTU, sampledata, TAX)

#Create new physeq object working at only the species level taxonomy. Collapses all ESVs identified as the same species
haka_VT_soil_physeq <- tax_glom(haka_soil_physeq,"Species")
ESV_rel_abund <- transform_sample_counts(haka_VT_soil_physeq,function(x)x/sum(x))

#Melt phyloseq object to make a dataframe for ggplot and bipartite
ESV_dataframe<-psmelt(ESV_rel_abund)


      #####################
      #### SPECIES RICH ###
      #####################

######################################################### 
###Phyloseq###
#Plots and analyses
par(mgp = c(3, 3, 0))

## By habitat type and plot
spec_rich_hab = plot_richness(haka_VT_soil_physeq, x ="Plot", measures="Observed", color="HabitatType") +
  geom_boxplot(col="black", aes(fill=HabitatType), alpha=0.8 , lwd=0.5) + 
  scale_color_manual(name= "Habitat Type", values=c("#88A550", "#336B87")) +
  scale_fill_manual(name= "Habitat Type", values=c("#88A550", "#336B87")) +
  scale_x_discrete(limits=c("RO1", "RO2", "RO3", "RO4", "RO5", "RO6",
                            "AK1", "AK2", "AK3", "AK4", "AK5", "AK6")) +
  scale_y_continuous(breaks=seq(0,45,by=10),limits=c(0,45)) +
  theme(text=element_text(colour="black",size=10)) + 
  ylab("AM fungal richness") + xlab("Habitat Type Plots") +
  theme(axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5,colour="black",size=8)) +
  theme(axis.text.y=element_text(colour="black",size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank(), legend.key=element_blank()) +
  facet_wrap(~SampleType) +
  theme(strip.text.x=element_blank(),strip.background = element_rect(fill="white")) 

plot(spec_rich_hab)
ggsave("figures/species_richness_by_plot.pdf", plot = spec_rich_hab, width=6,height=5)


#########################################################
# export richness from soil data and test for differences
soil_rich<-estimate_richness(haka_VT_soil_physeq, measures="Observed")
pairwise.wilcox.test(soil_rich$Observed, sample_data(haka_VT_soil_physeq)$Host) # not different
pairwise.wilcox.test(soil_rich$Observed, sample_data(haka_VT_soil_physeq)$Habitat) # not different 

# inspect alpha diversity a bit more
soil_rich$sampleID<-as.factor(rownames(soil_rich))
soil_rich$HabitatType <- substr(soil_rich$sampleID, 0, 2) # extract first 2 letters of ID
soil_rich$Plot <- substr(soil_rich$sampleID, 3, 3)
soil_rich$Host <- substr(soil_rich$sampleID, 4, 5)
soil_rich<-na.omit(soil_rich)
soil_rich$HabitatType<-revalue(soil_rich$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest"))


# GLM for soil richness
soil_rich.glm <- glm(Observed ~ HabitatType, data=soil_rich, family = poisson(link="log")) 
par(mfrow = c(2, 2))
plot(soil_rich.glm)
chi.anova<-anova(soil_rich.glm,test="Chisq")
chi.anova

# how does alpha diverstiy differ among Habitat Types? (18 AK vs. 18 RO)
soil_rich %>% 
  group_by(HabitatType) %>%
  dplyr::summarise(mean = mean(Observed))


#########################################################
#########################################################

         ######################
         ### Beta Diversity ###
         ######################

### Phyloseq ###
# Bray plot
bc_dist = as.matrix((vegdist(haka_otu, "bray")))
set.seed(520)
NMDS = metaMDS(bc_dist, k=3, trymax=100)
NMDS1=NMDS$points[,1]
NMDS2=NMDS$points[,2]
NMDS.plot.df=data.frame(NMDS1=NMDS1,NMDS2=NMDS2, Host=sample$Host, 
                     HabitatType=sample$HabitatType,Plot=sample$Plot)

NMDS.col<- c("#88A550","#336B87") # colors
Habitats<- NMDS.plot.df$HabitatType # habitat levels
groups.hab<-c("#88A550","#336B87")

par(mfrow = c(1, 1))
### make points semitransparent
# green
col2rgb("#88A550"); grn.60<-rgb(136, 165, 80, max=255, alpha=170)
# blue
col2rgb("#336B87"); blu.60<-rgb(51, 107, 135, max=255, alpha=170)
Hab.col<-c(grn.60, blu.60)



##### ##### ##### 
##### PLOTS ##### 
##### ##### ##### 


# Soil chemistry data
environmental_data <- read.csv("data/haka_chemistry.csv", header=TRUE, row.names=1)
all.equal(rownames(haka_otu), rownames(environmental_data))
rows_to_keep<-rownames(haka_otu)
environmental_data<-environmental_data[rows_to_keep,]
all.equal(rownames(haka_otu), rownames(environmental_data))

colnames(environmental_data)<-c("OM(%)","Total N", "P", "K", "Mg", "Ca", "Na", "S","pH", "H(meq/100g)", 
                                "CEC(meq/100g)", "K+", "Mg+2", "Ca+2", "H+", "Na+")


####################
####################
### summary of soil chemistry to be used in mean +/- SD table and plot level analyses
soil.chem<-environmental_data
colnames(soil.chem)<-c("OM.per","Total.N", "P", "K", "Mg", "Ca", "Na", "S","pH", "H.meq", 
                       "CEC.meq", "K.cat", "Mg.cat", "Ca.cat", "H.cat", "Na.cat")

# add factor levels
soil.chem$sampleID<-as.factor(rownames(soil.chem))
soil.chem$HabitatType <- as.factor(substr(soil.chem$sampleID, 0, 2)) # extract first 2 letters of ID
soil.chem$Plot <- as.factor(substr(soil.chem$sampleID, 0, 3))
soil.chem$Host <- as.factor(substr(soil.chem$sampleID, 4, 5))

# modify factor levels
soil.chem$HabitatType<-as.factor(revalue(soil.chem$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest")))

# remove NAs
# don't use 'soil.chem<-na.omit(soil.chem)' because there is one sample with NAs for some not all metrics
soil.chem<-soil.chem[!is.na(soil.chem$P),]

# remove dupicate rows where no additional analyses were run
# these duplicate rows represent "sampling intensity" under the same plant (DNA), but not soil chem reps 
soil.chem.dup.rem<-soil.chem %>% 
  distinct(OM.per, Total.N, P, K, Mg, Ca, Na, S, pH, .keep_all = TRUE)

soil.chem.dup.rem<-soil.chem.dup.rem[, c(18:20, 17, 1:16)]

# write csv
write.csv(soil.chem.dup.rem, "output/soil.chem.dup.rem.csv")

#################################
############### Sample size
# sample size for Plot (~14-25)
## table SD
chem.n.plot <- soil.chem.dup.rem %>% 
  group_by(Plot) %>%
  summarize(n =n())

# sample size for Habitats (99 RO, 128 AK)
## table SD
chem.n.habitat <- soil.chem.dup.rem %>% 
  group_by(HabitatType) %>%
  summarize(n =n())


#################################
###### soil chemistry box plots

# diagnostic plots for soil chemistry
for(i in c(5:20)){
  Y<-soil.chem.dup.rem[,i]
  full<-lm(Y~HabitatType, data=soil.chem.dup.rem, na.action=na.exclude)
  R <- resid(full) #save glm residuals
  op<-par(mfrow = c(2,2), mar=c(5,4,1,2), pty="sq")
  plot(full, add.smooth = FALSE, which=1)
  QQ <- qqnorm(R, main = colnames(soil.chem.dup.rem)[i]) 
  QQline <- qqline(R)
  hist(R, xlab="Residuals", main = colnames(soil.chem.dup.rem)[i])
  plot(soil.chem.dup.rem$HabitatType, R, xlab=colnames(soil.chem.dup.rem)[i], ylab="Residuals")
}


# make plots

chem<-soil.chem.dup.rem[,c(1:2,5:20)]
plot_list = list()
for(i in c(3:18)){
  p <-ggplot(chem, aes_string(y=chem[,i], x=chem[,2], group=chem[,2])) +
    ylab(colnames(chem)[i]) +
    xlab("Plot") +
    theme_bw() +
    theme(axis.text.x=element_text(size = 4)) +
    geom_boxplot(col="black", aes(fill=HabitatType), alpha=0.8 , lwd=0.5) + 
    scale_color_manual(name= "Habitat Type", values=c("#88A550", "#336B87")) +
    scale_fill_manual(name= "Habitat Type", values=c("#88A550", "#336B87")) +
    scale_x_discrete(limits=c("RO1", "RO2", "RO3", "RO4", "RO5", "RO6",
                              "AK1", "AK2", "AK3", "AK4", "AK5", "AK6"))
  plot_list[[i]] = p
  print(p)
}

# arrange list plots
soil.chem.plot1<-cowplot::plot_grid(plotlist = plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], plot_list[[11]])

soil.chem.plot2<-cowplot::plot_grid(plotlist= plot_list[[2]], plot_list[[12]], plot_list[[13]], plot_list[[14]], plot_list[[15]], plot_list[[16]], plot_list[[17]], plot_list[[18]])
  
    
pdf(file="output/soil.chem.plots.pdf", height=8, width=12)
soil.chem.plot1
soil.chem.plot2
dev.off()        

##########################
#####
#   means by replicate plots
library(dplyr)
sum_plot.df <- soil.chem.dup.rem %>% 
  group_by(Plot) %>%
  summarize(
    OM = mean(OM.per, na.rm=TRUE),
    Total.N = mean(Total.N, na.rm=TRUE),
    P = mean(P),
    K = mean(K),
    Mg = mean(Mg),
    Ca = mean(Ca),
    Na = mean(Na),
    S = mean(S),
    pH = mean(pH),
    H.meq = mean(H.meq),
    CEC.meq = mean(CEC.meq),
    K.cat = mean(K.cat),
    Mg.cat = mean(Mg.cat),
    Ca.cat = mean(Ca.cat),
    H.cat = mean(H.cat),
    Na.cat = mean(Na.cat)
  )

plot.soil.chem<-as.data.frame(sum_plot.df)
plot.soil.chem$HabitatType <- substr(plot.soil.chem$Plot, 0, 2)
plot.soil.chem$HabitatType<-as.factor(revalue(plot.soil.chem$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest")))

#####
# means by forest
forest.chem.df <- soil.chem.dup.rem %>% 
  group_by(HabitatType) %>%
  summarize(
    OM = mean(OM.per, na.rm=TRUE),
    Total.N = mean(Total.N, na.rm=TRUE),
    P = mean(P),
    K = mean(K),
    Mg = mean(Mg),
    Ca = mean(Ca),
    Na = mean(Na),
    S = mean(S),
    pH = mean(pH),
    H.meq = mean(H.meq),
    CEC.meq = mean(CEC.meq),
    K.cat = mean(K.cat),
    Mg.cat = mean(Mg.cat),
    Ca.cat = mean(Ca.cat),
    H.cat = mean(H.cat),
    Na.cat = mean(Na.cat)
  )

forest.chem.df<-as.data.frame(forest.chem.df)

forest.chem.mean<-as.data.frame(t(forest.chem.df)) # transpose
colnames(forest.chem.mean)<-c("AK.mean", "RO.mean")


#######
## table mean, n, SE for PLOT level


sum_plot.se <- soil.chem.dup.rem[,c(-1,-3,-4)] %>% 
  group_by(Plot) %>%
  dplyr::summarise_each(funs(mean(., na.rm=T), n = sum(!is.na(.)), se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))))

se = sd(.)/sqrt(n())

sum_plot.se <- soil.chem.dup.rem %>% 
  group_by(Plot) %>%
  plyr::summarize(
    OM = se(OM.per),
    Total.N = se(Total.N),
    P = se(P),
    K = se(K),
    Mg = se(Mg),
    Ca = se(Ca),
    Na = se(Na),
    S = se(S),
    pH = se(pH),
    H.meq = se(H.meq),
    CEC.meq = se(CEC.meq),
    K.cat = se(K.cat),
    Mg.cat = se(Mg.cat),
    Ca.cat = se(Ca.cat),
    H.cat = se(H.cat),
    Na.cat = se(Na.cat)
  )

soil.plot.chem.se<-as.data.frame(sum_plot.se)
soil.plot.chem.se$HabitatType <- substr(soil.plot.chem.se$Plot, 0, 2)
soil.plot.chem.se$HabitatType<-as.factor(revalue(soil.plot.chem.se$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest")))

######
# now se for habitat
## table SE
forest.chem.se <- soil.chem.dup.rem %>% 
  group_by(HabitatType) %>%
  summarize(
    OM = se(OM.per),
    Total.N = se(Total.N),
    P = se(P),
    K = se(K),
    Mg = se(Mg),
    Ca = se(Ca),
    Na = se(Na),
    S = se(S),
    pH = se(pH),
    H.meq = se(H.meq),
    CEC.meq = se(CEC.meq),
    K.cat = se(K.cat),
    Mg.cat = se(Mg.cat),
    Ca.cat = se(Ca.cat),
    H.cat = se(H.cat),
    Na.cat = se(Na.cat)
  )

forest.chem.se<-as.data.frame(forest.chem.se)
forest.chem.se.df<-as.data.frame(t(forest.chem.se)) # transpose
colnames(forest.chem.se.df)<-c("RO.SE", "AK.SE")

soil.chem.sum<- cbind(forest.chem.mean, forest.chem.se.df)
soil.chem.sum<-soil.chem.sum[, c(2,4,1,3)]# re-order

#export
write.csv(soil.chem.sum, "output/soil.chem.sum.csv")


################################
### t tests for forest soil chem
################################
# dataframe = 'plot.soil.chem' = averages at plot level

# OM -- p-value = 0.0003921
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$OM,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$OM,
       paired=FALSE, var.equal=FALSE)

# Total.N -- p-value = 0.0003909
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$Total.N,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$Total.N,
       paired=FALSE, var.equal=FALSE)

# P -- p-value = 0.5392
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$P,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$P,
       paired=FALSE, var.equal=FALSE)

# K -- p-value = 0.03977
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$K,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$K,
       paired=FALSE, var.equal=FALSE)

# Mg -- p-value = 0.7262
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$Mg,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$Mg,
       paired=FALSE, var.equal=FALSE)

# Ca -- p-value = 0.6473
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$Ca,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$Ca,
       paired=FALSE, var.equal=FALSE)

# Na -- p-value = 0.4954
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$Na,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$Na,
       paired=FALSE, var.equal=FALSE)

# S -- p-value = 1.986e-06
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$S,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$S,
       paired=FALSE, var.equal=FALSE)

# pH -- p-value = 0.0001515
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$pH,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$pH,
       paired=FALSE, var.equal=FALSE)

# H.meq -- p-value = 0.02044
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$H.meq,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$H.meq,
       paired=FALSE, var.equal=FALSE)

# CEC.meq -- p-value = 0.09223
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$CEC.meq,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$CEC.meq,
       paired=FALSE, var.equal=FALSE)

# K.cat -- p-value = 0.03011
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$K.cat,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$K.cat,
       paired=FALSE, var.equal=FALSE)

# Mg.cat -- p-value = 0.002985
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$Mg.cat,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$Mg.cat,
       paired=FALSE, var.equal=FALSE)

# Ca.cat -- p-value = 0.003329
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$Ca.cat,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$Ca.cat,
       paired=FALSE, var.equal=FALSE)

# H.cat -- p-value = 0.0001391
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$H.cat,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$H.cat,
       paired=FALSE, var.equal=FALSE)

# Na.cat -- p-value = 0.03775
t.test(subset(plot.soil.chem, HabitatType == "Remnant Forest")$Na.cat,
       subset(plot.soil.chem, HabitatType == "Restored Forest")$Na.cat,
       paired=FALSE, var.equal=FALSE)




##### ##### ##### ##### 
##### NMDS with vectors
##### ##### ##### ##### 

par(mfrow=c(1,1))

#Fit environmental vectors to ordination to see which environmental variables are correlated with the ordination
# ions wihtout "+" are in ppm
fit.env <- envfit(NMDS, environmental_data, na.rm=TRUE)


### make plot by habitat with vectors for environment
environm_plot<-ordiplot(NMDS, type="n", main=substitute(paste("")), cex.main=1, display="sites", xlim=c(-0.25, 0.8), cex.lab=0.8, cex.axis=0.8)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(NMDS, "sites", cex=0.8, pch=16, col=Hab.col[Habitats])
ordiellipse(NMDS, groups=Habitats, kind="sd", draw="polygon", conf=0.95, alpha=30, col=groups.hab, border=groups.hab)
legend("topright", legend=levels(Habitats), cex=1, pch=16, col=groups.hab, pt.cex=1, bty="n")
par.new=T
plot(fit.env, col="black", p.max=0.05, cex=0.9, lwd=1)

dev.copy(pdf, "figures/environm.NMDS.pdf", height=5.5, width=6)
dev.off() 


### Stats ###
# PERMANOVA for NMDS plot
# Bray PERMANOVA 
haka.bc.adonis <- adonis(bc_dist~HabitatType, data=sample, permutations = 9999)
haka.bc.adonis

betadisper(bc_dist, sample$HabitatType)


######## ######## ######## ######## ######## Make environmental data as PC1 and PC2
######## explore environmental data, use PCA to visualize
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
require(graphics)
library(plyr)


############### PCA

df.PCA<-environmental_data
df.PCA$sampleID<-as.factor(rownames(df.PCA))
df.PCA$HabitatType <- substr(df.PCA$sampleID, 0, 2) # extract first 2 letters of ID
df.PCA$Plot <- substr(df.PCA$sampleID, 3, 3)
df.PCA$Host <- substr(df.PCA$sampleID, 4, 5)
df.PCA<-na.omit(df.PCA)
df.PCA$HabitatType<-revalue(df.PCA$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest"))

# remove columns unnecessary for final analysis, few factors retained
env.PCA<-df.PCA[ , !names(df.PCA) %in% c("sampleID", "Plot", "Host", "HabitatType")]
Hak.env.PCA <- prcomp(env.PCA, center = TRUE, scale= TRUE) # with HabitatType in dataframe
PC.summary<-(summary(Hak.env.PCA))
ev<-Hak.env.PCA$sdev^2
newdat<-Hak.env.PCA$x[,1:4]
plot(Hak.env.PCA, type="lines", main="Hak.env.PCA eigenvalues")


############### PCA and NMDS combine
NMDS.otu<-as.data.frame(NMDS1) # NMDS Bray Curtis of OTUs
PCA.env<-as.data.frame(newdat) # the exported data from PCA of environment
PC.NMD<-merge(PCA.env, NMDS.otu, by = "row.names", all = TRUE) # merge dataframes
PC.NMD<-na.omit(PC.NMD) #drop NA columns that don't correspond
colnames(PC.NMD)[1]<-"sampleID"
PC.NMD$HabitatType <- as.factor(substr(PC.NMD$sampleID, 0, 2)) # extract first 2 letters of ID
PC.NMD$HabitatType<-revalue(PC.NMD$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest"))
PC.NMD$HabitatType<- factor(PC.NMD$HabitatType, levels=c("Remnant Forest", "Restored Forest"))
PC.NMD$Plot <- as.factor(substr(PC.NMD$sampleID, 3, 3)) # make plot ID
PC.NMD$Host <- as.factor(substr(PC.NMD$sampleID, 4, 5)) # make host ID


##### ##### ##### ##### ##### ##### 
## Tests for relationship pf PC x NMDS

# test RO/Remnant Forest relationship
RO.PCNMD<-PC.NMD[(PC.NMD$HabitatType=="Remnant Forest"),]
mod.RO<-lm(NMDS1~PC1, data=RO.PCNMD); print(anova(mod.RO), digits=6)

# test AK/Restored Forest relationship
AK.PCNMD<-PC.NMD[(PC.NMD$HabitatType=="Restored Forest"),]
mod.AK<-lm(NMDS1~PC1, data=AK.PCNMD); print(anova(mod.AK), digits=6)

# full model, used in plotting
mod<-lm(NMDS1~PC1*HabitatType, data=PC.NMD)
print(anova(mod), digits=6)

int<-coef(summary(mod))[1] # intercept
PC1.coef<-coef(summary(mod))[2] # PC1
HabAK.coef<-coef(summary(mod))[3] # HabitatType Restored Forest
PC1.AK.coef<-coef(summary(mod))[4] # PC1:HabitatTypeRestored Forest

library(plotrix)
plot(NMDS1~PC1, data=PC.NMD, pch=16,
     col=Hab.col[as.factor(HabitatType)],
     xlab="PC1 (soil chemistry)", 
     ylab="NMDS1 (AMF Bray-Curtis)")
ablineclip(int, (PC1.coef) , col=NMDS.col[1], lwd=2, 
           x1 = min(PC.NMD$PC1[PC.NMD$HabitatType=="Remnant Forest"], na.rm=T), 
           x2 = max(PC.NMD$PC1[PC.NMD$HabitatType=="Remnant Forest"], na.rm=T)) # RO model
ablineclip((int+HabAK.coef), (PC1.coef + PC1.AK.coef), col=NMDS.col[2], lwd=2, 
           x1 = min(PC.NMD$PC1[PC.NMD$HabitatType=="Restored Forest"], na.rm=T), 
           x2 = max(PC.NMD$PC1[PC.NMD$HabitatType=="Restored Forest"], na.rm=T)) # AK model
legend("topleft", c("Remnant Forest", "Restored Forest"), lty=c(1,1), lwd=c(2,2), col=NMDS.col, cex=0.8, pch=16, y.intersp = 0.5, bty="n")

#dev.copy(pdf, "figures/PC.NMD.slopes.pdf", height=5, width=6)
#dev.off()

################
##### using ggplot
BW.back<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.2))

PC.NMDS.plot<-ggplot(data=PC.NMD, aes(x=PC1, y=NMDS1, color=HabitatType))+
  geom_point(aes(color=HabitatType), alpha=.7) + 
  scale_color_manual(values=Hab.col) + 
  scale_fill_manual(values=Hab.col) +
  xlab("PC1 (soil chemistry)") +
  ylab("NMDS1 (AMF Bray-Curtis)") +
  geom_smooth(aes(color=HabitatType, fill=HabitatType), values=Hab.col, method = "lm", alpha = .2, size=0.6) +
  theme(legend.title = element_blank(),
        axis.text=element_text(size=6),
        axis.title=element_text(size=10)) +
  theme_bw() + BW.back

PC.NMDS.plot
dev.copy(pdf, "figures/PC.NMD.slopes2.pdf", height=3, width=6)
dev.off()




########################
### Spatial Analyses ###
########################

 #####################################################
 ## Mantel tests on distance-dissimilarity patterns ##
 ####################################################

#To run a Mantel test, we will need to generate two distance matrices: one containing spatial 
#distances and one containing Bray-Curtis distances between samples at the given points.  In the 
#spatial distance matrix, entries for pairs of points that are close together are lower than for 
#pairs of points that are far apart.  In the BC matrix, entries for pairs of locations 
#with similar outcomes are lower than for pairs of points with dissimilar outcomes.  We do this 
#using the dist function.  The Mantel test function will require objects of this "distance" class.

######################
# entire dataset
haka_meta # metadata
hak_otu<-t(otu) # transpose

# merge to have all data together
Study.dat<-merge(haka_meta, hak_otu, by = "row.names", all = TRUE)
Study.dat<-Study.dat[!is.na(Study.dat[35]),]

#separate dataframes
haka_met_mant<-Study.dat[,c(1:34)] # the meta data
hak_otu<-Study.dat[,c(35:1468)] # OTU data

# geographic distance
dists_km<-earth.dist(haka_met_mant[,31:32], dist=TRUE)
dists_m <- (dists_km * 1000)
dists_m <- as.matrix(dists_m)
rownames(dists_m) <- rownames(haka_met_mant)
colnames(dists_m)<- t(rownames(haka_met_mant))

geomat<-dists_m # rename
dists_m<-as.dist(dists_m) # make distance matrix

study_bc_dist = as.dist((vegdist(hak_otu, "bray"))); study_bc_mat<-as.matrix(vegdist(hak_otu, "bray"))

# mantel
study_mantel_test<-mantel.rtest(dists_m, study_bc_dist, nrepet=9999)
study_mantel_test # significant

# new dataframe
Study.df<-data.frame(Distance=geomat[lower.tri(geomat)],
                  BrayCurtis= study_bc_mat[lower.tri(study_bc_mat)])

##############
#Generate plots for all samples combined
Study_dist_plot<- ggplot(Study.df, aes(x=log(Distance+1), y=BrayCurtis)) +
  geom_point(size=0.8,alpha=0.5) +
  geom_smooth(method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')), lwd=0.7)+ 
  ggtitle("Hakalau")+ 
  scale_y_continuous(name="Bray-Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,2),limits=c(0,8)) +
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.x=element_text(colour="black",size=7), 
        axis.title.x=element_blank(),
        axis.text.y=element_text(colour="black",size=7),
        axis.title=element_text(colour="black",size=9),
        legend.title=element_text(colour="black",size=7,face="bold"),
        legend.text=element_text(colour="black",size=7),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

Study_dist_plot
### use this later to combine with Habitat-Specific Mantels ###

mantel.mod<-lm(BrayCurtis~log(Distance+1), data=Study.df) # linear model fit, 
mantel.mod[1] # coefficients


#####################
#####################  Separated by sites
#####################

# Subset ESV table and metadata into RO and AK habitat types
RO_meta_data <- subset(haka_meta, HabitatType == "Remnant Forest")

##################### RO

## subset the OTU file to give only RO
ROdat<-t(otu[,grep("^RO", colnames(otu))]) # transposed data with only RO
samples.RO<-sample[grep("^RO", rownames(sample)),] # only rows with "RO"
all.equal(rownames(ROdat), rownames(samples.RO)) # rows match

RO<- merge(RO_meta_data, ROdat, by = "row.names", all = TRUE) # all data merged so that row names match
#drop NAs from samples that have no sequence data
RO<-RO[!is.na(RO[35]),]

RO_meta<-RO[,c(1:34)] # the meta data
RO_otu<-RO[,c(35:1468)] # OTU data

#Calculate geographic distances among sample points on subset data
RO_dists_km<-earth.dist(RO_meta[,31:32], dist=TRUE)
RO_dists_m <- (RO_dists_km * 1000)
RO_dists_m <- as.matrix(RO_dists_m)
rownames(RO_dists_m) <- rownames(RO_meta)
colnames(RO_dists_m)<- t(rownames(RO_meta))

RO_geomat<-RO_dists_m
RO_dists_m<-as.dist(RO_dists_m) # make distance matrix


######################## AK 

## subset the OTU file to give only AK
AK_meta_data <-subset(haka_meta, HabitatType== "Restored Forest")

AKdat<-t(otu[,grep("^AK", colnames(otu))]) # transposed data with only AK
samples.AK<-sample[grep("^AK", rownames(sample)),] # only rows with "AK"
all.equal(rownames(AKdat), rownames(samples.AK)) # rows match

AK<- merge(AK_meta_data, AKdat, by = "row.names", all = TRUE) # all data merged so that row names match
#drop NAs from samples that have no sequence data
AK<-AK[!is.na(AK[35]),]

AK_meta<-AK[,c(1:34)] # the meta data
AK_otu<-AK[,c(35:1468)] # OTU data

AK_dists_km<-earth.dist(AK_meta[,31:32],dist=TRUE)
AK_dists_m <- (AK_dists_km * 1000)
AK_dists_m <- as.matrix(AK_dists_m)
rownames(AK_dists_m) <- rownames(AK_meta)
colnames(AK_dists_m)<- t(rownames(AK_meta))

AK_geomat<-AK_dists_m
AK_dists_m<-as.dist(AK_dists_m) # make distance matrix
######################

#Calculate Bray-Curtis distances on subset data
RO_bc_dist = as.dist((vegdist(RO_otu, "bray"))); RO_bc_mat<-as.matrix(vegdist(RO_otu, "bray"))
AK_bc_dist = as.dist((vegdist(AK_otu, "bray"))); AK_bc_mat<-as.matrix(vegdist(AK_otu, "bray"))


#Use distance matrices to test for a correlation. The test consists of calculating the 
#correlation of the entries in the matrices, then permuting the matrices and calculating the 
#same test statistic under each permutation and comparing the original test statistic to the 
#distribution of test statistics from the permutations to generate a p-value. The number of 
#permutations defines the precision with which the p-value can be calculated.  The function 
#to perform the Mantel test is mantel.rtest and the required arguments are the two distance 
#matrices. The number of permutations used are 9999.

###########
RO_mantel_test<-mantel.rtest(RO_dists_m, RO_bc_dist, nrepet=9999)
RO_mantel_test
RO_df<-data.frame(Distance=RO_geomat[lower.tri(RO_geomat)],
                    BrayCurtis= RO_bc_mat[lower.tri(RO_bc_mat)])
                    #Host=RO_meta$Host)
summary(RO_df)


###########
AK_mantel_test <- mantel.rtest(AK_dists_m, AK_bc_dist,nrepet=9999)
AK_mantel_test
AK_df<-data.frame(Distance=AK_geomat[lower.tri(AK_geomat)],
                    BrayCurtis=AK_bc_mat[lower.tri(AK_bc_mat)])
                    #Host=AK_meta$Host)
summary(AK_df)

############

#Export dataframe, remove NAs, bring dataframe back into R
write.csv(RO_df,file="output/RO_dist_dataframe.csv")
write.csv(AK_df,file="output/AK_dist_dataframe.csv")

RO_dist_df<-read.csv("output/RO_dist_dataframe.csv", header=T, as.is=T)
AK_dist_df<-read.csv("output/AK_dist_dataframe.csv",header=T, as.is=T)


##############
#Generate plots for all elevations combined

##### AK plot
AK_dist_plot<- ggplot(AK_df,aes(x=log(Distance+1),y=BrayCurtis)) +
  geom_point(size=1,alpha=0.5) +
  geom_smooth(method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')), lwd=0.7) +
  ggtitle("Restored Forest")+ 
  scale_y_continuous(name="Bray-Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,2),limits=c(0,8)) +
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.text.x=element_text(colour="black",size=7),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())


mantel.mod.AK<-lm(BrayCurtis~log(Distance+1), data=AK_df) # linear model fit, 
mantel.mod.AK[1] # coefficients

###############
##### AK plot

# not significant, so no line for RO: trend +
RO_dist_plot<- ggplot(RO_df,aes(x=log(Distance+1), y=BrayCurtis)) +
  geom_point(size=1,alpha=0.5) +
  ggtitle("Remnant Forest")+
  scale_y_continuous(name="Bray-Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,2),limits=c(0,8)) +
  theme(plot.title = element_text(size=10, hjust = 0.5),
        axis.title.x=element_blank(), 
        axis.text.x=element_text(colour="black",size=7),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())


#####################################################################
# combine 3 mantel plots
dist.plots<-plot_grid(Study_dist_plot, RO_dist_plot, AK_dist_plot,
              nrow = 1, rel_widths = c(1, 0.85, 0.85), labels = c("(A)", "(B)", "(C)"), label_size = 10) 
ggdraw(add_sub(dist.plots, "           log(Distance (m)+1)", hjust = 0.5, y=0.8, size=10))

dev.copy(png, "figures/All.Hak_RO.AK.distdecay.png",  width= 2000, height=1500, res=300)
dev.off()


#####################################################################
#####################################################################
## AK-RO proximity Mantels
#####################################################################
#####################################################################

# entire dataset
haka_meta # metadata
hak_otu<-t(otu) # transpose


#separate dataframes
haka_met_mant<-Study.dat[,c(1:32)] # the meta data
hak_otu<-Study.dat[,c(33:1468)] # OTU data

### distance from plots (coded)
plot.dist # all dist
plot.m.dist # mean dist
plot.m.dist$Transect<-c(8,8,8,9,9,9,8,8,8,9,9,9) #add Transect

# set lat-long to mean value for full dataframe


haka_meta$Lat.mean<-ifelse(haka_meta$Plot=="AK1", plot.m.dist[1,3],
                    ifelse(haka_meta$Plot=="AK2", plot.m.dist[2,3],
                    ifelse(haka_meta$Plot=="AK3", plot.m.dist[3,3],
                    ifelse(haka_meta$Plot=="AK4", plot.m.dist[4,3],
                    ifelse(haka_meta$Plot=="AK5", plot.m.dist[5,3],
                    ifelse(haka_meta$Plot=="AK6", plot.m.dist[6,3],
                    ifelse(haka_meta$Plot=="RO1", plot.m.dist[7,3],
                    ifelse(haka_meta$Plot=="RO2", plot.m.dist[8,3],
                    ifelse(haka_meta$Plot=="RO3", plot.m.dist[9,3],
                    ifelse(haka_meta$Plot=="RO4", plot.m.dist[10,3],
                    ifelse(haka_meta$Plot=="RO5", plot.m.dist[11,3],
                              plot.m.dist[12,3])  # last for RO6
                               ))))))))))

haka_meta$Long.mean<-ifelse(haka_meta$Plot=="AK1", plot.m.dist[1,4],
                     ifelse(haka_meta$Plot=="AK2", plot.m.dist[2,4],
                     ifelse(haka_meta$Plot=="AK3", plot.m.dist[3,4],
                     ifelse(haka_meta$Plot=="AK4", plot.m.dist[4,4],
                     ifelse(haka_meta$Plot=="AK5", plot.m.dist[5,4],
                     ifelse(haka_meta$Plot=="AK6", plot.m.dist[6,4],
                     ifelse(haka_meta$Plot=="RO1", plot.m.dist[7,4],
                     ifelse(haka_meta$Plot=="RO2", plot.m.dist[8,4],
                     ifelse(haka_meta$Plot=="RO3", plot.m.dist[9,4],
                     ifelse(haka_meta$Plot=="RO4", plot.m.dist[10,4],
                     ifelse(haka_meta$Plot=="RO5", plot.m.dist[11,4],
                               plot.m.dist[12,4]) # last for RO6
                               ))))))))))


########## separate distance matrix transect
##### Transect 8

haka_met_mant.T8<-haka_meta[(haka_meta$Transect==8),]

########### separate OTUs accordingly
## subset the OTU file to give only T8 data = AK 1-3 and RO 1-3
# transposed data with only Transect 8
T8dat<-t(otu[, grepl("^AK1", colnames(otu)) | 
               grepl("^AK2", colnames(otu)) |
               grepl("^AK3", colnames(otu)) |
               grepl("^RO1", colnames(otu)) | 
               grepl("^RO2", colnames(otu)) |
               grepl("^RO3", colnames(otu))])

# only rows from T8
samples.T8<-sample[grepl("^AK1", rownames(sample)) |
                     grepl("^AK2", rownames(sample)) |
                     grepl("^AK3", rownames(sample)) |
                     grepl("^RO1", rownames(sample)) |
                     grepl("^RO2", rownames(sample)) |
                     grepl("^RO3", rownames(sample)),]

all.equal(rownames(T8dat), rownames(samples.T8)) # rows match
T8.df<- merge(haka_met_mant.T8, T8dat, by = "row.names", all = TRUE) # all data merged so that row names match
T8.df<-na.omit(T8.df) # remove NAs

T8_meta<-T8.df[,c(1:34)] # the meta data
T8_otu<-T8.df[,c(35:1470)] # OTU data

# Distance to neighboring plots
g_dists_km.T8<- earth.dist(T8_meta[,33:34],dist=TRUE)
g_dists_m.T8<- (g_dists_km.T8 * 1000)
g_dists_m.T8 <- as.matrix(g_dists_m.T8)
rownames(g_dists_m.T8) <- rownames(T8_meta)
colnames(g_dists_m.T8)<- t(rownames(T8_meta))

# rename and make distance matrix
geomat.T8<-g_dists_m.T8
g_dists_m.T8<-as.dist(g_dists_m.T8) 


######## Transect 9
haka_met_mant.T9<-haka_meta[(haka_meta$Transect==9),]

##### subset the OTU file to give only T9 data = AK 4-6 and RO 4-6
# transposed data with only Transect 9
T9dat<-t(otu[, grepl("^AK4", colnames(otu)) | 
               grepl("^AK5", colnames(otu)) |
               grepl("^AK6", colnames(otu)) |
               grepl("^RO4", colnames(otu)) | 
               grepl("^RO5", colnames(otu)) |
               grepl("^RO6", colnames(otu))])

# only rows from T9
samples.T9<-sample[grepl("^AK4", rownames(sample)) |
                     grepl("^AK5", rownames(sample)) |
                     grepl("^AK6", rownames(sample)) |
                     grepl("^RO4", rownames(sample)) |
                     grepl("^RO5", rownames(sample)) |
                     grepl("^RO6", rownames(sample)),]

all.equal(rownames(T9dat), rownames(samples.T9)) # rows match
T9.df<- merge(haka_met_mant.T9, T9dat, by = "row.names", all = TRUE, na.rm=T) # all data merged so that row names match
T9.df<-na.omit(T9.df) # remove NAs
T9_meta<-T9.df[,c(1:34)] # the meta data
T9_otu<-T9.df[,c(35:1470)] # OTU data


# Distance to neighboring plots
g_dists_km.T9<- earth.dist(T9_meta[,33:34],dist=TRUE)
g_dists_m.T9<- (g_dists_km.T9 * 1000)
g_dists_m.T9 <- as.matrix(g_dists_m.T9)
rownames(g_dists_m.T9) <- rownames(T9_meta)
colnames(g_dists_m.T9)<- t(rownames(T9_meta))

# rename and make distance matrix
geomat.T9<-g_dists_m.T9
g_dists_m.T9<-as.dist(g_dists_m.T9) 


############################
############################
#Calculate Bray-Curtis distances on T8 and T9 data
T8_bc_dist = as.dist((vegdist(T8_otu, "bray"))); T8_bc_mat<-as.matrix(vegdist(T8_otu, "bray"))
T9_bc_dist = as.dist((vegdist(T9_otu, "bray"))); T9_bc_mat<-as.matrix(vegdist(T9_otu, "bray"))


# Mantel Test
T8_mantel_test<-mantel.rtest(g_dists_m.T8, T8_bc_dist, nrepet=9999)
T8_mantel_test

T9_mantel_test<-mantel.rtest(g_dists_m.T9, T9_bc_dist, nrepet=9999)
T9_mantel_test

# make df
T8_bcdf<-data.frame(Distance=geomat.T8[lower.tri(geomat.T8)],
                    BrayCurtis= T8_bc_mat[lower.tri(T8_bc_mat)])
summary(T8_bcdf)

# make df
T9_bcdf<-data.frame(Distance=geomat.T9[lower.tri(geomat.T9)],
                    BrayCurtis= T9_bc_mat[lower.tri(T9_bc_mat)])
summary(T9_bcdf)

#### mantel modesl
mantel.mod.T8<-lm(BrayCurtis~log(Distance+1), data=T8_bcdf) # linear model fit
mantel.mod.T9<-lm(BrayCurtis~log(Distance+1), data=T9_bcdf) # linear model fit


### plots
T8_dist_plot<- ggplot(T8_bcdf,aes(x=log(Distance+1), y=BrayCurtis)) +
  geom_point(size=1,alpha=0.5) +
  theme(text=element_text(colour="black",size=12)) + 
  ggtitle("transect 8") +
  geom_smooth(method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')), lwd=0.7) +
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=7), 
        axis.text.y=element_text(colour="black",size=7),
        axis.title=element_text(colour="black",size=9),
        legend.title=element_text(colour="black",size=7,face="bold"),
        legend.text=element_text(colour="black",size=7)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

plot(T8_dist_plot)
ggsave("figures/T8_dist_plot.pdf", 
       plot = T8_dist_plot, width = 5, height = 5)


#### T9 plot
T9_dist_plot<- ggplot(T9_bcdf,aes(x=log(Distance+1), y=BrayCurtis)) +
  geom_point(size=1,alpha=0.5) +
  theme(text=element_text(colour="black",size=12)) + 
  geom_smooth(method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')), lwd=0.7) +
  ggtitle("transect 9") +
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=7), 
        axis.text.y=element_text(colour="black",size=7),
        axis.title=element_text(colour="black",size=9),
        legend.title=element_text(colour="black",size=7,face="bold"),
        legend.text=element_text(colour="black",size=7)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

plot(T9_dist_plot)
ggsave("figures/T9_dist_plot.pdf", 
       plot = T9_dist_plot, width = 5, height = 5)


##################
####################################
######################################################
# what about plot-plot comparisons within a transect. 
# how does each RO compare to each AK within a transect
# Use AK close, AK further, AK furthest 
# compare to close/further/furthest RO
# generate Bray-Curtis dissimiliarity
######################################################
####################################
##################

# Transect 9 : North 
# plot distance: W --> E (AK4, AK5, AK6, RO4, RO5, RO6)

# Transect 8 : South 
# plot distance: W --> E (AK1, AK2, AK3, RO1, RO2, RO3)


# T8 dat = transposed OTUs
# sample.T8 = samples
# T8.df = dataframe 
# T8_meta = the meta data
# T8_otu = otus

T8dat<-t(otu[, grepl("^AK1", colnames(otu)) | 
               grepl("^AK2", colnames(otu)) |
               grepl("^AK3", colnames(otu)) |
               grepl("^RO1", colnames(otu)) | 
               grepl("^RO2", colnames(otu)) |
               grepl("^RO3", colnames(otu))])

# only rows from T8
samples.T8<-sample[grepl("^AK1", rownames(sample)) |
                     grepl("^AK2", rownames(sample)) |
                     grepl("^AK3", rownames(sample)) |
                     grepl("^RO1", rownames(sample)) |
                     grepl("^RO2", rownames(sample)) |
                     grepl("^RO3", rownames(sample)),]

all.equal(rownames(T8dat), rownames(samples.T8)) # rows match
T8.df<- merge(haka_met_mant.T8, T8dat, by = "row.names", all = TRUE) # all data merged so that row names match
T8.df<-na.omit(T8.df) # remove NAs

T8_meta<-T8.df[,c(1:34)] # the meta data
T8_otu<-T8.df[,c(35:1470)] # OTU data

## T8 single plot dataframes
T8.AK1.df<-T8.df[(T8.df$Plot=="AK1"),]; T8.AK1.meta<- T8.AK1.df[,c(1:34)]; T8.AK1.otu<-T8.AK1.df[,c(35:1470)]
T8.AK2.df<-T8.df[(T8.df$Plot=="AK2"),]; T8.AK2.meta<- T8.AK2.df[,c(1:34)]; T8.AK2.otu<-T8.AK2.df[,c(35:1470)]
T8.AK3.df<-T8.df[(T8.df$Plot=="AK3"),]; T8.AK3.meta<- T8.AK3.df[,c(1:34)]; T8.AK3.otu<-T8.AK3.df[,c(35:1470)]
T8.RO1.df<-T8.df[(T8.df$Plot=="RO1"),]; T8.RO1.meta<- T8.RO1.df[,c(1:34)]; T8.RO1.otu<-T8.RO1.df[,c(35:1470)]
T8.RO2.df<-T8.df[(T8.df$Plot=="RO2"),]; T8.RO2.meta<- T8.RO2.df[,c(1:34)]; T8.RO2.otu<-T8.RO2.df[,c(35:1470)]
T8.RO3.df<-T8.df[(T8.df$Plot=="RO3"),]; T8.RO3.meta<- T8.RO3.df[,c(1:34)]; T8.RO3.otu<-T8.RO3.df[,c(35:1470)]

# combine for comparisons
T8.RO1.AK1.df<-rbind(T8.RO1.df, T8.AK1.df)
T8.RO1.AK1.meta <-rbind(T8.RO1.meta,T8.AK1.meta)
T8.RO1.AK1.otu<-rbind(T8.RO1.otu,T8.AK1.otu)

T8.RO1.AK2.df<-rbind(T8.RO1.df, T8.AK2.df)
T8.RO1.AK2.meta <-rbind(T8.RO1.meta,T8.AK2.meta)
T8.RO1.AK2.otu<-rbind(T8.RO1.otu,T8.AK2.otu)

T8.RO1.AK3.df<-rbind(T8.RO1.df, T8.AK3.df)
T8.RO1.AK3.meta <-rbind(T8.RO1.meta,T8.AK3.meta)
T8.RO1.AK3.otu<-rbind(T8.RO1.otu,T8.AK3.otu)

T8.RO2.AK1.df<-rbind(T8.RO2.df, T8.AK1.df)
T8.RO2.AK1.meta <-rbind(T8.RO2.meta,T8.AK1.meta)
T8.RO2.AK1.otu<-rbind(T8.RO2.otu,T8.AK1.otu)

T8.RO2.AK2.df<-rbind(T8.RO2.df, T8.AK2.df)
T8.RO2.AK2.meta <-rbind(T8.RO2.meta,T8.AK2.meta)
T8.RO2.AK2.otu<-rbind(T8.RO2.otu,T8.AK2.otu)

T8.RO2.AK3.df<-rbind(T8.RO2.df, T8.AK3.df)
T8.RO2.AK3.meta <-rbind(T8.RO2.meta,T8.AK3.meta)
T8.RO2.AK3.otu<-rbind(T8.RO2.otu,T8.AK3.otu)

T8.RO3.AK1.df<-rbind(T8.RO3.df, T8.AK1.df)
T8.RO3.AK1.meta <-rbind(T8.RO3.meta,T8.AK1.meta)
T8.RO3.AK1.otu<-rbind(T8.RO3.otu,T8.AK1.otu)

T8.RO3.AK2.df<-rbind(T8.RO3.df, T8.AK2.df)
T8.RO3.AK2.meta <-rbind(T8.RO3.meta,T8.AK2.meta)
T8.RO3.AK2.otu<-rbind(T8.RO3.otu,T8.AK2.otu)

T8.RO3.AK3.df<-rbind(T8.RO3.df, T8.AK3.df)
T8.RO3.AK3.meta <-rbind(T8.RO3.meta,T8.AK3.meta)
T8.RO3.AK3.otu<-rbind(T8.RO3.otu,T8.AK3.otu)


# Distance to neighboring plots
######### RO1
# RO1 - AK1
g_dists_km.T8.RO1.AK1<- earth.dist(T8.RO1.AK1.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO1.AK1<- (g_dists_km.T8.RO1.AK1 * 1000)
g_dists_m.T8.RO1.AK1 <- as.matrix(g_dists_m.T8.RO1.AK1)
rownames(g_dists_m.T8.RO1.AK1) <- rownames(T8.RO1.AK1.meta)
colnames(g_dists_m.T8.RO1.AK1)<- t(rownames(T8.RO1.AK1.meta))
# rename and make distance matrix
geomat.T8.RO1.AK1<-g_dists_m.T8.RO1.AK1
g_dists_m.T8.RO1.AK1<-as.dist(g_dists_m.T8.RO1.AK1) 


#### RO1 - AK2
g_dists_km.T8.RO1.AK2<- earth.dist(T8.RO1.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO1.AK2<- (g_dists_km.T8.RO1.AK2 * 1000)
g_dists_m.T8.RO1.AK2 <- as.matrix(g_dists_m.T8.RO1.AK2)
rownames(g_dists_m.T8.RO1.AK2) <- rownames(T8.RO1.AK2.meta)
colnames(g_dists_m.T8.RO1.AK2)<- t(rownames(T8.RO1.AK2.meta))
# rename and make distance matrix
geomat.T8.RO1.AK2<-g_dists_m.T8.RO1.AK2
g_dists_m.T8.RO1.AK2<-as.dist(g_dists_m.T8.RO1.AK2) 


#### RO1 - AK3
g_dists_km.T8.RO1.AK3<- earth.dist(T8.RO1.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO1.AK3<- (g_dists_km.T8.RO1.AK3 * 1000)
g_dists_m.T8.RO1.AK3 <- as.matrix(g_dists_m.T8.RO1.AK3)
rownames(g_dists_m.T8.RO1.AK3) <- rownames(T8.RO1.AK3.meta)
colnames(g_dists_m.T8.RO1.AK3)<- t(rownames(T8.RO1.AK3.meta))
# rename and make distance matrix
geomat.T8.RO1.AK3<-g_dists_m.T8.RO1.AK3
g_dists_m.T8.RO1.AK3<-as.dist(g_dists_m.T8.RO1.AK3) 

################ RO2
# RO2 - AK1
g_dists_km.T8.RO2.AK1<- earth.dist(T8.RO2.AK1.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO2.AK1<- (g_dists_km.T8.RO2.AK1 * 1000)
g_dists_m.T8.RO2.AK1 <- as.matrix(g_dists_m.T8.RO2.AK1)
rownames(g_dists_m.T8.RO2.AK1) <- rownames(T8.RO2.AK1.meta)
colnames(g_dists_m.T8.RO2.AK1)<- t(rownames(T8.RO2.AK1.meta))
# rename and make distance matrix
geomat.T8.RO2.AK1<-g_dists_m.T8.RO2.AK1
g_dists_m.T8.RO2.AK1<-as.dist(g_dists_m.T8.RO2.AK1) 


#### RO2 - AK2
g_dists_km.T8.RO2.AK2<- earth.dist(T8.RO2.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO2.AK2<- (g_dists_km.T8.RO2.AK2 * 1000)
g_dists_m.T8.RO2.AK2 <- as.matrix(g_dists_m.T8.RO2.AK2)
rownames(g_dists_m.T8.RO2.AK2) <- rownames(T8.RO2.AK2.meta)
colnames(g_dists_m.T8.RO2.AK2)<- t(rownames(T8.RO2.AK2.meta))
# rename and make distance matrix
geomat.T8.RO2.AK2<-g_dists_m.T8.RO2.AK2
g_dists_m.T8.RO2.AK2<-as.dist(g_dists_m.T8.RO2.AK2) 


#### RO2 - AK3
g_dists_km.T8.RO2.AK3<- earth.dist(T8.RO2.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO2.AK3<- (g_dists_km.T8.RO2.AK3 * 1000)
g_dists_m.T8.RO2.AK3 <- as.matrix(g_dists_m.T8.RO2.AK3)
rownames(g_dists_m.T8.RO2.AK3) <- rownames(T8.RO2.AK3.meta)
colnames(g_dists_m.T8.RO2.AK3)<- t(rownames(T8.RO2.AK3.meta))
# rename and make distance matrix
geomat.T8.RO2.AK3<-g_dists_m.T8.RO2.AK3
g_dists_m.T8.RO2.AK3<-as.dist(g_dists_m.T8.RO2.AK3) 

################ RO3
# RO3 - AK1
g_dists_km.T8.RO3.AK1<- earth.dist(T8.RO3.AK1.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO3.AK1<- (g_dists_km.T8.RO3.AK1 * 1000)
g_dists_m.T8.RO3.AK1 <- as.matrix(g_dists_m.T8.RO3.AK1)
rownames(g_dists_m.T8.RO3.AK1) <- rownames(T8.RO3.AK1.meta)
colnames(g_dists_m.T8.RO3.AK1)<- t(rownames(T8.RO3.AK1.meta))
# rename and make distance matrix
geomat.T8.RO3.AK1<-g_dists_m.T8.RO3.AK1
g_dists_m.T8.RO3.AK1<-as.dist(g_dists_m.T8.RO3.AK1) 


#### RO3 - AK2
g_dists_km.T8.RO3.AK2<- earth.dist(T8.RO3.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO3.AK2<- (g_dists_km.T8.RO3.AK2 * 1000)
g_dists_m.T8.RO3.AK2 <- as.matrix(g_dists_m.T8.RO3.AK2)
rownames(g_dists_m.T8.RO3.AK2) <- rownames(T8.RO3.AK2.meta)
colnames(g_dists_m.T8.RO3.AK2)<- t(rownames(T8.RO3.AK2.meta))
# rename and make distance matrix
geomat.T8.RO3.AK2<-g_dists_m.T8.RO3.AK2
g_dists_m.T8.RO3.AK2<-as.dist(g_dists_m.T8.RO3.AK2) 


#### RO3 - AK3
g_dists_km.T8.RO3.AK3<- earth.dist(T8.RO3.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO3.AK3<- (g_dists_km.T8.RO3.AK3 * 1000)
g_dists_m.T8.RO3.AK3 <- as.matrix(g_dists_m.T8.RO3.AK3)
rownames(g_dists_m.T8.RO3.AK3) <- rownames(T8.RO3.AK3.meta)
colnames(g_dists_m.T8.RO3.AK3)<- t(rownames(T8.RO3.AK3.meta))
# rename and make distance matrix
geomat.T8.RO3.AK3<-g_dists_m.T8.RO3.AK3
g_dists_m.T8.RO3.AK3<-as.dist(g_dists_m.T8.RO3.AK3) 


############################
#Calculate Bray-Curtis distances on T8 between RO to AK

#### RO1
T8_RO1.AK1.bc_dist = as.dist((vegdist(T8.RO1.AK1.otu, "bray")))
T8_RO1.AK1.bc_mat<-as.matrix(vegdist(T8.RO1.AK1.otu, "bray"))

T8_RO1.AK2.bc_dist = as.dist((vegdist(T8.RO1.AK2.otu, "bray")))
T8_RO1.AK2.bc_mat<-as.matrix(vegdist(T8.RO1.AK2.otu, "bray"))

T8_RO1.AK3.bc_dist = as.dist((vegdist(T8.RO1.AK3.otu, "bray")))
T8_RO1.AK3.bc_mat<-as.matrix(vegdist(T8.RO1.AK3.otu, "bray"))

#### RO2
T8_RO2.AK1.bc_dist = as.dist((vegdist(T8.RO2.AK1.otu, "bray")))
T8_RO2.AK1.bc_mat<-as.matrix(vegdist(T8.RO2.AK1.otu, "bray"))

T8_RO2.AK2.bc_dist = as.dist((vegdist(T8.RO2.AK2.otu, "bray")))
T8_RO2.AK2.bc_mat<-as.matrix(vegdist(T8.RO2.AK2.otu, "bray"))

T8_RO2.AK3.bc_dist = as.dist((vegdist(T8.RO2.AK3.otu, "bray")))
T8_RO2.AK3.bc_mat<-as.matrix(vegdist(T8.RO2.AK3.otu, "bray"))

#### RO3
T8_RO3.AK1.bc_dist = as.dist((vegdist(T8.RO3.AK1.otu, "bray")))
T8_RO3.AK1.bc_mat<-as.matrix(vegdist(T8.RO3.AK1.otu, "bray"))

T8_RO3.AK2.bc_dist = as.dist((vegdist(T8.RO3.AK2.otu, "bray")))
T8_RO3.AK2.bc_mat<-as.matrix(vegdist(T8.RO3.AK2.otu, "bray"))

T8_RO3.AK3.bc_dist = as.dist((vegdist(T8.RO3.AK3.otu, "bray")))
T8_RO3.AK3.bc_mat<-as.matrix(vegdist(T8.RO3.AK3.otu, "bray"))

#################
# Mantel Test
T8.RO1.AK1_mantel_test<-mantel.rtest(g_dists_m.T8.RO1.AK1, T8_RO1.AK1.bc_dist, nrepet=9999)
T8.RO1.AK1_mantel_test

T8.RO1.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.RO1.AK2, T8_RO1.AK2.bc_dist, nrepet=9999)
T8.RO1.AK2_mantel_test

T8.RO1.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.RO1.AK3, T8_RO1.AK3.bc_dist, nrepet=9999)
T8.RO1.AK3_mantel_test


### RO2
T8.RO2.AK1_mantel_test<-mantel.rtest(g_dists_m.T8.RO2.AK1, T8_RO2.AK1.bc_dist, nrepet=9999)
T8.RO2.AK1_mantel_test

T8.RO2.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.RO2.AK2, T8_RO2.AK2.bc_dist, nrepet=9999)
T8.RO2.AK2_mantel_test

T8.RO2.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.RO2.AK3, T8_RO2.AK3.bc_dist, nrepet=9999)
T8.RO2.AK3_mantel_test

### RO3
T8.RO3.AK1_mantel_test<-mantel.rtest(g_dists_m.T8.RO3.AK1, T8_RO3.AK1.bc_dist, nrepet=9999)
T8.RO3.AK1_mantel_test

T8.RO3.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.RO3.AK2, T8_RO3.AK2.bc_dist, nrepet=9999)
T8.RO3.AK2_mantel_test

T8.RO3.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.RO3.AK3, T8_RO3.AK3.bc_dist, nrepet=9999)
T8.RO3.AK3_mantel_test

###############
# make df
### RO1
T8_RO1.AK1.bcdf<-data.frame(Distance=geomat.T8.RO1.AK1[lower.tri(geomat.T8.RO1.AK1)],
                            BrayCurtis= T8_RO1.AK1.bc_mat[lower.tri(T8_RO1.AK1.bc_mat)])

T8_RO1.AK2.bcdf<-data.frame(Distance=geomat.T8.RO1.AK2[lower.tri(geomat.T8.RO1.AK2)],
                            BrayCurtis= T8_RO1.AK2.bc_mat[lower.tri(T8_RO1.AK2.bc_mat)])

T8_RO1.AK3.bcdf<-data.frame(Distance=geomat.T8.RO1.AK3[lower.tri(geomat.T8.RO1.AK3)],
                            BrayCurtis= T8_RO1.AK3.bc_mat[lower.tri(T8_RO1.AK3.bc_mat)])

### RO2
T8_RO2.AK1.bcdf<-data.frame(Distance=geomat.T8.RO2.AK1[lower.tri(geomat.T8.RO2.AK1)],
                            BrayCurtis= T8_RO2.AK1.bc_mat[lower.tri(T8_RO2.AK1.bc_mat)])

T8_RO2.AK2.bcdf<-data.frame(Distance=geomat.T8.RO2.AK2[lower.tri(geomat.T8.RO2.AK2)],
                            BrayCurtis= T8_RO2.AK2.bc_mat[lower.tri(T8_RO2.AK2.bc_mat)])

T8_RO2.AK3.bcdf<-data.frame(Distance=geomat.T8.RO2.AK3[lower.tri(geomat.T8.RO2.AK3)],
                            BrayCurtis= T8_RO2.AK3.bc_mat[lower.tri(T8_RO2.AK3.bc_mat)])

### RO3
T8_RO3.AK1.bcdf<-data.frame(Distance=geomat.T8.RO3.AK1[lower.tri(geomat.T8.RO3.AK1)],
                            BrayCurtis= T8_RO3.AK1.bc_mat[lower.tri(T8_RO3.AK1.bc_mat)])

T8_RO3.AK2.bcdf<-data.frame(Distance=geomat.T8.RO3.AK2[lower.tri(geomat.T8.RO3.AK2)],
                            BrayCurtis= T8_RO3.AK2.bc_mat[lower.tri(T8_RO3.AK2.bc_mat)])

T8_RO3.AK3.bcdf<-data.frame(Distance=geomat.T8.RO3.AK3[lower.tri(geomat.T8.RO3.AK3)],
                            BrayCurtis= T8_RO3.AK3.bc_mat[lower.tri(T8_RO3.AK3.bc_mat)])

# mantel test df
prox.mantT8<-as.data.frame(rbind(T8.RO1.AK1_mantel_test[1], T8.RO1.AK2_mantel_test[1],T8.RO1.AK3_mantel_test[1],
                                 T8.RO2.AK1_mantel_test[1], T8.RO2.AK2_mantel_test[1], T8.RO2.AK3_mantel_test[1],
                                 T8.RO3.AK1_mantel_test[1], T8.RO3.AK2_mantel_test[1], T8.RO3.AK3_mantel_test[1]))
prox.mantT8$obs<-as.numeric(prox.mantT8$obs)

prox.mantT8$pvalue<-as.numeric(rbind(
  T8.RO1.AK1_mantel_test[5], T8.RO1.AK2_mantel_test[5], T8.RO1.AK3_mantel_test[5],
  T8.RO2.AK1_mantel_test[5], T8.RO2.AK2_mantel_test[5], T8.RO2.AK3_mantel_test[5],
  T8.RO3.AK1_mantel_test[5], T8.RO3.AK2_mantel_test[5], T8.RO3.AK3_mantel_test[5]))
prox.mantT8$comparison<-as.factor(c("T8.RO1.AK1", "T8.RO1.AK2", "T8.RO1.AK3",
                                    "T8.RO2.AK1", "T8.RO2.AK2", "T8.RO2.AK3",
                                    "T8.RO3.AK1", "T8.RO3.AK2", "T8.RO3.AK3"))
prox.mantT8$AK.group<-as.factor(c("AK far", "AK mid", "AK near", 
                                  "AK far", "AK mid", "AK near", 
                                  "AK far", "AK mid", "AK near"))
prox.mantT8$Transect<-"8"

################# 
################# 
################# Transect 9
# Transect 9 : North 
# plot distance: W --> E (AK4, AK5, AK6, RO4, RO5, RO6)


# T9 dat = transposed OTUs
# sample.T9 = samples
# T9.df = dataframe 
# T9_meta = the meta data
# T9_otu = otus


T9dat<-t(otu[, grepl("^AK4", colnames(otu)) | 
               grepl("^AK5", colnames(otu)) |
               grepl("^AK6", colnames(otu)) |
               grepl("^RO4", colnames(otu)) | 
               grepl("^RO5", colnames(otu)) |
               grepl("^RO6", colnames(otu))])

# only rows from T9
samples.T9<-sample[grepl("^AK4", rownames(sample)) |
                     grepl("^AK5", rownames(sample)) |
                     grepl("^AK6", rownames(sample)) |
                     grepl("^RO4", rownames(sample)) |
                     grepl("^RO5", rownames(sample)) |
                     grepl("^RO6", rownames(sample)),]

all.equal(rownames(T9dat), rownames(samples.T9)) # rows match
T9.df<- merge(haka_met_mant.T9, T9dat, by = "row.names", all = TRUE) # all data merged so that row names match
T9.df<-na.omit(T9.df) # remove NAs

T9_meta<-T9.df[,c(1:34)] # the meta data
T9_otu<-T9.df[,c(35:1470)] # OTU data

## T9 single plot dataframes
T9.AK4.df<-T9.df[(T9.df$Plot=="AK4"),]; T9.AK4.meta<- T9.AK4.df[,c(1:34)]; T9.AK4.otu<-T9.AK4.df[,c(35:1470)]
T9.AK5.df<-T9.df[(T9.df$Plot=="AK5"),]; T9.AK5.meta<- T9.AK5.df[,c(1:34)]; T9.AK5.otu<-T9.AK5.df[,c(35:1470)]
T9.AK6.df<-T9.df[(T9.df$Plot=="AK6"),]; T9.AK6.meta<- T9.AK6.df[,c(1:34)]; T9.AK6.otu<-T9.AK6.df[,c(35:1470)]
T9.RO4.df<-T9.df[(T9.df$Plot=="RO4"),]; T9.RO4.meta<- T9.RO4.df[,c(1:34)]; T9.RO4.otu<-T9.RO4.df[,c(35:1470)]
T9.RO5.df<-T9.df[(T9.df$Plot=="RO5"),]; T9.RO5.meta<- T9.RO5.df[,c(1:34)]; T9.RO5.otu<-T9.RO5.df[,c(35:1470)]
T9.RO6.df<-T9.df[(T9.df$Plot=="RO6"),]; T9.RO6.meta<- T9.RO6.df[,c(1:34)]; T9.RO6.otu<-T9.RO6.df[,c(35:1470)]

# combine for comparisons
T9.RO4.AK4.df<-rbind(T9.RO4.df, T9.AK4.df)
T9.RO4.AK4.meta <-rbind(T9.RO4.meta,T9.AK4.meta)
T9.RO4.AK4.otu<-rbind(T9.RO4.otu,T9.AK4.otu)

T9.RO4.AK5.df<-rbind(T9.RO4.df, T9.AK5.df)
T9.RO4.AK5.meta <-rbind(T9.RO4.meta,T9.AK5.meta)
T9.RO4.AK5.otu<-rbind(T9.RO4.otu,T9.AK5.otu)

T9.RO4.AK6.df<-rbind(T9.RO4.df, T9.AK6.df)
T9.RO4.AK6.meta <-rbind(T9.RO4.meta,T9.AK6.meta)
T9.RO4.AK6.otu<-rbind(T9.RO4.otu,T9.AK6.otu)

T9.RO5.AK4.df<-rbind(T9.RO5.df, T9.AK4.df)
T9.RO5.AK4.meta <-rbind(T9.RO5.meta,T9.AK4.meta)
T9.RO5.AK4.otu<-rbind(T9.RO5.otu,T9.AK4.otu)

T9.RO5.AK5.df<-rbind(T9.RO5.df, T9.AK5.df)
T9.RO5.AK5.meta <-rbind(T9.RO5.meta,T9.AK5.meta)
T9.RO5.AK5.otu<-rbind(T9.RO5.otu,T9.AK5.otu)

T9.RO5.AK6.df<-rbind(T9.RO5.df, T9.AK6.df)
T9.RO5.AK6.meta <-rbind(T9.RO5.meta,T9.AK6.meta)
T9.RO5.AK6.otu<-rbind(T9.RO5.otu,T9.AK6.otu)

T9.RO6.AK4.df<-rbind(T9.RO6.df, T9.AK4.df)
T9.RO6.AK4.meta <-rbind(T9.RO6.meta,T9.AK4.meta)
T9.RO6.AK4.otu<-rbind(T9.RO6.otu,T9.AK4.otu)

T9.RO6.AK5.df<-rbind(T9.RO6.df, T9.AK5.df)
T9.RO6.AK5.meta <-rbind(T9.RO6.meta,T9.AK5.meta)
T9.RO6.AK5.otu<-rbind(T9.RO6.otu,T9.AK5.otu)

T9.RO6.AK6.df<-rbind(T9.RO6.df, T9.AK6.df)
T9.RO6.AK6.meta <-rbind(T9.RO6.meta,T9.AK6.meta)
T9.RO6.AK6.otu<-rbind(T9.RO6.otu,T9.AK6.otu)


# Distance to neighboring plots
######### RO4
# RO4 - AK4
g_dists_km.T9.RO4.AK4<- earth.dist(T9.RO4.AK4.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO4.AK4<- (g_dists_km.T9.RO4.AK4 * 1000)
g_dists_m.T9.RO4.AK4 <- as.matrix(g_dists_m.T9.RO4.AK4)
rownames(g_dists_m.T9.RO4.AK4) <- rownames(T9.RO4.AK4.meta)
colnames(g_dists_m.T9.RO4.AK4)<- t(rownames(T9.RO4.AK4.meta))
# rename and make distance matrix
geomat.T9.RO4.AK4<-g_dists_m.T9.RO4.AK4
g_dists_m.T9.RO4.AK4<-as.dist(g_dists_m.T9.RO4.AK4) 


#### RO4 - AK5
g_dists_km.T9.RO4.AK5<- earth.dist(T9.RO4.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO4.AK5<- (g_dists_km.T9.RO4.AK5 * 1000)
g_dists_m.T9.RO4.AK5 <- as.matrix(g_dists_m.T9.RO4.AK5)
rownames(g_dists_m.T9.RO4.AK5) <- rownames(T9.RO4.AK5.meta)
colnames(g_dists_m.T9.RO4.AK5)<- t(rownames(T9.RO4.AK5.meta))
# rename and make distance matrix
geomat.T9.RO4.AK5<-g_dists_m.T9.RO4.AK5
g_dists_m.T9.RO4.AK5<-as.dist(g_dists_m.T9.RO4.AK5) 


#### RO4 - AK6
g_dists_km.T9.RO4.AK6<- earth.dist(T9.RO4.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO4.AK6<- (g_dists_km.T9.RO4.AK6 * 1000)
g_dists_m.T9.RO4.AK6 <- as.matrix(g_dists_m.T9.RO4.AK6)
rownames(g_dists_m.T9.RO4.AK6) <- rownames(T9.RO4.AK6.meta)
colnames(g_dists_m.T9.RO4.AK6)<- t(rownames(T9.RO4.AK6.meta))
# rename and make distance matrix
geomat.T9.RO4.AK6<-g_dists_m.T9.RO4.AK6
g_dists_m.T9.RO4.AK6<-as.dist(g_dists_m.T9.RO4.AK6) 

################ RO5
# RO5 - AK4
g_dists_km.T9.RO5.AK4<- earth.dist(T9.RO5.AK4.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO5.AK4<- (g_dists_km.T9.RO5.AK4 * 1000)
g_dists_m.T9.RO5.AK4 <- as.matrix(g_dists_m.T9.RO5.AK4)
rownames(g_dists_m.T9.RO5.AK4) <- rownames(T9.RO5.AK4.meta)
colnames(g_dists_m.T9.RO5.AK4)<- t(rownames(T9.RO5.AK4.meta))
# rename and make distance matrix
geomat.T9.RO5.AK4<-g_dists_m.T9.RO5.AK4
g_dists_m.T9.RO5.AK4<-as.dist(g_dists_m.T9.RO5.AK4) 


#### RO5 - AK5
g_dists_km.T9.RO5.AK5<- earth.dist(T9.RO5.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO5.AK5<- (g_dists_km.T9.RO5.AK5 * 1000)
g_dists_m.T9.RO5.AK5 <- as.matrix(g_dists_m.T9.RO5.AK5)
rownames(g_dists_m.T9.RO5.AK5) <- rownames(T9.RO5.AK5.meta)
colnames(g_dists_m.T9.RO5.AK5)<- t(rownames(T9.RO5.AK5.meta))
# rename and make distance matrix
geomat.T9.RO5.AK5<-g_dists_m.T9.RO5.AK5
g_dists_m.T9.RO5.AK5<-as.dist(g_dists_m.T9.RO5.AK5) 


#### RO5 - AK6
g_dists_km.T9.RO5.AK6<- earth.dist(T9.RO5.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO5.AK6<- (g_dists_km.T9.RO5.AK6 * 1000)
g_dists_m.T9.RO5.AK6 <- as.matrix(g_dists_m.T9.RO5.AK6)
rownames(g_dists_m.T9.RO5.AK6) <- rownames(T9.RO5.AK6.meta)
colnames(g_dists_m.T9.RO5.AK6)<- t(rownames(T9.RO5.AK6.meta))

# rename and make distance matrix
geomat.T9.RO5.AK6<-g_dists_m.T9.RO5.AK6
g_dists_m.T9.RO5.AK6<-as.dist(g_dists_m.T9.RO5.AK6) 

################ RO6
# RO6 - AK4
g_dists_km.T9.RO6.AK4<- earth.dist(T9.RO6.AK4.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO6.AK4<- (g_dists_km.T9.RO6.AK4 * 1000)
g_dists_m.T9.RO6.AK4 <- as.matrix(g_dists_m.T9.RO6.AK4)
rownames(g_dists_m.T9.RO6.AK4) <- rownames(T9.RO6.AK4.meta)
colnames(g_dists_m.T9.RO6.AK4)<- t(rownames(T9.RO6.AK4.meta))
# rename and make distance matrix
geomat.T9.RO6.AK4<-g_dists_m.T9.RO6.AK4
g_dists_m.T9.RO6.AK4<-as.dist(g_dists_m.T9.RO6.AK4) 


#### RO6 - AK5
g_dists_km.T9.RO6.AK5<- earth.dist(T9.RO6.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO6.AK5<- (g_dists_km.T9.RO6.AK5 * 1000)
g_dists_m.T9.RO6.AK5 <- as.matrix(g_dists_m.T9.RO6.AK5)
rownames(g_dists_m.T9.RO6.AK5) <- rownames(T9.RO6.AK5.meta)
colnames(g_dists_m.T9.RO6.AK5)<- t(rownames(T9.RO6.AK5.meta))
# rename and make distance matrix
geomat.T9.RO6.AK5<-g_dists_m.T9.RO6.AK5
g_dists_m.T9.RO6.AK5<-as.dist(g_dists_m.T9.RO6.AK5) 


#### RO6 - AK6
g_dists_km.T9.RO6.AK6<- earth.dist(T9.RO6.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO6.AK6<- (g_dists_km.T9.RO6.AK6 * 1000)
g_dists_m.T9.RO6.AK6 <- as.matrix(g_dists_m.T9.RO6.AK6)
rownames(g_dists_m.T9.RO6.AK6) <- rownames(T9.RO6.AK6.meta)
colnames(g_dists_m.T9.RO6.AK6)<- t(rownames(T9.RO6.AK6.meta))
# rename and make distance matrix
geomat.T9.RO6.AK6<-g_dists_m.T9.RO6.AK6
g_dists_m.T9.RO6.AK6<-as.dist(g_dists_m.T9.RO6.AK6) 


############################
#Calculate Bray-Curtis distances on T9 between RO to AK

#### RO4
T9_RO4.AK4.bc_dist = as.dist((vegdist(T9.RO4.AK4.otu, "bray")))
T9_RO4.AK4.bc_mat<-as.matrix(vegdist(T9.RO4.AK4.otu, "bray"))

T9_RO4.AK5.bc_dist = as.dist((vegdist(T9.RO4.AK5.otu, "bray")))
T9_RO4.AK5.bc_mat<-as.matrix(vegdist(T9.RO4.AK5.otu, "bray"))

T9_RO4.AK6.bc_dist = as.dist((vegdist(T9.RO4.AK6.otu, "bray")))
T9_RO4.AK6.bc_mat<-as.matrix(vegdist(T9.RO4.AK6.otu, "bray"))

#### RO5
T9_RO5.AK4.bc_dist = as.dist((vegdist(T9.RO5.AK4.otu, "bray")))
T9_RO5.AK4.bc_mat<-as.matrix(vegdist(T9.RO5.AK4.otu, "bray"))

T9_RO5.AK5.bc_dist = as.dist((vegdist(T9.RO5.AK5.otu, "bray")))
T9_RO5.AK5.bc_mat<-as.matrix(vegdist(T9.RO5.AK5.otu, "bray"))

T9_RO5.AK6.bc_dist = as.dist((vegdist(T9.RO5.AK6.otu, "bray")))
T9_RO5.AK6.bc_mat<-as.matrix(vegdist(T9.RO5.AK6.otu, "bray"))

#### RO6
T9_RO6.AK4.bc_dist = as.dist((vegdist(T9.RO6.AK4.otu, "bray")))
T9_RO6.AK4.bc_mat<-as.matrix(vegdist(T9.RO6.AK4.otu, "bray"))

T9_RO6.AK5.bc_dist = as.dist((vegdist(T9.RO6.AK5.otu, "bray")))
T9_RO6.AK5.bc_mat<-as.matrix(vegdist(T9.RO6.AK5.otu, "bray"))

T9_RO6.AK6.bc_dist = as.dist((vegdist(T9.RO6.AK6.otu, "bray")))
T9_RO6.AK6.bc_mat<-as.matrix(vegdist(T9.RO6.AK6.otu, "bray"))

#################
# Mantel Test
T9.RO4.AK4_mantel_test<-mantel.rtest(g_dists_m.T9.RO4.AK4, T9_RO4.AK4.bc_dist, nrepet=9999)
T9.RO4.AK4_mantel_test

T9.RO4.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.RO4.AK5, T9_RO4.AK5.bc_dist, nrepet=9999)
T9.RO4.AK5_mantel_test

T9.RO4.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.RO4.AK6, T9_RO4.AK6.bc_dist, nrepet=9999)
T9.RO4.AK6_mantel_test


### RO5
T9.RO5.AK4_mantel_test<-mantel.rtest(g_dists_m.T9.RO5.AK4, T9_RO5.AK4.bc_dist, nrepet=9999)
T9.RO5.AK4_mantel_test

T9.RO5.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.RO5.AK5, T9_RO5.AK5.bc_dist, nrepet=9999)
T9.RO5.AK5_mantel_test

T9.RO5.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.RO5.AK6, T9_RO5.AK6.bc_dist, nrepet=9999)
T9.RO5.AK6_mantel_test

### RO6
T9.RO6.AK4_mantel_test<-mantel.rtest(g_dists_m.T9.RO6.AK4, T9_RO6.AK4.bc_dist, nrepet=9999)
T9.RO6.AK4_mantel_test

T9.RO6.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.RO6.AK5, T9_RO6.AK5.bc_dist, nrepet=9999)
T9.RO6.AK5_mantel_test

T9.RO6.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.RO6.AK6, T9_RO6.AK6.bc_dist, nrepet=9999)
T9.RO6.AK6_mantel_test

mantelPower(T9.RO6.AK6_mantel_test, effect.size = 0.5)

###############
# make df
### RO4
T9_RO4.AK4.bcdf<-data.frame(Distance=geomat.T9.RO4.AK4[lower.tri(geomat.T9.RO4.AK4)],
                            BrayCurtis= T9_RO4.AK4.bc_mat[lower.tri(T9_RO4.AK4.bc_mat)])

T9_RO4.AK5.bcdf<-data.frame(Distance=geomat.T9.RO4.AK5[lower.tri(geomat.T9.RO4.AK5)],
                            BrayCurtis= T9_RO4.AK5.bc_mat[lower.tri(T9_RO4.AK5.bc_mat)])

T9_RO4.AK6.bcdf<-data.frame(Distance=geomat.T9.RO4.AK6[lower.tri(geomat.T9.RO4.AK6)],
                            BrayCurtis= T9_RO4.AK6.bc_mat[lower.tri(T9_RO4.AK6.bc_mat)])

### RO5
T9_RO5.AK4.bcdf<-data.frame(Distance=geomat.T9.RO5.AK4[lower.tri(geomat.T9.RO5.AK4)],
                            BrayCurtis= T9_RO5.AK4.bc_mat[lower.tri(T9_RO5.AK4.bc_mat)])

T9_RO5.AK5.bcdf<-data.frame(Distance=geomat.T9.RO5.AK5[lower.tri(geomat.T9.RO5.AK5)],
                            BrayCurtis= T9_RO5.AK5.bc_mat[lower.tri(T9_RO5.AK5.bc_mat)])

T9_RO5.AK6.bcdf<-data.frame(Distance=geomat.T9.RO5.AK6[lower.tri(geomat.T9.RO5.AK6)],
                            BrayCurtis= T9_RO5.AK6.bc_mat[lower.tri(T9_RO5.AK6.bc_mat)])

### RO6
T9_RO6.AK4.bcdf<-data.frame(Distance=geomat.T9.RO6.AK4[lower.tri(geomat.T9.RO6.AK4)],
                            BrayCurtis= T9_RO6.AK4.bc_mat[lower.tri(T9_RO6.AK4.bc_mat)])

T9_RO6.AK5.bcdf<-data.frame(Distance=geomat.T9.RO6.AK5[lower.tri(geomat.T9.RO6.AK5)],
                            BrayCurtis= T9_RO6.AK5.bc_mat[lower.tri(T9_RO6.AK5.bc_mat)])

T9_RO6.AK6.bcdf<-data.frame(Distance=geomat.T9.RO6.AK6[lower.tri(geomat.T9.RO6.AK6)],
                            BrayCurtis= T9_RO6.AK6.bc_mat[lower.tri(T9_RO6.AK6.bc_mat)])

# mantel test df
prox.mantT9<-as.data.frame(rbind(
  T9.RO4.AK4_mantel_test[1], T9.RO4.AK5_mantel_test[1],T9.RO4.AK6_mantel_test[1],
  T9.RO5.AK4_mantel_test[1], T9.RO5.AK5_mantel_test[1], T9.RO5.AK6_mantel_test[1],
  T9.RO6.AK4_mantel_test[1], T9.RO6.AK5_mantel_test[1], T9.RO6.AK6_mantel_test[1]))
prox.mantT9$obs<-as.numeric(prox.mantT9$obs)

prox.mantT9$pvalue<-as.numeric(rbind(
  T9.RO4.AK4_mantel_test[5], T9.RO4.AK5_mantel_test[5], T9.RO4.AK6_mantel_test[5],
  T9.RO5.AK4_mantel_test[5], T9.RO5.AK5_mantel_test[5], T9.RO5.AK6_mantel_test[5],
  T9.RO6.AK4_mantel_test[5], T9.RO6.AK5_mantel_test[5], T9.RO6.AK6_mantel_test[5]))

prox.mantT9$comparison<-as.factor(c("T9.RO4.AK4", "T9.RO4.AK5", "T9.RO4.AK6",
                                    "T9.RO5.AK4", "T9.RO5.AK5", "T9.RO5.AK6",
                                    "T9.RO6.AK4", "T9.RO6.AK5", "T9.RO6.AK6"))

prox.mantT9$AK.group<-as.factor(c("AK far", "AK mid", "AK near", 
                                  "AK far", "AK mid", "AK near", 
                                  "AK far", "AK mid", "AK near"))
prox.mantT9$Transect<-"9"

##############
####################### plots

#combined df
prox.mantT8.T9<-rbind(prox.mantT8, prox.mantT9)
prox.mantT8.T9$AK.group<-factor(prox.mantT8.T9$AK.group, levels = c("AK near", "AK mid", "AK far"))
write.csv(prox.mantT8.T9, "output/proximity.test.csv")

prox.mantT8.T9$Transect<-factor(prox.mantT8.T9$Transect, levels=c("9", "8")) # order for north to south

#rename to add north-south
prox.mantT8.T9$Transect<-revalue(prox.mantT8.T9$Transect,
                                 c("9"="9-northern",
                                 "8"="8-southern"))

mantel.proximity<-ggplot(prox.mantT8.T9, aes(x=AK.group, y=obs, color=Transect)) + 
  geom_boxplot(aes(color = Transect), width = 0.5, size = 0.4,
               position = position_dodge(0.7)) +
  geom_dotplot(aes(fill = Transect, color = Transect),
               binaxis='y', stackdir='center', dotsize = 0.5, alpha=0.5,
               position = position_dodge(0.7))+
  scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
  scale_color_manual(values=c("#00AFBB", "#E7B800"))+
  ylab("Kendall rank correlation coefficient (r)") +
  xlab("Remnant-to-Restored Forest Proximity") +
  theme(axis.text.x=element_text(colour="black",size=7), 
        axis.text.y=element_text(colour="black",size=7),
        axis.title=element_text(colour="black",size=9),
        legend.key=element_blank())+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

mantel.proximity
dev.copy(png, "figures/mantel.proximity.png", width=1400, height=1500, res=300)
dev.off()



##################################################################
############################################
######################
# Test Mantel distance between the 3 AKS in a single transect pairwise
# Just AKs 
######################
################# Transect 8
# plot distance: W --> E (AK1, AK2, AK3, RO4, RO5, RO6)


# T8 dat = transposed OTUs
# sample.T8 = samples
# T8.df = dataframe 
# T8_meta = the meta data
# T8_otu = otus

## T8 single plot dataframes
T8.AK1.df<-T8.df[(T8.df$Plot=="AK1"),]; T8.AK1.meta<- T8.AK1.df[,c(1:34)]; T8.AK1.otu<-T8.AK1.df[,c(35:1470)]
T8.AK2.df<-T8.df[(T8.df$Plot=="AK2"),]; T8.AK2.meta<- T8.AK2.df[,c(1:34)]; T8.AK2.otu<-T8.AK2.df[,c(35:1470)]
T8.AK3.df<-T8.df[(T8.df$Plot=="AK3"),]; T8.AK3.meta<- T8.AK3.df[,c(1:34)]; T8.AK3.otu<-T8.AK3.df[,c(35:1470)]

# combine for comparisons
T8.AK1.AK2.df<-rbind(T8.AK1.df,T8.AK2.df)
T8.AK1.AK2.meta <-rbind(T8.AK1.meta,T8.AK2.meta)
T8.AK1.AK2.otu<-rbind(T8.AK1.otu,T8.AK2.otu)

T8.AK1.AK3.df<-rbind(T8.AK1.df,T8.AK3.df)
T8.AK1.AK3.meta <-rbind(T8.AK1.meta,T8.AK3.meta)
T8.AK1.AK3.otu<-rbind(T8.AK1.otu,T8.AK3.otu)

T8.AK2.AK3.df<-rbind(T8.AK2.df,T8.AK3.df)
T8.AK2.AK3.meta <-rbind(T8.AK2.meta,T8.AK3.meta)
T8.AK2.AK3.otu<-rbind(T8.AK2.otu,T8.AK3.otu)


# Distance to neighboring plots
######### 
# AK1 - AK2
g_dists_km.T8.AK1.AK2<- earth.dist(T8.AK1.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.AK1.AK2<- (g_dists_km.T8.AK1.AK2 * 1000)
g_dists_m.T8.AK1.AK2 <- as.matrix(g_dists_m.T8.AK1.AK2)
rownames(g_dists_m.T8.AK1.AK2) <- rownames(T8.AK1.AK2.meta)
colnames(g_dists_m.T8.AK1.AK2)<- t(rownames(T8.AK1.AK2.meta))
# rename and make distance matrix
geomat.T8.AK1.AK2<-g_dists_m.T8.AK1.AK2
g_dists_m.T8.AK1.AK2<-as.dist(g_dists_m.T8.AK1.AK2) 


#### AK1-AK3
g_dists_km.T8.AK1.AK3<- earth.dist(T8.AK1.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.AK1.AK3<- (g_dists_km.T8.AK1.AK3 * 1000)
g_dists_m.T8.AK1.AK3 <- as.matrix(g_dists_m.T8.AK1.AK3)
rownames(g_dists_m.T8.AK1.AK3) <- rownames(T8.AK1.AK3.meta)
colnames(g_dists_m.T8.AK1.AK3)<- t(rownames(T8.AK1.AK3.meta))
# rename and make distance matrix
geomat.T8.AK1.AK3<-g_dists_m.T8.AK1.AK3
g_dists_m.T8.AK1.AK3<-as.dist(g_dists_m.T8.AK1.AK3) 


#### AK2-AK3
g_dists_km.T8.AK2.AK3<- earth.dist(T8.AK2.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.AK2.AK3<- (g_dists_km.T8.AK2.AK3 * 1000)
g_dists_m.T8.AK2.AK3 <- as.matrix(g_dists_m.T8.AK2.AK3)
rownames(g_dists_m.T8.AK2.AK3) <- rownames(T8.AK2.AK3.meta)
colnames(g_dists_m.T8.AK2.AK3)<- t(rownames(T8.AK2.AK3.meta))
# rename and make distance matrix
geomat.T8.AK2.AK3<-g_dists_m.T8.AK2.AK3
g_dists_m.T8.AK2.AK3<-as.dist(g_dists_m.T8.AK2.AK3) 


############################
#Calculate Bray-Curtis distances on T8 between RO to AK

#### AK1-2-3
T8_AK1.AK2.bc_dist = as.dist((vegdist(T8.AK1.AK2.otu, "bray")))
T8_AK1.AK2.bc_mat<-as.matrix(vegdist(T8.AK1.AK2.otu, "bray"))

T8_AK1.AK3.bc_dist = as.dist((vegdist(T8.AK1.AK3.otu, "bray")))
T8_AK1.AK3.bc_mat<-as.matrix(vegdist(T8.AK1.AK3.otu, "bray"))

T8_AK2.AK3.bc_dist = as.dist((vegdist(T8.AK2.AK3.otu, "bray")))
T8_AK2.AK3.bc_mat<-as.matrix(vegdist(T8.AK2.AK3.otu, "bray"))


#################
# Mantel Test
T8.AK1.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.AK1.AK2, T8_AK1.AK2.bc_dist, nrepet=9999)
T8.AK1.AK2_mantel_test

T8.AK1.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.AK1.AK3, T8_AK1.AK3.bc_dist, nrepet=9999)
T8.AK1.AK3_mantel_test

T8.AK2.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.AK2.AK3, T8_AK2.AK3.bc_dist, nrepet=9999)
T8.AK2.AK3_mantel_test


######################
######################
######################
################# Transect 9
# Just AKs 
# plot distance: W --> E (AK4, AK5, AK6, RO4, RO5, RO6)


# T9 dat = transposed OTUs
# sample.T9 = samples
# T9.df = dataframe 
# T9_meta = the meta data
# T9_otu = otus

## T9 single plot dataframes
T9.AK4.df<-T9.df[(T9.df$Plot=="AK4"),]; T9.AK4.meta<- T9.AK4.df[,c(1:34)]; T9.AK4.otu<-T9.AK4.df[,c(35:1470)]
T9.AK5.df<-T9.df[(T9.df$Plot=="AK5"),]; T9.AK5.meta<- T9.AK5.df[,c(1:34)]; T9.AK5.otu<-T9.AK5.df[,c(35:1470)]
T9.AK6.df<-T9.df[(T9.df$Plot=="AK6"),]; T9.AK6.meta<- T9.AK6.df[,c(1:34)]; T9.AK6.otu<-T9.AK6.df[,c(35:1470)]

# combine for comparisons
T9.AK4.AK5.df<-rbind(T9.AK4.df,T9.AK5.df)
T9.AK4.AK5.meta <-rbind(T9.AK4.meta,T9.AK5.meta)
T9.AK4.AK5.otu<-rbind(T9.AK4.otu,T9.AK5.otu)

T9.AK4.AK6.df<-rbind(T9.AK4.df,T9.AK6.df)
T9.AK4.AK6.meta <-rbind(T9.AK4.meta,T9.AK6.meta)
T9.AK4.AK6.otu<-rbind(T9.AK4.otu,T9.AK6.otu)

T9.AK5.AK6.df<-rbind(T9.AK5.df,T9.AK6.df)
T9.AK5.AK6.meta <-rbind(T9.AK5.meta,T9.AK6.meta)
T9.AK5.AK6.otu<-rbind(T9.AK5.otu,T9.AK6.otu)


# Distance to neighboring plots
######### 
# AK4 - AK5
g_dists_km.T9.AK4.AK5<- earth.dist(T9.AK4.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.AK4.AK5<- (g_dists_km.T9.AK4.AK5 * 1000)
g_dists_m.T9.AK4.AK5 <- as.matrix(g_dists_m.T9.AK4.AK5)
rownames(g_dists_m.T9.AK4.AK5) <- rownames(T9.AK4.AK5.meta)
colnames(g_dists_m.T9.AK4.AK5)<- t(rownames(T9.AK4.AK5.meta))
# rename and make distance matrix
geomat.T9.AK4.AK5<-g_dists_m.T9.AK4.AK5
g_dists_m.T9.AK4.AK5<-as.dist(g_dists_m.T9.AK4.AK5) 


#### AK4-AK6
g_dists_km.T9.AK4.AK6<- earth.dist(T9.AK4.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.AK4.AK6<- (g_dists_km.T9.AK4.AK6 * 1000)
g_dists_m.T9.AK4.AK6 <- as.matrix(g_dists_m.T9.AK4.AK6)
rownames(g_dists_m.T9.AK4.AK6) <- rownames(T9.AK4.AK6.meta)
colnames(g_dists_m.T9.AK4.AK6)<- t(rownames(T9.AK4.AK6.meta))
# rename and make distance matrix
geomat.T9.AK4.AK6<-g_dists_m.T9.AK4.AK6
g_dists_m.T9.AK4.AK6<-as.dist(g_dists_m.T9.AK4.AK6) 


#### AK5-AK6
g_dists_km.T9.AK5.AK6<- earth.dist(T9.AK5.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.AK5.AK6<- (g_dists_km.T9.AK5.AK6 * 1000)
g_dists_m.T9.AK5.AK6 <- as.matrix(g_dists_m.T9.AK5.AK6)
rownames(g_dists_m.T9.AK5.AK6) <- rownames(T9.AK5.AK6.meta)
colnames(g_dists_m.T9.AK5.AK6)<- t(rownames(T9.AK5.AK6.meta))
# rename and make distance matrix
geomat.T9.AK5.AK6<-g_dists_m.T9.AK5.AK6
g_dists_m.T9.AK5.AK6<-as.dist(g_dists_m.T9.AK5.AK6) 


############################
#Calculate Bray-Curtis distances on T9 between RO to AK

#### AK4-5-6
T9_AK4.AK5.bc_dist = as.dist((vegdist(T9.AK4.AK5.otu, "bray")))
T9_AK4.AK5.bc_mat<-as.matrix(vegdist(T9.AK4.AK5.otu, "bray"))

T9_AK4.AK6.bc_dist = as.dist((vegdist(T9.AK4.AK6.otu, "bray")))
T9_AK4.AK6.bc_mat<-as.matrix(vegdist(T9.AK4.AK6.otu, "bray"))

T9_AK5.AK6.bc_dist = as.dist((vegdist(T9.AK5.AK6.otu, "bray")))
T9_AK5.AK6.bc_mat<-as.matrix(vegdist(T9.AK5.AK6.otu, "bray"))


#################
# Mantel Test
T9.AK4.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.AK4.AK5, T9_AK4.AK5.bc_dist, nrepet=9999)
T9.AK4.AK5_mantel_test

T9.AK4.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.AK4.AK6, T9_AK4.AK6.bc_dist, nrepet=9999)
T9.AK4.AK6_mantel_test

T9.AK5.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.AK5.AK6, T9_AK5.AK6.bc_dist, nrepet=9999)
T9.AK5.AK6_mantel_test

#####################################################################
#####################################################################

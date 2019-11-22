###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

###### Hakalau Metacommunity 
###### Code and heavy lifting by Cameron Egan
###### "Calvary" coding by Chris Wall
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

pacman::p_load("ade4", "multtest","phyloseq","rhdf5","ggplot2","colorspace","stringi", "geosphere", 
               "ggplot2", "ggmap", "dplyr", "gridExtra", "geosphere", "sf", "raster", "spData",
               "tmap", "leaflet", "mapview", "shiny", "fossil", "RgoogleMaps", "devtools", "ggsn", "vegan", "multcomp",
               "dplyr", "grid", "scales", "gridExtra", "emmeans", "multcompView", "ggpubr", "Rmisc", "purrr",
               "RVAideMemoire", "RColorBrewer", "vegan")


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
hakalau_map_zoom <-get_map(location=c(-155.320,19.83),zoom=14,maptype="satellite",color="color")

haka_metadata <-read.csv("data/haka_soil_metadata.csv", header=TRUE, row.names=1)

sampling_plots <- ggmap(hakalau_map_zoom) + 
  geom_point(data=haka_metadata,aes(x=lon,y=lat,fill=HabitatType),pch=21, stroke=0.3,colour="white",size=3) +
  scale_fill_manual(values=c("#88A550","#336B87")) +
  theme(legend.text=element_text(size=6),legend.title = element_text(size=8), legend.key=element_blank()) +
  xlab("Longitude") + ylab("Latitude") +
  scale_y_continuous(limits=c(19.81, 19.84)) +
  scale_x_continuous(limits=c(-155.345, -155.295)) +
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

#Calculate geographic distances between samples
library(fossil)
haka_gps <- cbind(new_lat_lon)
geog_dists_km<-earth.dist(haka_gps[,1:2],dist=TRUE)
geog_dists_m <- (geog_dists_km * 1000)
geog_dists_m <- as.matrix(geog_dists_m)
rownames(geog_dists_m) <- rownames(haka_meta)
colnames(geog_dists_m)<- t(rownames(haka_meta))

write.csv(geog_dists_m,file="output/haka_dists.csv")

haka_dists <- lower.tri(geog_dists_m)

         

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

sample$Hab.by.hot<-interaction(sample$HabitatType, sample$Host)

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

spec_rich_host = plot_richness(haka_VT_soil_physeq, x ="Host", measures="Observed", color="Host")  +
  geom_boxplot(col="black", aes(fill=Host), alpha=0.8 , lwd=0.5, outlier.colour = "gray50") +
  scale_color_manual(name= "Host species", values=c("#D73027","#FC8D59","goldenrod","#008000","darkslategray3","#008B8B","#4575B4")) +
  scale_fill_manual(name= "Host species", values=c("#D73027","#FC8D59","goldenrod","#008000","darkslategray3","#008B8B","#4575B4"))  + 
  xlab("Habitat Type") + ylab("AM fungal richness") + 
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(colour="black",size=8)) +
  scale_y_continuous(breaks=seq(0,75,by=10),limits=c(0,75)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  facet_wrap(~HabitatType) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(strip.text.x=element_text(size=10,face="bold"),strip.background = element_rect(fill="gray90")) +
  theme(legend.text = element_text(face="italic")) 

plot(spec_rich_host)
ggsave("figures/host_spec_rich.tiff", plot = spec_rich_host, width=7,height=5)


##By habitat type alone
spec_rich_hab = plot_richness(haka_VT_soil_physeq, x ="Plot", measures="Observed", color="HabitatType") +
  geom_boxplot(col="black", aes(fill=HabitatType), alpha=0.8 , lwd=0.5) + 
  scale_color_manual(name= "Habitat Type", values=c("#88A550", "#336B87")) +
  scale_fill_manual(name= "Habitat Type", values=c("#88A550", "#336B87")) +
  scale_x_discrete(limits=c("RO1", "RO2", "RO3", "RO4", "RO5", "RO6",
                            "AK1", "AK2", "AK3", "AK4", "AK5", "AK6")) +
  ylim(0,70) + 
  theme(text=element_text(colour="black",size=10)) + 
  ylab("AM fungal richness") + xlab("Habitat Type Plots") +
  theme(axis.text.x=element_text(angle=0,hjust=0.5, vjust=0.5,colour="black",size=8)) +
  theme(axis.text.y=element_text(colour="black",size=8)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  facet_wrap(~SampleType) +
  theme(strip.text.x=element_text(size=8,face="bold"),strip.background = element_rect(fill="white")) 

plot(spec_rich_hab)
ggsave("figures/species_richness_by_plot.tiff", plot = spec_rich_hab, width=6,height=5)


#########################################################
# export richness from soil data and test for differences
soil_rich<-estimate_richness(haka_VT_soil_physeq, measures="Observed")
pairwise.wilcox.test(soil_rich$Observed, sample_data(haka_VT_soil_physeq)$Host) # not different
pairwise.wilcox.test(soil_rich$Observed, sample_data(haka_VT_soil_physeq)$Habitat) # different 

# inspect alpha diversity a bit more
soil_rich$sampleID<-as.factor(rownames(soil_rich))
soil_rich$HabitatType <- substr(soil_rich$sampleID, 0, 2) # extract first 2 letters of ID
soil_rich$Plot <- substr(soil_rich$sampleID, 3, 3)
soil_rich$Host <- substr(soil_rich$sampleID, 4, 5)
soil_rich<-na.omit(soil_rich)
soil_rich$HabitatType<-revalue(soil_rich$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest"))

# how does alpha diverstiy differ among Habitat Types? (32 AK vs. 30 RO, so minor...)
soil_rich %>% 
  group_by(HabitatType) %>%
  dplyr::summarise(mean = mean(Observed))

#GLM with poisson error distribution to test for differences in Species richness in soil 
soil_glm <- glm(Observed ~ HabitatType*Host, data=soil_rich, family = poisson(link="log")) 
par(mfrow = c(2, 2))
plot(soil_glm)
chi.anova<-anova(soil_glm ,test="Chisq")
chi.anova # Habitat:Host interaction

#Create a unique level for every combination of site and treatment to do a post-hoc test on
soil_rich$Hab.by.Host<- interaction(soil_rich$HabitatType, soil_rich$Host)
soil_glm.2 <- glm(Observed ~ Hab.by.Host, data=soil_rich, family = poisson(link="log"))

#Calculate EMM on interactions
rich_emmeans <- emmeans(soil_glm.2, specs="Hab.by.Host")
rich_posthoc.pairs = pairs(rich_emmeans)

#Create letters for interaction differences
rich_mc_letters<-cld(rich_emmeans,Letters="abcdefg")
rich_mc_letters

# What differences exist with species and habitat interaction?
# less alpha diversity for Cheirodendron trigynum (CH) in RO (i.e, RO.CH) compared to other groups.
# slightly higher in Metrosideros polymorpha from RO (i.e., RO.ME) 
# slightly highest in Grass from AK (i.e., AK.GR)

#########################################################
#########################################################


         ######################
         ### Beta Diversity ###
         ######################

###Phyloseq###
#Bray plot
bc_dist = as.matrix((vegdist(haka_otu, "bray")))
NMDS = metaMDS(bc_dist)
NMDS1=NMDS$points[,1]
NMDS2=NMDS$points[,2]
NMDS.plot.df=data.frame(NMDS1=NMDS1,NMDS2=NMDS2, Host=sample$Host, 
                     HabitatType=sample$HabitatType,Plot=sample$Plot)
NMDS.plot.df$Sp.Hab<-interaction(NMDS.plot.df$Host, NMDS.plot.df$HabitatType)

NMDS.col<- c("#88A550","#336B87") # colors
Habitats<- NMDS.plot.df$HabitatType # habitat levels
groups.sp<-c("coral", "seagreen", "mediumpurple", "goldenrod", "dodgerblue", "gray50", "navajowhite2")
groups.hab<-c("#88A550","#336B87")


### make points semitransparent
# green
col2rgb("#88A550"); grn.60<-rgb(136, 165, 80, max=255, alpha=170)
# blue
col2rgb("#336B87"); blu.60<-rgb(51, 107, 135, max=255, alpha=170)
Hab.col<-c(grn.60, blu.60)

##### PLOTS

###### NMDS plot by Species
ordiplot(NMDS, type="n", main=substitute(paste("Species NMDS")), cex.main=1, display="sites", xlim=c(-0.25, 0.8), cex.lab=0.8, cex.axis=0.8)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(NMDS, "sites", cex=0.8, pch=16, col=groups.sp[NMDS.plot.df$Host])
ordihull(NMDS, groups=NMDS.plot.df$Host, draw="polygon", alpha=20, col=groups.sp, border=groups.sp)
#ordiellipse(NMDS, groups=NMDS.plot.df$Host, draw="polygon", kind="sd", alpha=20, conf=0.95, col=groups.sp, border=groups.sp)
legend("topright", legend=levels(NMDS.plot.df$Host),  text.font=3, cex=1, pch=16, col=groups.sp, pt.cex=1, bty="n")

dev.copy(pdf, "figures/species.NMDS.pdf", height=8, width=8)
dev.off() 


###### NMDS plot by Habitat
ordiplot(NMDS, type="n", main=substitute(paste("Habitat NMDS")), cex.main=1, display="sites", xlim=c(-0.25, 0.8), cex.lab=0.8, cex.axis=0.8)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(NMDS, "sites", cex=0.8, pch=16, col=Hab.col[Habitats])
ordiellipse(NMDS, groups=Habitats, kind="sd", draw="polygon", conf=0.95, alpha=30, col=groups.hab, border=groups.hab)
legend("topright", legend=levels(Habitats), cex=1, pch=16, col=groups.hab, pt.cex=1, bty="n")

dev.copy(pdf, "figures/habitat.NMDS.pdf", height=8, width=8)
dev.off() 


###### NMDS plot of Species x Habitat
NMDS_host_plot = ggplot(NMDS.plot.df, aes(x=NMDS1,y=NMDS2,shape=HabitatType)) +
  geom_point(size=1.5, stroke=1, alpha=0.7, aes(fill=HabitatType,colour=Host)) +
  scale_fill_manual(values=c("black","white")) +
  scale_shape_manual(values=c(16,16)) +
  stat_ellipse(size=0.8,alpha=0.75,aes(color=Host,linetype=HabitatType),type='t') +
  scale_linetype_manual(values=c(1,6)) +
  scale_color_manual(values=c("#D73027","#FC8D59","#FEE090","#008000","#E0F3F8","#008B8B","#4575B4")) +
  theme(axis.text.x=element_text(colour="black",size=10)) +
  theme(axis.text.y=element_text(colour="black",size=10)) +
  theme(axis.title=element_text(colour="black",size=10)) +
  theme(panel.border = element_rect(fill=NA),
        panel.background=element_blank()) +
  theme(strip.text.x=element_text(size=10,face="bold"),strip.background = element_rect(fill=NA)) +
  theme(legend.text=element_text(size=8, face='italic')) +
  theme(legend.title=element_text(size=8,face='bold')) 

plot(NMDS_host_plot)
ggsave("figures/NMDS_host_plot.tiff",width= 8,height=8, plot=NMDS_host_plot)


##### ##### ##### ##### ##### ##### #####

#Examination of environmental data
environmental_data <- read.csv("data/haka_chemistry.csv", header=TRUE, row.names=1)
all.equal(rownames(haka_otu), rownames(environmental_data))
rows_to_keep<-rownames(haka_otu)
environmental_data<-environmental_data[rows_to_keep,]
all.equal(rownames(haka_otu), rownames(environmental_data))

#Fit environmental vectors to ordination to see which environmental variables are correlated with the ordination

colnames(environmental_data)<-c("OM(%)","Total N", "P", "K", "Mg", "Ca", "Na", "S","pH", "H(meq/100g)", "CEC(meq/100g)", "K+", "Mg+2", "Ca+2", "H+", "Na+")


##########  NMDS

fit.env <- envfit(NMDS, environmental_data, na.rm=TRUE)

######### make plot by habitat with vectors for environment
environm_plot<-ordiplot(NMDS, type="n", main=substitute(paste("Habitat NMDS")), cex.main=1, display="sites", xlim=c(-0.25, 0.8), cex.lab=0.8, cex.axis=0.8)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(NMDS, "sites", cex=0.8, pch=16, col=Hab.col[Habitats])
ordihull(NMDS, groups=Habitats, draw="polygon", alpha=20, col=groups.hab, border=groups.hab)
legend("topright", legend=levels(Habitats), cex=1, pch=16, col=groups.hab, pt.cex=1, bty="n")
par.new=T
plot(fit.env, col="black", p.max=0.05, cex=0.9, lwd=1)


dev.copy(pdf, "figures/environm.NMDS.pdf", height=8, width=8)
dev.off() 


##### breakup data

######## Only RO site
ROdat<-t(otu[,grep("^RO", colnames(otu))]) # transposed data with only RO
samples.RO<-sample[grep("^RO", rownames(sample)),] # only rows with "RO"
all.equal(rownames(ROdat), rownames(samples.RO)) # rows match


RO_dist = as.matrix((vegdist(ROdat, "bray")))
RO.NMDS = metaMDS(RO_dist)
RO.NMDS1=RO.NMDS$points[,1]
RO.NMDS2=RO.NMDS$points[,2]
RO.NMDS.plot.df=data.frame(NMDS1=RO.NMDS1,NMDS2=RO.NMDS2, Host=samples.RO$Host, 
                        HabitatType=samples.RO$HabitatType,Plot=samples.RO$Plot)


###### RESTORED NMDS plot by Species
RO.species_plot<-ordiplot(RO.NMDS, type="n", main=substitute(paste("Restored Forest: Species NMDS")), cex.main=1, display="sites", xlim=c(-0.7, 0.8), cex.lab=0.8, cex.axis=0.8)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(RO.NMDS, "sites", cex=0.8, pch=16, col=groups.sp[RO.NMDS.plot.df$Host])
ordihull(RO.NMDS, groups=RO.NMDS.plot.df$Host, draw="polygon", alpha=20, col=groups.sp, border=groups.sp)
#ordiellipse(NMDS, groups=NMDS.plot.df$Host, draw="polygon", kind="sd", alpha=20, conf=0.95, col=groups.sp, border=groups.sp)
legend("topright", legend=levels(RO.NMDS.plot.df$Host),  text.font=3, cex=1, pch=16, col=groups.sp, pt.cex=1, bty="n")

dev.copy(pdf, "figures/RO.NMDS.pdf", height=8, width=8)
dev.off() 

haka.RO.bc.adonis <- adonis(RO_dist~Host, data=samples.RO, permutations = 9999)
haka.RO.bc.adonis


######## Only AK site
AKdat<-t(otu[,grep("^AK", colnames(otu))]) # transposed data with only AK
samples.AK<-sample[grep("^AK", rownames(sample)),] # only AKws with "AK"
all.equal(rownames(AKdat), rownames(samples.AK)) # AKws match


AK_dist = as.matrix((vegdist(AKdat, "bray")))
AK.NMDS = metaMDS(AK_dist)
AK.NMDS1=AK.NMDS$points[,1]
AK.NMDS2=AK.NMDS$points[,2]
AK.NMDS.plot.df=data.frame(NMDS1=AK.NMDS1,NMDS2=AK.NMDS2, Host=samples.AK$Host, 
                        HabitatType=samples.AK$HabitatType, Plot=samples.AK$Plot)
 

###### REMNANT NMDS plot by Species
AK.species_plot<-ordiplot(AK.NMDS, type="n", main=substitute(paste("Remnant Forest: Species NMDS")), cex.main=1, display="sites", xlim=c(-0.7, 1.1), cex.lab=0.8, cex.axis=0.8)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(AK.NMDS, "sites", cex=0.8, pch=16, col=groups.sp[AK.NMDS.plot.df$Host])
ordihull(AK.NMDS, groups=AK.NMDS.plot.df$Host, draw="polygon", alpha=20, col=groups.sp, border=groups.sp)
#ordiellipse(NMDS, groups=NMDS.plot.df$Host, draw="polygon", kind="sd", alpha=20, conf=0.95, col=groups.sp, border=groups.sp)
legend("topright", legend=levels(AK.NMDS.plot.df$Host),  text.font=3, cex=1, pch=16, col=groups.sp, pt.cex=1, bty="n")

dev.copy(pdf, "figures/AK.NMDS.pdf", height=8, width=8)
dev.off() 

haka.AK.bc.adonis <- adonis(AK_dist~Host, data=samples.AK, permutations = 9999)
haka.AK.bc.adonis


###Stats###
#Bray PERMANOVA 
haka.bc.adonis <- adonis(bc_dist~HabitatType*Host, data=sample, permutations = 9999)
haka.bc.adonis

#Pairwise permanova
haka.pair.perm <- pairwise.perm.manova(bc_dist, sample$Hab.by.host,
                                       nperm=9999,p.method="fdr")
haka.pair.perm

##Bray Betadisper
#Total difference
haka.beta.disper <- betadisper(bc_dist, root_metadata$Hab.by.Host)
haka.beta.disper.results <- permutest(haka.beta.disper, pairwise = TRUE, iter=9999)
haka.beta.disper.results


######## ######## ######## ######## ######## Make environmental data as PC1 and PC2
######## explore environmental data, use PCA to visualize
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
require(graphics)
library(plyr)

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

####### by HabitatType
HabitatType<-df.PCA$HabitatType
PC.habtype.fig <- ggbiplot(Hak.env.PCA, choices = 1:2, obs.scale = 1, var.scale = 1, 
                           groups= HabitatType, ellipse = TRUE, ellipse.prob = 0.90,
                           circle = FALSE, alpha=0, PC.Site.fig) +
  scale_color_manual(name = '', values=NMDS.col) +
  geom_point(aes(colour=HabitatType), shape=17, size = 1, alpha=6/10)+
  theme_bw() +
  theme(axis.ticks.length=unit(-0.25, "cm"), axis.text.y=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")), axis.text.x=element_text(margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))) +
  theme(legend.text=element_text(size=10)) +
  theme(panel.background = element_rect(colour = "black", size=1))+
  theme(legend.key = element_blank())+
  theme(legend.direction = 'horizontal', legend.position = 'top') + theme(aspect.ratio=0.7)

PC.habtype.fig
dev.copy(pdf, "figures/environm.PCA.pdf", height=5, width=6)
dev.off() 


###############
###############
NMDS.otu<-as.data.frame(NMDS1) # NMDS Bray Curtis of OTUs
PCA.env<-as.data.frame(newdat) # the exported data from PCA of environment
PC.NMD<-merge(PCA.env, NMDS.otu, by = "row.names", all = TRUE) # merge dataframes
PC.NMD<-na.omit(PC.NMD) #drop NA columns that don't correspond
colnames(PC.NMD)[1]<-"sampleID"
PC.NMD$HabitatType <- as.factor(substr(PC.NMD$sampleID, 0, 2)) # extract first 2 letters of ID
PC.NMD$HabitatType<-revalue(PC.NMD$HabitatType, c("AK"="Restored Forest", "RO"="Remnant Forest"))
PC.NMD$Plot <- as.factor(substr(PC.NMD$sampleID, 3, 3)) # make plot ID
PC.NMD$Host <- as.factor(substr(PC.NMD$sampleID, 4, 5)) # make host ID

PC.NMD$HabitatType<-factor(PC.NMD$HabitatType, levels=c("Remnant Forest", "Restored Forest"))
ggplot(PC.NMD, aes(x=NMDS1, y=PC1)) + geom_point(size=1,alpha=0.5, aes(color=HabitatType)) +
  scale_color_manual(values=NMDS.col) +
  stat_smooth(method = "lm", size = 1, se=T, col="orchid")+
  xlab("NMDS1 (ESV Bray-Curtis)") +
  ylab("PC1 (environment)") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())
dev.copy(pdf, "figures/PC.NMD.pdf", height=5, width=6)
dev.off()

mod<-lm(PC1~NMDS1:HabitatType, data=PC.NMD)
int<-coef(summary(mod))[1]
slope.AK<-coef(summary(mod))[2] # NMDS1:HabitatTypeRestored Forest
slope.RO<-coef(summary(mod))[3] #NMDS1:HabitatTypeRemnant Forest

library(plotrix)
plot(PC1~NMDS1, data=PC.NMD, pch=16,
     col=Hab.col[as.factor(HabitatType)],
     ylab="PC1 (environment)", 
     xlab="NMDS1 (ESV Bray-Curtis)")
ablineclip(int, slope.RO, col=NMDS.col[2], lwd=2, 
           x1 = min(PC.NMD$NMDS1[PC.NMD$HabitatType=="Remnant Forest"], na.rm=T), 
           x2 = max(PC.NMD$NMDS1[PC.NMD$HabitatType=="Remnant Forest"], na.rm=T)) # RK model
ablineclip(int, slope.AK, col=NMDS.col[1], lwd=2, 
           x1 = min(PC.NMD$NMDS1[PC.NMD$HabitatType=="Restored Forest"], na.rm=T), 
           x2 = max(PC.NMD$NMDS1[PC.NMD$HabitatType=="Restored Forest"], na.rm=T)) # AK model
legend("topleft", c("Remnant Forest", "Restored Forest"), lty=c(1,1), lwd=c(2,2), col=NMDS.col, cex=0.8, pch=16, y.intersp = 0.5, bty="n")
dev.copy(pdf, "figures/PC.NMD.slopes.pdf", height=5, width=6)
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
Study.dat<-Study.dat[!is.na(Study.dat[34]),]

#separate dataframes
haka_met_mant<-Study.dat[,c(1:32)] # the meta data
hak_otu<-Study.dat[,c(33:1468)] # OTU data

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
mantel_test<-mantel.rtest(dists_m, study_bc_dist, nrepet=9999)
mantel_test # significant

# new dataframe
Study.df<-data.frame(Distance=geomat[lower.tri(geomat)],
                  BrayCurtis= study_bc_mat[lower.tri(study_bc_mat)])

##############
#Generate plots for all samples combined
Study_dist_plot<- ggplot(Study.df, aes(x=log(Distance+1), y=BrayCurtis)) +
  geom_point(size=1,alpha=0.5) +
  stat_smooth(method = "lm", size = 1, se=F) +  
  theme(text=element_text(colour="black",size=15)) + 
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=12)) +
  theme(axis.text.y=element_text(colour="black",size=12)) +
  theme(legend.title=element_text(colour="black",size=12,face="bold")) +
  theme(legend.text=element_text(colour="black",size=12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

plot(Study_dist_plot)
ggsave("figures/log.study_dist_decay.tiff", 
       plot = Study_dist_plot, width = 5, height = 5)


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
RO<-RO[!is.na(RO[34]),]

RO_meta<-RO[,c(1:32)] # the meta data
RO_otu<-RO[,c(33:1468)] # OTU data

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
AK<-AK[!is.na(AK[34]),]

AK_meta<-AK[,c(1:32)] # the meta data
AK_otu<-AK[,c(33:1468)] # OTU data

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
AK_mantel_test <- mantel.rtest(AK_dists_m,AK_bc_dist,nrepet=9999)
AK_mantel_test
AK_df<-data.frame(Distance=AK_geomat[lower.tri(AK_geomat)],
                    BrayCurtis=AK_bc_dist[lower.tri(AK_bc_dist)])
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

RO_dist_plot<- ggplot(RO_df,aes(x=log(Distance+1), y=BrayCurtis)) + #color=Host
  #scale_color_manual(values=c(A.millefolium="#0A191E",D.fruticosa="#D8B65C",F.idahoensis="#4A9878")) +
  geom_point(size=1,alpha=0.5) +
  stat_smooth(method = "lm", size = 1, se=F) +  
  theme(text=element_text(colour="black",size=12)) + 
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=12)) +
  theme(axis.text.y=element_text(colour="black",size=12)) +
  theme(legend.title=element_text(colour="black",size=12,face="bold")) +
  theme(legend.text=element_text(colour="black",size=12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

plot(RO_dist_plot)
ggsave("figures/log.RO_dist_decay.pdf", 
       plot = RO_dist_plot, width = 5, height = 5)


AK_dist_plot<- ggplot(AK_df,aes(x=log(Distance+1),y=BrayCurtis))+ #color=Host
  #scale_color_manual(values=c(A.millefolium="#0A191E",D.fruticosa="#D8B65C",F.idahoensis="#4A9878")) +
  geom_point(size=1,alpha=0.5) +
  stat_smooth(method = "lm", size = 1, se=F) +  
  theme(text=element_text(colour="black",size=12)) + 
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=12)) +
  theme(axis.text.y=element_text(colour="black",size=12)) +
  theme(legend.title=element_text(colour="black",size=12,face="bold")) +
  theme(legend.text=element_text(colour="black",size=12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

plot(AK_dist_plot)
ggsave("figures/log.AK_dist_decay.pdf", 
       plot = AK_dist_plot, width = 5, height = 5)

#Analysis of slopes
# ALL hosts, RO soil
RO_dist_model <- lm(BrayCurtis~poly(Distance,2,raw=TRUE),data=RO_dist_df)
summary(RO_dist_model)
par(mfrow = c(2,2))
plot(RO_dist_model)

# ALL hosts, AK roots
AK_dist_model <- lm(BrayCurtis~poly(Distance,2,raw=TRUE),data=AK_dist_df)
summary(AK_dist_model)
par(mfrow = c(2,2))
plot(AK_dist_model)


#Calculate the difference in slop between overall regression lines 
library(simba)
diff_in_slope<-diffslope(RO_dist_df$Distance,RO_dist_df$BrayCurtis,AK_dist_df$Distance,AK_dist_df$BrayCurtis,
                         permutations=9999)
diff_in_slope



 #########################################
 ### MetaCommunity Simulations (MCSim) ###
 #########################################
# From tutorial: http://rstudio-pubs-static.s3.amazonaws.com/159425_80725873417e42fdb13821c10a198281.html
# -- Install the current dev version of MCSim
install.packages("devtools")
devtools::install_github('sokole/MCSim')
devtools::install_github('sokole/MCSim@v0.4.1.9001')
install.packages("MCSim"); library(MCSim)

# 1. Make a "landscape"
# The landscape is the "game board" on which the simulation plays out
# It is created using the fn.make.landscape function. 
# The landscape can be made from ESV_dataframe (non spatial points version).
# When making a landscape, we have to embed additional information about the simulation in the landscape 
# so that MCSim can keep track of site characteristics, in addition to their location. 
# For example, m can be used to specify an immigration rate and JM can be used to specify the metacommunity 
# size (see Hubbell 2001), where JM is the number of individual organisms that inhabit the metacommunity.
ESV_dataframe %>% as.data.frame %>% 
  ggplot(aes(Longitude, Latitude)) + geom_point(aes(colour=HabitatType), alpha=3/4) + 
  scale_color_manual(values=c("#88A550","#336B87")) +
  ggtitle("Plot Locations") + coord_equal() + theme_bw()

xy.cord<-data.frame(ESV_dataframe$Longitude, ESV_dataframe$Latitude); colnames(xy.cord)<-c("Longitude", "Latitude")

haka.landscape <- MCSim::fn.make.landscape(
  site.coords = xy.cord,
  m = 0.5,
  JM = 1000000)

 ####################################
 ### Kriging and Variogram Plots ###
 ###################################
# Taken mostly from https://rpubs.com/nabilabd/118172
install.packages('gstat');require('gstat')

# Will be working with dataframe created from phyloseq object title ESV_dataframe
str(ESV_dataframe)
theme_set(theme_classic())

#Can do fun things like plot where sample sites were 
ESV_dataframe %>% as.data.frame %>% 
  ggplot(aes(Longitude, Latitude)) + geom_point(aes(colour=HabitatType), alpha=3/4) + 
  scale_color_manual(values=c("#88A550","#336B87")) +
  ggtitle("Plot Locations") + coord_equal() + theme_bw()

### Need to convert dataframe to a spatial dataframe (SPDF)
#To convert it to a spatial dataframe, we must first specify which of the columns contain the coordinates of the data. 
#This is done by using R's formula notation as follows:
coordinates(ESV_dataframe) <- ~ xCoords + yCoords
class(ESV_dataframe)
str(ESV_dataframe)

      #############################################################
      ### Subset dataframe to get at plot-level spatial sorting  ##
      #############################################################
# Remnant plots
rem_all <- subset(ESV_dataframe, HabitatType == "Remnant Forest")

rem_RO1_all <- subset(rem_all, Plot == "RO1")
rem_RO1_koa <-subset(rem_RO1_all,Host=="Acacia koa")
rem_RO1_cop <-subset(rem_RO1_all,Host=="Corposma rhynchocarpa")
rem_RO1_chei <-subset(rem_RO1_all,Host=="Cheirodendron trigynum")
rem_RO1_grass <-subset(rem_RO1_all,Host=="Grass")
rem_RO1_ohia <- subset(rem_RO1_all,Host=="Metrosideros polymorpha")
rem_RO1_myr <- subset(rem_RO1_all,Host=="Myrsine lessertiana")
rem_RO1_rub <- subset(rem_RO1_all,Host=="Rubus hawaiiensis")

rem_RO2_all <- subset(rem_all, Plot == "RO2")
rem_RO2_koa <-subset(rem_RO2_all,Host=="Acacia koa")
rem_RO2_cop <-subset(rem_RO2_all,Host=="Corposma rhynchocarpa")
rem_RO2_chei <-subset(rem_RO2_all,Host=="Cheirodendron trigynum")
rem_RO2_grass <-subset(rem_RO2_all,Host=="Grass")
rem_RO2_ohia <- subset(rem_RO2_all,Host=="Metrosideros polymorpha")
rem_RO2_myr <- subset(rem_RO2_all,Host=="Myrsine lessertiana")
rem_RO2_rub <- subset(rem_RO2_all,Host=="Rubus hawaiiensis")

rem_RO3_all <- subset(rem_all, Plot == "RO3")
rem_RO3_koa <-subset(rem_RO3_all,Host=="Acacia koa")
rem_RO3_cop <-subset(rem_RO3_all,Host=="Corposma rhynchocarpa")
rem_RO3_chei <-subset(rem_RO3_all,Host=="Cheirodendron trigynum")
rem_RO3_grass <-subset(rem_RO3_all,Host=="Grass")
rem_RO3_ohia <- subset(rem_RO3_all,Host=="Metrosideros polymorpha")
rem_RO3_myr <- subset(rem_RO3_all,Host=="Myrsine lessertiana")
rem_RO3_rub <- subset(rem_RO3_all,Host=="Rubus hawaiiensis")

rem_RO4_all <- subset(rem_all, Plot == "RO4")
rem_RO4_koa <-subset(rem_RO4_all,Host=="Acacia koa")
rem_RO4_cop <-subset(rem_RO4_all,Host=="Corposma rhynchocarpa")
rem_RO4_chei <-subset(rem_RO4_all,Host=="Cheirodendron trigynum")
rem_RO4_grass <-subset(rem_RO4_all,Host=="Grass")
rem_RO4_ohia <- subset(rem_RO4_all,Host=="Metrosideros polymorpha")
rem_RO4_myr <- subset(rem_RO4_all,Host=="Myrsine lessertiana")
rem_RO4_rub <- subset(rem_RO4_all,Host=="Rubus hawaiiensis")

rem_RO5_all <- subset(rem_all, Plot == "RO5")
rem_RO5_koa <-subset(rem_RO5_all,Host=="Acacia koa")
rem_RO5_cop <-subset(rem_RO5_all,Host=="Corposma rhynchocarpa")
rem_RO5_chei <-subset(rem_RO5_all,Host=="Cheirodendron trigynum")
rem_RO5_grass <-subset(rem_RO5_all,Host=="Grass")
rem_RO5_ohia <- subset(rem_RO5_all,Host=="Metrosideros polymorpha")
rem_RO5_myr <- subset(rem_RO5_all,Host=="Myrsine lessertiana")
rem_RO5_rub <- subset(rem_RO5_all,Host=="Rubus hawaiiensis")

rem_RO6_all <- subset(rem_all, Plot == "RO6")
rem_RO6_koa <-subset(rem_RO6_all,Host=="Acacia koa")
rem_RO6_cop <-subset(rem_RO6_all,Host=="Corposma rhynchocarpa")
rem_RO6_chei <-subset(rem_RO6_all,Host=="Cheirodendron trigynum")
rem_RO6_grass <-subset(rem_RO6_all,Host=="Grass")
rem_RO6_ohia <- subset(rem_RO6_all,Host=="Metrosideros polymorpha")
rem_RO6_myr <- subset(rem_RO6_all,Host=="Myrsine lessertiana")
rem_RO6_rub <- subset(rem_RO6_all,Host=="Rubus hawaiiensis")

# Restored plots
rest_all <- subset(ESV_dataframe, HabitatType == "Restored Forest")

rest_AK1_all <- subset(rem_all, Plot == "AK1")
rest_AK1_koa <-subset(rest_AK1_all,Host=="Acacia koa")
rest_AK1_cop <-subset(rest_AK1_all,Host=="Corposma rhynchocarpa")
rest_AK1_chei <-subset(rest_AK1_all,Host=="Cheirodendron trigynum")
rest_AK1_grass <-subset(rest_AK1_all,Host=="Grass")
rest_AK1_ohia <- subset(rest_AK1_all,Host=="Metrosideros polymorpha")
rest_AK1_myr <- subset(rest_AK1_all,Host=="Myrsine lessertiana")
rest_AK1_rub <- subset(rest_AK1_all,Host=="Rubus hawaiiensis")

rest_AK2_all <- subset(rem_all, Plot == "AK2")
rest_AK2_koa <-subset(rest_AK2_all,Host=="Acacia koa")
rest_AK2_cop <-subset(rest_AK2_all,Host=="Corposma rhynchocarpa")
rest_AK2_chei <-subset(rest_AK2_all,Host=="Cheirodendron trigynum")
rest_AK2_grass <-subset(rest_AK2_all,Host=="Grass")
rest_AK2_ohia <- subset(rest_AK2_all,Host=="Metrosideros polymorpha")
rest_AK2_myr <- subset(rest_AK2_all,Host=="Myrsine lessertiana")
rest_AK2_rub <- subset(rest_AK2_all,Host=="Rubus hawaiiensis")

rest_AK3_all <- subset(rem_all, Plot == "AK3")
rest_AK3_koa <-subset(rest_AK3_all,Host=="Acacia koa")
rest_AK3_cop <-subset(rest_AK3_all,Host=="Corposma rhynchocarpa")
rest_AK3_chei <-subset(rest_AK3_all,Host=="Cheirodendron trigynum")
rest_AK3_grass <-subset(rest_AK3_all,Host=="Grass")
rest_AK3_ohia <- subset(rest_AK3_all,Host=="Metrosideros polymorpha")
rest_AK3_myr <- subset(rest_AK3_all,Host=="Myrsine lessertiana")
rest_AK3_rub <- subset(rest_AK3_all,Host=="Rubus hawaiiensis")

rest_AK4_all <- subset(rem_all, Plot == "AK4")
rest_AK4_koa <-subset(rest_AK4_all,Host=="Acacia koa")
rest_AK4_cop <-subset(rest_AK4_all,Host=="Corposma rhynchocarpa")
rest_AK4_chei <-subset(rest_AK4_all,Host=="Cheirodendron trigynum")
rest_AK4_grass <-subset(rest_AK4_all,Host=="Grass")
rest_AK4_ohia <- subset(rest_AK4_all,Host=="Metrosideros polymorpha")
rest_AK4_myr <- subset(rest_AK4_all,Host=="Myrsine lessertiana")
rest_AK4_rub <- subset(rest_AK4_all,Host=="Rubus hawaiiensis")

rest_AK5_all <- subset(rem_all, Plot == "AK5")
rest_AK5_koa <-subset(rest_AK5_all,Host=="Acacia koa")
rest_AK5_cop <-subset(rest_AK5_all,Host=="Corposma rhynchocarpa")
rest_AK5_chei <-subset(rest_AK5_all,Host=="Cheirodendron trigynum")
rest_AK5_grass <-subset(rest_AK5_all,Host=="Grass")
rest_AK5_ohia <- subset(rest_AK5_all,Host=="Metrosideros polymorpha")
rest_AK5_myr <- subset(rest_AK5_all,Host=="Myrsine lessertiana")
rest_AK5_rub <- subset(rest_AK5_all,Host=="Rubus hawaiiensis")

rest_AK6_all <- subset(rem_all, Plot == "AK6")
rest_AK6_koa <-subset(rest_AK6_all,Host=="Acacia koa")
rest_AK6_cop <-subset(rest_AK6_all,Host=="Corposma rhynchocarpa")
rest_AK6_chei <-subset(rest_AK6_all,Host=="Cheirodendron trigynum")
rest_AK6_grass <-subset(rest_AK6_all,Host=="Grass")
rest_AK6_ohia <- subset(rest_AK6_all,Host=="Metrosideros polymorpha")
rest_AK6_myr <- subset(rest_AK6_all,Host=="Myrsine lessertiana")
rest_AK6_rub <- subset(rest_AK6_all,Host=="Rubus hawaiiensis")

       ###########################
       #### Fitting variograms ###
       ###########################
        ################
        ## Plot-level ##
        ################
# To perform kriging, you must first have a variogram model, from which the data can be interpolated. 
# There are a couple steps involved:
# Calculate the sample variogram for each plot using the variogram function.

# RO1
rem_RO1_all.vgm <- variogram(Abundance~1, rem_RO1_all) 
rem_RO1_all.fit <- fit.variogram(rem_RO1_all.vgm, model=vgm(0.00015, "Sph", 0.0001, 1))

rem_RO1_koa.vgm <- variogram(Abundance~1, rem_RO1_koa)
rem_RO1_chei.vgm <- variogram(Abundance~1, rem_RO1_chei)
rem_RO1_myr.vgm <- variogram(Abundance~1, rem_RO1_myr)
rem_RO1_ohia.vgm <- variogram(Abundance~1, rem_RO1_ohia)
rem_RO1_grass.vgm <- variogram(Abundance~1, rem_RO1_grass)
rem_RO1_rub.vgm <- variogram(Abundance~1, rem_RO1_rub)

# Combine all vgm from the same plot into a single dataframe
RO1_host_labs<-rep(c("All","Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                     "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                   times=c(15,5,3,3,6,10,4))
rem_RO1_vgm <- rbind(rem_RO1_all.vgm,rem_RO1_ohia.vgm,rem_RO1_koa.vgm,rem_RO1_myr.vgm,
                     rem_RO1_chei.vgm,rem_RO1_rub.vgm,rem_RO1_grass.vgm)
rem_RO1_vgm <- cbind(rem_RO1_vgm,RO1_host_labs)
colnames(rem_RO1_vgm) <- c("np","dist","gamma","dir.hor","dir.ver","id","Host")
rem_RO1_vgm$Host<- factor(rem_RO1_vgm$Host,levels=c("All","Metrosideros polymorpha", "Acacia koa", 
                                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                          labels=c("All","Metrosideros polymorpha", "Acacia koa",  
                                   "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Plot variogram results
RO1_vgm_plot <- ggplot(rem_RO1_vgm,aes(x=dist,y=gamma)) + geom_point(alpha=0.75,size=5,aes(colour=Host)) +
  scale_colour_manual(values=c("black","#D73027","#FC8D59","#008000","#E0F3F8","#008B8B","#4575B4")) +
  ylab("Semivariance") + xlab("Distance")
plot(RO1_vgm_plot)

ggplot(data.frame(rem_RO1_vgm),aes(x=map.dx,y=map.dy,fill=map.z))+ 
  geom_raster() + scale_fill_gradient() +
  facet_grid(~Host)

# RO2
rem_RO2_all.vgm <- variogram(Abundance~1, rem_RO2_all) 
head(rem_RO2_all.vgm)
rem_RO2_koa.vgm <- variogram(Abundance~1, rem_RO2_koa)
rem_RO2_chei.vgm <- variogram(Abundance~1, rem_RO2_chei)
rem_RO2_myr.vgm <- variogram(Abundance~1, rem_RO2_myr)
rem_RO2_ohia.vgm <- variogram(Abundance~1, rem_RO2_ohia)
rem_RO2_grass.vgm <- variogram(Abundance~1, rem_RO2_grass)
rem_RO2_rub.vgm <- variogram(Abundance~1, rem_RO2_rub)

# Combine all vgm from the same plot into a single dataframe
RO2_host_labs<-rep(c("All","Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                     "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                   times=c(12,5,9,5,8,1,5))
rem_RO2_vgm <- rbind(rem_RO2_all.vgm,rem_RO2_ohia.vgm,rem_RO2_koa.vgm,rem_RO2_myr.vgm,
                     rem_RO2_chei.vgm,rem_RO2_rub.vgm,rem_RO2_grass.vgm)
rem_RO2_vgm <- cbind(rem_RO2_vgm,RO2_host_labs)
colnames(rem_RO2_vgm) <- c("np","dist","gamma","dir.hor","dir.ver","id","Host")
rem_RO2_vgm$Host<- factor(rem_RO2_vgm$Host,levels=c("All","Metrosideros polymorpha", "Acacia koa", 
                                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                          labels=c("All","Metrosideros polymorpha", "Acacia koa",  
                                   "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Plot variogram results
RO2_vgm_plot <- ggplot(rem_RO2_vgm,aes(x=dist,y=gamma)) + geom_point(alpha=0.75,size=5,aes(colour=Host)) +
  scale_colour_manual(values=c("black","#D73027","#FC8D59","#008000","#E0F3F8","#008B8B","#4575B4")) +
  ylab("Semivariance") + xlab("Distance")
plot(RO2_vgm_plot)

# RO3
rem_RO3_all.vgm <- variogram(Abundance~1, rem_RO3_all) 
rem_RO3_koa.vgm <- variogram(Abundance~1, rem_RO3_koa)
rem_RO3_chei.vgm <- variogram(Abundance~1, rem_RO3_chei)
rem_RO3_myr.vgm <- variogram(Abundance~1, rem_RO3_myr)
rem_RO3_ohia.vgm <- variogram(Abundance~1, rem_RO3_ohia)
rem_RO3_grass.vgm <- variogram(Abundance~1, rem_RO3_grass)
rem_RO3_rub.vgm <- variogram(Abundance~1, rem_RO3_rub)

# Combine all vgm from the same plot into a single dataframe
RO3_host_labs<-rep(c("All","Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                     "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                   times=c(15,2,4,1,5,6,4))
rem_RO3_vgm <- rbind(rem_RO3_all.vgm,rem_RO3_ohia.vgm,rem_RO3_koa.vgm,rem_RO3_myr.vgm,
                     rem_RO3_chei.vgm,rem_RO3_rub.vgm,rem_RO3_grass.vgm)
rem_RO3_vgm <- cbind(rem_RO3_vgm,RO3_host_labs)
colnames(rem_RO3_vgm) <- c("np","dist","gamma","dir.hor","dir.ver","id","Host")
rem_RO3_vgm$Host<- factor(rem_RO3_vgm$Host,levels=c("All","Metrosideros polymorpha", "Acacia koa", 
                                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                          labels=c("All","Metrosideros polymorpha", "Acacia koa",  
                                   "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Plot variogram results
RO3_vgm_plot <- ggplot(rem_RO3_vgm,aes(x=dist,y=gamma)) + geom_point(alpha=0.75,size=5,aes(colour=Host)) +
  scale_colour_manual(values=c("black","#D73027","#FC8D59","#008000","#E0F3F8","#008B8B","#4575B4")) +
  ylab("Semivariance") + xlab("Distance")
plot(RO3_vgm_plot)

# RO4
rem_RO4_all.vgm <- variogram(Abundance~1, rem_RO4_all) 
head(rem_RO4_all.vgm)
rem_RO4_koa.vgm <- variogram(Abundance~1, rem_RO4_koa)
rem_RO4_chei.vgm <- variogram(Abundance~1, rem_RO4_chei)
rem_RO4_myr.vgm <- variogram(Abundance~1, rem_RO4_myr)
rem_RO4_ohia.vgm <- variogram(Abundance~1, rem_RO4_ohia)
rem_RO4_grass.vgm <- variogram(Abundance~1, rem_RO4_grass)
rem_RO4_rub.vgm <- variogram(Abundance~1, rem_RO4_rub)

# Combine all vgm from the same plot into a single dataframe
RO4_host_labs<-rep(c("All","Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                     "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                   times=c(15,6,1,1,3,5,6))
rem_RO4_vgm <- rbind(rem_RO4_all.vgm,rem_RO4_ohia.vgm,rem_RO4_koa.vgm,rem_RO4_myr.vgm,
                     rem_RO4_chei.vgm,rem_RO4_rub.vgm,rem_RO4_grass.vgm)
rem_RO4_vgm <- cbind(rem_RO4_vgm,RO4_host_labs)
colnames(rem_RO4_vgm) <- c("np","dist","gamma","dir.hor","dir.ver","id","Host")
rem_RO4_vgm$Host<- factor(rem_RO4_vgm$Host,levels=c("All","Metrosideros polymorpha", "Acacia koa", 
                                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                          labels=c("All","Metrosideros polymorpha", "Acacia koa",  
                                   "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Plot variogram results
RO4_vgm_plot <- ggplot(rem_RO4_vgm,aes(x=dist,y=gamma)) + geom_point(alpha=0.75,size=5,aes(colour=Host)) +
  scale_colour_manual(values=c("black","#D73027","#FC8D59","#008000","#E0F3F8","#008B8B","#4575B4")) +
  ylab("Semivariance") + xlab("Distance")
plot(RO4_vgm_plot)

# RO5
rem_RO5_all.vgm <- variogram(Abundance~1, rem_RO5_all) 
head(rem_RO5_all.vgm)
rem_RO5_koa.vgm <- variogram(Abundance~1, rem_RO5_koa)
rem_RO5_chei.vgm <- variogram(Abundance~1, rem_RO5_chei)
rem_RO5_myr.vgm <- variogram(Abundance~1, rem_RO5_myr)
rem_RO5_ohia.vgm <- variogram(Abundance~1, rem_RO5_ohia)
rem_RO5_grass.vgm <- variogram(Abundance~1, rem_RO5_grass)
rem_RO5_rub.vgm <- variogram(Abundance~1, rem_RO5_rub)

# Combine all vgm from the same plot into a single dataframe
RO5_host_labs<-rep(c("All","Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                     "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                   times=c(15,1,1,1,4,4,5))
rem_RO5_vgm <- rbind(rem_RO5_all.vgm,rem_RO5_ohia.vgm,rem_RO5_koa.vgm,rem_RO5_myr.vgm,
                     rem_RO5_chei.vgm,rem_RO5_rub.vgm,rem_RO5_grass.vgm)
rem_RO5_vgm <- cbind(rem_RO5_vgm,RO5_host_labs)
colnames(rem_RO5_vgm) <- c("np","dist","gamma","dir.hor","dir.ver","id","Host")
rem_RO5_vgm$Host<- factor(rem_RO5_vgm$Host,levels=c("All","Metrosideros polymorpha", "Acacia koa", 
                                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                          labels=c("All","Metrosideros polymorpha", "Acacia koa",  
                                   "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Plot variogram results
RO5_vgm_plot <- ggplot(rem_RO5_vgm,aes(x=dist,y=gamma)) + geom_point(alpha=0.75,size=5,aes(colour=Host)) +
  scale_colour_manual(values=c("black","#D73027","#FC8D59","#008000","#E0F3F8","#008B8B","#4575B4")) +
  ylab("Semivariance") + xlab("Distance")
plot(RO5_vgm_plot)

# RO6
rem_RO6_all.vgm <- variogram(Abundance~1, rem_RO6_all) 
head(rem_RO6_all.vgm)
rem_RO6_koa.vgm <- variogram(Abundance~1, rem_RO6_koa)
rem_RO6_chei.vgm <- variogram(Abundance~1, rem_RO6_chei)
rem_RO6_myr.vgm <- variogram(Abundance~1, rem_RO6_myr)
rem_RO6_ohia.vgm <- variogram(Abundance~1, rem_RO6_ohia)
rem_RO6_grass.vgm <- variogram(Abundance~1, rem_RO6_grass)
rem_RO6_rub.vgm <- variogram(Abundance~1, rem_RO6_rub)

# Combine all vgm from the same plot into a single dataframe
RO6_host_labs<-rep(c("All","Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                     "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                   times=c(15,1,2,6,5,1,2))
rem_RO6_vgm <- rbind(rem_RO6_all.vgm,rem_RO6_ohia.vgm,rem_RO6_koa.vgm,rem_RO6_myr.vgm,
                     rem_RO6_chei.vgm,rem_RO6_rub.vgm,rem_RO6_grass.vgm)
rem_RO6_vgm <- cbind(rem_RO6_vgm,RO6_host_labs)
colnames(rem_RO6_vgm) <- c("np","dist","gamma","dir.hor","dir.ver","id","Host")
rem_RO6_vgm$Host<- factor(rem_RO6_vgm$Host,levels=c("All","Metrosideros polymorpha", "Acacia koa", 
                                                    "Myrsine lessertiana", 
                                                    "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                          labels=c("All","Metrosideros polymorpha", "Acacia koa", 
                                   "Myrsine lessertiana", 
                                   "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))


#Plot variogram results
RO6_vgm_plot <- ggplot(rem_RO6_vgm,aes(x=dist,y=gamma)) + geom_point(alpha=0.75,size=5,aes(colour=Host)) +
  scale_colour_manual(values=c("black","#D73027","#FC8D59","#008000","#E0F3F8","#008B8B","#4575B4")) +
  ylab("Semivariance") + xlab("Distance")
plot(RO6_vgm_plot)
         ###################
         ## Habitat-level ##
         ###################
# Remnant forest
remnant_vgm_df <- rbind(rem_RO1_vgm,rem_RO2_vgm,rem_RO3_vgm,rem_RO4_vgm,rem_RO5_vgm,rem_RO6_vgm)
remnant_forest_vgm_plot <- ggplot(remnant_vgm_df,aes(x=dist,y=gamma)) + 
  geom_point(alpha=0.75,size=5,aes(colour=Host)) +
  ylab("Semivariance") + xlab("Distance") +
  scale_colour_manual(values=c("black","#D73027","#FC8D59","#008000","#E0F3F8","#008B8B","#4575B4")) 

plot(remnant_forest_vgm_plot)


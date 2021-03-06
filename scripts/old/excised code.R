# removed code


hakalau_map_zoom <-get_map(location=c(-155.320,19.83),zoom=14,maptype="satellite",color="color")
haka_map <- ggmap(hakalau_map_zoom) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text=element_text(colour="black",size=8)) +
  theme(axis.title=element_text(colour="black",size=8))
plot(haka_map)

groups.sp<-c("coral", "seagreen", "mediumpurple", "goldenrod", "dodgerblue", "gray50", "navajowhite2")

##########################
####   Haka hab types ####
##########################

haka_hab_types <-read.csv("data/haka_big_transect_data.csv", header=TRUE, row.names=1)
haka_hab_types$NewHabType<- factor(haka_hab_types$NewHabType,levels=c("Open pasture","Koa + Grass","Koa + Understory",
                                                                      "Restored ohia","Remnant koa forest",
                                                                      "Remnant ohia forest"),
                                   labels=c("Open pasture","Koa + Grass","Koa + Understory",
                                            "Restored ohia forest","Remnant koa dominated forest",
                                            "Remnant ohia dominated forest"))

hab_types_plot <- haka_map +
  geom_point(data=haka_hab_types,aes(x=Longitude,y=Latitude,fill=NewHabType),pch=21,stroke=0.3,colour="white",size=2) +
  scale_fill_manual(values=c("#F7AF51","#97D8E5","#1C84B5","#336B87","#88A550","dark green")) +
  theme(legend.text=element_text(size=6),legend.title = element_text(size=8)) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text=element_text(colour="black",size=8)) +
  theme(axis.title=element_text(colour="black",size=8)) +
  
  plot(hab_types_plot)
ggsave("figures/Haka_hab_types.tiff",width= 6,height=5,plot=hab_types_plot)

hab_types_plot_no_legend <- haka_map +
  geom_point(data=haka_hab_types,aes(x=Longitude,y=Latitude,fill=NewHabType),pch=21,stroke=0.3,colour="white",size=2) +
  scale_fill_manual(values=c("#F7AF51","#97D8E5","#1C84B5","#336B87","#88A550","dark green")) +
  theme(legend.position="none") +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text=element_text(colour="black",size=8)) +
  theme(axis.title=element_text(colour="black",size=8))

plot(hab_types_plot_no_legend)
ggsave("figures/Haka_hab_types_no_legend.tiff",width= 6,height=5,plot=hab_types_plot_no_legend)


############
############



spec_rich_host = plot_richness(haka_VT_soil_physeq, x ="Host", measures="Observed", color="Host")  +
  geom_boxplot(col="black", aes(fill=Host), alpha=0.8 , lwd=0.5, outlier.colour = "gray50") +
  scale_color_manual(name= "Host species", values=c("#D73027","#FC8D59","goldenrod","#008000","darkslategray3","#008B8B","#4575B4")) +
  scale_fill_manual(name= "Host species", values=c("#D73027","#FC8D59","goldenrod","#008000","darkslategray3","#008B8B","#4575B4"))  + 
  xlab("Habitat Type") + ylab("AM fungal richness") + 
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_text(colour="black",size=8)) +
  scale_y_continuous(breaks=seq(0,45,by=10),limits=c(0,45)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  facet_wrap(~HabitatType) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(strip.text.x=element_text(size=10,face="bold"),strip.background = element_rect(fill="gray90")) +
  theme(legend.text = element_text(face="italic")) 

plot(spec_rich_host)
ggsave("figures/host_spec_rich.tiff", plot = spec_rich_host, width=7,height=5)





###### NMDS plot by Habitat
ordiplot(NMDS, type="n", main=substitute(paste("")), cex.main=1, display="sites", xlim=c(-0.25, 0.8), cex.lab=0.8, cex.axis=0.8)
abline(h = 0, lty = "dotted")
abline(v = 0, lty = "dotted")
points(NMDS, "sites", cex=0.8, pch=16, col=Hab.col[Habitats])
ordiellipse(NMDS, groups=Habitats, kind="sd", draw="polygon", conf=0.95, alpha=30, col=groups.hab, border=groups.hab)
legend("topright", legend=levels(Habitats), cex=1, pch=16, col=groups.hab, pt.cex=1, bty="n")

dev.copy(pdf, "figures/habitat.NMDS.pdf", height=5.5, width=7)
dev.off() 


##### ##### ##### ##### ##### ##### #####

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


###### 
###################





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




#### T-tests
# Welch Tests
centrality_between_hab <- t.test(ig.RO_between,ig.AK_between, paired=FALSE, var.equal=FALSE)
centrality_between_hab # signific, but barely. Diff 'betweenness' by habitat

connect_degree_hab <- t.test(ig.RO_degree,ig.AK_degree, paired=FALSE, var.equal=FALSE)
connect_degree_hab # NS similar 'degree nodes' by habitat



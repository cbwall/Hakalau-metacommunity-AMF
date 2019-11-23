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


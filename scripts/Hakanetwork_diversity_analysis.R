install.packages("BiocManager")
BiocManager::install(c("multtest","phyloseq","rhdf5","ggplot2","colorspace","stringi"))

setwd("C:/Users/egan_/Documents/R work/Hakalau - network")
library(ggplot2); library(phyloseq);library(vegan);library(multcomp);library(dplyr);library(grid);library(scales)
require(gridExtra); library(emmeans);library(multcompView); library(ggpubr); library(Rmisc); library(RColorBrewer)
library(purrr); library(RVAideMemoire)

                                                    ###################################
                                                    ####Import files and pre-process###
                                                    ###################################
#Import files generated in QIIME
otu <- as.matrix(read.csv("haka_feb_large_species_table.csv", header = TRUE,row.names = 1))
haka_otu <- t(otu)
taxmat <- as.matrix(read.csv("haka_feb_large_tax.csv", header = TRUE,row.names = 1))
sample <-read.csv("haka_metadata.csv", header=TRUE, row.names=1)
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
                                    Na.cation= sample$Na.cation, stringsAsFactors = FALSE))
row.names(sampledata) <- row.names(sample)
#Change each file to the phyloseq format
OTU = otu_table(haka_otu, taxa_are_rows = FALSE)
physeq = phyloseq(OTU)
TAX = tax_table(taxmat)
haka_physeq = merge_phyloseq(OTU, sampledata, TAX)
#Create new physeq object working at only the species level taxonomy. Collapses all ESVs identified
#as the same species, and subset to only include root samples
vt_physeq <- tax_glom(haka_physeq,"Species")
vt_physeq <- subset_samples(vt_physeq,SampleType=="roots")
vt_rel_abund <- transform_sample_counts(vt_physeq,function(x)x/sum(x))

#Melt phyloseq object to make a dataframe for ggplot and bipartite
vt.melt<-psmelt(vt_rel_abund)
haka_root_df<-vt.melt[which(vt.melt$SampleType=='roots'),]
#Subset sample metadata and otu table to only include rows with root samples
root_metadata <-as(sample_data(vt_rel_abund),"data.frame")
root_otu <- as(otu_table(vt_rel_abund),"matrix")


  
                                                      ################################
                                                ####      Family Relative Abundance     ###
                                                      ################################
#Family Relative Abundance
haka_fam<-tax_glom(vt_rel_abund,taxrank="Family")
ps.melt <- psmelt(haka_fam)

log_hab_type_family_plot = ggplot(ps.melt,aes(x=HabitatType,y=Abundance,fill=HabitatType)) +
  geom_boxplot(size=0.75, outlier.alpha = 0.5) +
  ylab("Relative Abundance") + xlab("Habitat Type") +
  scale_fill_manual(values=c("#88A550","#336B87")) +
  theme(text=element_text(colour="black",size=15)) + 
  scale_y_continuous(trans='log10') +
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title=element_text(size=20,face='bold')) +
  facet_wrap(~Family) +
  theme(strip.text.x=element_text(size=20,face="bold"),strip.background = element_rect(fill="white")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5)

plot(log_hab_type_family_plot)
ggsave("Fig_1_haka_fam_log_rel_abund.tiff", 
       plot = log_hab_type_family_plot, width= 20,height=15,limitsize=FALSE)

## Stats ##
#Acaulosporaceae
Acaulosporaceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                    Family=="Acaulosporaceae")$Abundance,
                            subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                     Family=="Acaulosporaceae")$Abundance,
                            paired=FALSE,var.equal=FALSE)
Acaulosporaceae_root_welch
Acaulosporaceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                              Family=="Acaulosporaceae")$Abundance,
                                     subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                              Family=="Acaulosporaceae")$Abundance,
                                     paired=FALSE,var.equal=FALSE)
Acaulosporaceae_soil_welch
#Ambisporaceae
Ambisporaceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                              Family=="Ambisporaceae")$Abundance,
                                     subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                              Family=="Ambisporaceae")$Abundance,
                                     paired=FALSE,var.equal=FALSE)
Ambisporaceae_root_welch
Ambisporaceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                              Family=="Ambisporaceae")$Abundance,
                                     subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                              Family=="Ambisporaceae")$Abundance,
                                     paired=FALSE,var.equal=FALSE)
Ambisporaceae_soil_welch
#Archaeosporaceae
Archaeosporaceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                            Family=="Archaeosporaceae")$Abundance,
                                   subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                            Family=="Archaeosporaceae")$Abundance,
                                   paired=FALSE,var.equal=FALSE)
Archaeosporaceae_root_welch
Archaeosporaceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                            Family=="Archaeosporaceae")$Abundance,
                                   subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                            Family=="Archaeosporaceae")$Abundance,
                                   paired=FALSE,var.equal=FALSE)
Archaeosporaceae_soil_welch
#Claroideoglomeraceae
Claroideoglomeraceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                               Family=="Claroideoglomeraceae")$Abundance,
                                      subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                               Family=="Claroideoglomeraceae")$Abundance,
                                      paired=FALSE,var.equal=FALSE)
Claroideoglomeraceae_root_welch
Claroideoglomeraceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                               Family=="Claroideoglomeraceae")$Abundance,
                                      subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                               Family=="Claroideoglomeraceae")$Abundance,
                                      paired=FALSE,var.equal=FALSE)
Claroideoglomeraceae_soil_welch
#Diversisporaceae
Diversisporaceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                                   Family=="Diversisporaceae")$Abundance,
                                          subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                                   Family=="Diversisporaceae")$Abundance,
                                          paired=FALSE,var.equal=FALSE)
Diversisporaceae_root_welch
Diversisporaceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                                   Family=="Diversisporaceae")$Abundance,
                                          subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                                   Family=="Diversisporaceae")$Abundance,
                                          paired=FALSE,var.equal=FALSE)
Diversisporaceae_soil_welch
#Geosiphonaceae
Geosiphonaceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                               Family=="Geosiphonaceae")$Abundance,
                                      subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                               Family=="Geosiphonaceae")$Abundance,
                                      paired=FALSE,var.equal=FALSE)
Geosiphonaceae_root_welch
Geosiphonaceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                               Family=="Geosiphonaceae")$Abundance,
                                      subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                               Family=="Geosiphonaceae")$Abundance,
                                      paired=FALSE,var.equal=FALSE)
Geosiphonaceae_soil_welch
#Gigasporaceae
Gigasporaceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                             Family=="Gigasporaceae")$Abundance,
                                    subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                             Family=="Gigasporaceae")$Abundance,
                                    paired=FALSE,var.equal=FALSE)
Gigasporaceae_root_welch
Gigasporaceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                             Family=="Gigasporaceae")$Abundance,
                                    subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                             Family=="Gigasporaceae")$Abundance,
                                    paired=FALSE,var.equal=FALSE)
Gigasporaceae_soil_welch
#Glomeraceae
Glomeraceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                            Family=="Glomeraceae")$Abundance,
                                   subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                            Family=="Glomeraceae")$Abundance,
                                   paired=FALSE,var.equal=FALSE)
Glomeraceae_root_welch
Glomeraceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                            Family=="Glomeraceae")$Abundance,
                                   subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                            Family=="Glomeraceae")$Abundance,
                                   paired=FALSE,var.equal=FALSE)
Glomeraceae_soil_welch
#Paraglomeraceae
Paraglomeraceae_root_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="roots" & 
                                          Family=="Paraglomeraceae")$Abundance,
                                 subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="roots" &
                                          Family=="Paraglomeraceae")$Abundance,
                                 paired=FALSE,var.equal=FALSE)
Paraglomeraceae_root_welch
Paraglomeraceae_soil_welch <- t.test(subset(ps.melt,HabitatType == "Afforested koa" & SampleType=="soil" & 
                                          Family=="Paraglomeraceae")$Abundance,
                                 subset(ps.melt,HabitatType == "Remnant ohia" & SampleType=="soil" &
                                          Family=="Paraglomeraceae")$Abundance,
                                 paired=FALSE,var.equal=FALSE)
Paraglomeraceae_soil_welch
                                                        #######################
                                                  ####     SPECIES RICHNESS      ###
                                                        #######################
###Phyloseq###
#Richness summary stats
library(doBy)
spec_rich <- estimate_richness(vt_physeq,measures="Observed")
all.equal(rownames(spec_rich),rownames(root_metadata))
row.names(spec_rich) <- row.names(root_metadata)
spec_rich <- cbind(root_metadata,spec_rich)
rich_sum_stats <- summarySE(spec_rich,measurevar="Observed",groupvars=c("HabitatType","Host"))
rich_sum_stats

hab_sum_stats <- summarySE(spec_rich,measurevar="Observed",groupvars=c("HabitatType"))
hab_sum_stats


#Plots and analyses
spec_rich_host = plot_richness(vt_physeq, x ="HabitatType", measures="Observed")  +
  geom_boxplot(aes(fill=Host),size=0.75, outlier.alpha = 0.5) +  
  theme(text=element_text(colour="black",size=15)) + 
  xlab("Habitat Type") + ylab("AM fungal richness") + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15,vjust=1)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  scale_y_continuous(breaks=seq(0,45,by=10),limits=c(0,45)) +
  scale_fill_manual(values=c("#D73027","#FC8D59","#FEE090","#008000","#E0F3F8","#008B8B","#4575B4")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(strip.text.x=element_text(size=20,face="bold"),strip.background = element_rect(fill="white")) +
  theme(legend.text = element_text(face="italic")) 
  
plot(spec_rich_host)
ggsave("Fig_1_host_spec_rich.tiff", 
       plot = spec_rich_host,width=7,height=6)

#GLM with poisson error distribution to test for differences in Species richness in roots and soil 
#Large GLM 
host_glm1 <- glm(spec_rich$Observed ~ Host*HabitatType, data=root_metadata, family = poisson(link="log")) 
par(mfrow = c(2, 2))
plot(host_glm1)
chi.anova<-anova(host_glm1,test="Chisq")
chi.anova
#Create a unique level for every combination of site and treatment to do a post-hoc test on
root_metadata$Hab.by.Host<- interaction(root_metadata$HabitatType,root_metadata$Host)
host_glm2 <- glm(spec_rich$Observed ~ Hab.by.Host, data=root_metadata, family = poisson(link="log"))
#Calculate EMM on interactions
host_spec_rich_emmeans <- emmeans(host_glm2, specs="Hab.by.Host")
host_spec_rich_emmeans
host_spec_rich_posthoc.pairs = pairs(host_spec_rich_emmeans)
host_spec_rich_posthoc.pairs
#Create letters for interaction differences
host_spec_rich_mc_letters<-cld(host_spec_rich_emmeans,Letters="abcdefg")
host_spec_rich_mc_letters

##By habitat type alone
spec_rich_hab = plot_richness(vt_physeq, x ="Plot", measures="Observed") +
  geom_boxplot(aes(fill=HabitatType),size=0.75, outlier.alpha = 0.5) +  
  geom_jitter(shape=21,aes(fill=HabitatType)) +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  ylim(0,45) + 
  theme(text=element_text(colour="black",size=15)) + 
  ylab("AM fungal richness") + xlab("") +
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  facet_wrap(~SampleType) +
  theme(strip.text.x=element_text(size=20,face="bold"),strip.background = element_rect(fill="white")) 

plot(spec_rich_hab)
ggsave("species_richness_by_plot.tiff", 
       plot = spec_rich_hab, width=10,height=6)




################################
#    Nestedness Temperature    #
################################
library(analytics)
remnant.ps <- subset_samples(vt_rel_abund,HabitatType=="Remnant Forest")
remnant_metadata <- as(sample_data(remnant.ps),"data.frame")
remnant_otu <- as(otu_table(remnant.ps),"matrix")
rownames(remnant_otu)<- remnant_metadata$Host
binned_rem_otu <- rowmean(remnant_otu)

restored.ps <- subset_samples(vt_rel_abund,HabitatType=="Restored Forest")
restored_metadata <-as(sample_data(restored.ps),"data.frame")
restored_otu<-as(otu_table(restored.ps),"matrix")
rownames(restored_otu)<-restored_metadata$Host
binned_rest_otu <- rowmean(restored_otu)

remnant_nest_temp <- nestedtemp(binned_rem_otu)
par(mar=c(2.1,8,4.1,2.1))
plot(remnant_nest_temp,names=TRUE, kind="incidence")
set.seed(96822)
remnant_oecosimu <- oecosimu(binned_rem_otu,nestfun=nestedtemp,"r0", nsimul=9999)
saveRDS(remnant_oecosimu,"remnant_oecosimu.RData")
remnant_oecosimu <- readRDS("remnant_oecosimu.RData")
remnant_oecosimu

restored_nest_temp<- nestedtemp(binned_rest_otu)
plot(restored_nest_temp,names=TRUE,kind="incidence")
set.seed(96822)
restored_oecosim <- oecosimu(binned_rest_otu,nestfun=nestedtemp,"r0", nsimul=9999)
saveRDS(restored_oecosim,"restored_oecosimu.RData")
restored_oecosimu <- readRDS("restored_oecosimu.RData")
restored_oecosimu



                                                          #######################
                                                    ####       Beta Diversity       ###
                                                          #######################
###Phyloseq###
#Bray plot
bc_dist = as.matrix((vegdist(root_otu, "bray")))
NMDS = metaMDS(bc_dist)
NMDS1=NMDS$points[,1]
NMDS2=NMDS$points[,2]
NMDS.plot=data.frame(NMDS1=NMDS1,NMDS2=NMDS2, Host=root_metadata$Host, 
                     HabitatType=root_metadata$HabitatType,Plot=root_metadata$Plot)

set.seed(605)
NMDS_hab_plot = ggplot(NMDS.plot, aes(x=NMDS1,y=NMDS2, col=HabitatType, shape=HabitatType)) +
  geom_point(size=3, stroke=3) +
  scale_color_manual(values=c("#88A550","#336B87")) +
  scale_shape_manual(values=c(21,23)) +
  stat_ellipse(size=2) +
  theme(axis.text.x=element_text(colour="black",size=17)) +
  theme(axis.text.y=element_text(colour="black",size=17)) +
  theme(axis.title=element_text(colour="black",size=20)) +
  theme(panel.border = element_rect(fill=NA),axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(strip.text.x=element_text(size=20,face="bold"),strip.background = element_rect(fill="white")) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title=element_text(size=20,face='bold')) 

plot(NMDS_hab_plot)
ggsave("Extra_NMDS_hab_plot_no_env_data.tiff",width= 17,height=15,limitsize=FALSE,plot=NMDS_hab_plot)

set.seed(605)
NMDS_host_plot = ggplot(NMDS.plot, aes(x=NMDS1,y=NMDS2,shape=HabitatType)) +
  geom_point(size=4, stroke=3, aes(fill=HabitatType,colour=Host)) +
  scale_fill_manual(values=c("black","white")) +
  scale_shape_manual(values=c(21,23)) +
  stat_ellipse(size=2,alpha=0.75,aes(color=Host,linetype=HabitatType),type='t') +
  scale_linetype_manual(values=c(1,6)) +
  scale_color_manual(values=c("#D73027","#FC8D59","#FEE090","#008000","#E0F3F8","#008B8B","#4575B4")) +
  theme(axis.text.x=element_text(colour="black",size=17)) +
  theme(axis.text.y=element_text(colour="black",size=17)) +
  theme(axis.title=element_text(colour="black",size=20)) +
  theme(panel.border = element_rect(fill=NA),axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(strip.text.x=element_text(size=20,face="bold"),strip.background = element_rect(fill="white")) +
  theme(legend.text=element_text(size=20,face='italic')) +
  theme(legend.title=element_text(size=20,face='bold')) 

plot(NMDS_host_plot)
ggsave("Fig_S6_NMDS_host_plot.tiff",width= 20,height=15,limitsize=FALSE,plot=NMDS_host_plot)

#Examination of environmental data
environmental_data <- read.csv("haka_chemistry.csv", header=TRUE, row.names=1)
all.equal(rownames(root_otu), rownames(environmental_data))
rows_to_keep<-rownames(root_otu)
environmental_data<-environmental_data[rows_to_keep,]
all.equal(rownames(root_otu), rownames(environmental_data))

#Fit environmental vectors to ordination to see which environmental variables are correlated with the ordination
fit <- envfit(NMDS,environmental_data,na.rm=TRUE)
arrow<-data.frame(fit$vectors$arrows,R=fit$vectors$r, P=fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p <- dplyr::filter(arrow, P <=0.05)
label_map <- aes(x=1.3*(NMDS1/2),y=1.3*(NMDS2/2),shape=NULL,color=NULL,label=FG)

arrow.p$FG<- as.factor(arrow.p$FG)
arrow.p$FG<- c("OM (%)","Total N", "K", "Mg", "Ca", "S","pH", "H(meq/100g)", "CEC(meq/100g)", "K+", 
                                        "Mg+2", "Ca+2", "H+", "Na+")

environmental_plot = ggplot(NMDS.plot, aes(x=NMDS1,y=NMDS2)) +
  geom_point(size=3, stroke=3, aes(color=HabitatType,shape=HabitatType),position=position_jitter(0.1)) +
  scale_color_manual(values=c("#88A550","#336B87")) +
  scale_shape_manual(values=c(21,23)) +
  stat_ellipse(aes(color=HabitatType,lty=HabitatType), type='t',size =2)+
  theme(axis.text.x=element_text(colour="black",size=17)) +
  theme(axis.text.y=element_text(colour="black",size=17)) +
  theme(axis.title=element_text(colour="black",size=20)) +
  theme(panel.border = element_rect(fill=NA),axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(strip.text.x=element_text(size=20,face="bold"),strip.background = element_rect(fill="white")) +
  theme(legend.text=element_text(size=20)) +
  theme(legend.title=element_text(size=20,face='bold')) +
  geom_segment(data=arrow.p,size=2,aes(x=0,y=0,xend=NMDS1/2,yend=NMDS2/2,shape=NULL),
               arrow=arrow(length=unit(3,"cm")*arrow.p$R),colour="black") +
  geom_text(mapping=label_map,data=arrow.p,show.legend=FALSE,size=10)

plot(environmental_plot)
ggsave("Fig_6_NMDS_environmental_plot.tiff",width= 20,height=15,limitsize=FALSE,plot=environmental_plot)



###Stats###
#Bray PERMANOVA 
haka.bc.adonis <- adonis(bc_dist~HabitatType*Host,data=root_metadata,permutations = 9999)
haka.bc.adonis

#Pairwise permanova
bc_dist <- distance(vt_rel_abund,method="bray")

haka.pair.perm <- pairwise.perm.manova(bc_dist,root_metadata$Hab.by.Host,
                                                    nperm=9999,p.method="fdr")
haka.pair.perm

##Bray Betadisper
#Total difference
haka.beta.disper <- betadisper(bc_dist,root_metadata$Hab.by.Host)
haka.beta.disper.results <- permutest(haka.beta.disper, pairwise = TRUE, iter=9999)
haka.beta.disper.results

####################################
#    Indicator Species Analysis    #
####################################
library(indicspecies)
host = as.vector(root_metadata$Host)
comm_mat <- as.data.frame(root_otu)
host_indval = multipatt(comm_mat, host, control=how(nperm=9999))
host_indval_sum <- summary(host_indval,indvalcomp=TRUE)
host_sign_assoc = signassoc(comm_mat, cluster=host, control=how(nperm=9999))
host_sign_assoc


hab = as.vector(root_metadata$HabitatType)
hab_indval = multipatt(comm_mat, hab, control=how(nperm=9999))
hab_indval_sum <- summary(hab_indval,indvalcomp=TRUE)
hab_sign_assoc = signassoc(comm_mat, cluster=hab, control=how(nperm=9999))

                    #########################################
                    #    Species abundance distributions    #
                    #########################################
# Generate and analyze species abundance curves of using sequence abundance on rarefied OTU table. Remember that each
# sample is rarefied to 
install.packages('sads')
library(sads);library(dplyr)
remnant_spec_abund <- subset_samples(vt_physeq, HabitatType=="Remnant Forest")
remnant_otu <- otu_table(remnant_spec_abund)
remnant_abund <- as.data.frame(rowSums(t(remnant_otu)))
remnant_abund <- remnant_abund[remnant_abund$`rowSums(t(remnant_otu))` !=0, ]
remnant_abund <- as.integer(remnant_abund)

restored_spec_abund <- subset_samples(vt_physeq, HabitatType=="Restored Forest")
restored_otu <- otu_table(restored_spec_abund)
restored_abund <- as.data.frame(rowSums(t(restored_otu)))
restored_abund <- restored_abund[restored_abund$`rowSums(t(restored_otu))` !=0, ]
restored_abund <- as.integer(restored_abund)

#Octav tabulates the number of speceis in classes of logarithm of abundances at base 2 
# (Preston's octaves) and returns a dataframe
remnant.oc <- octav(remnant_abund)
restored.oc <- octav(restored_abund)
# Plot sads
par(mfrow=c(2,1))
par(mar=c(5,5,2,2))
plot(remnant.oc,main="Remnant forest")
plot(restored.oc,main="Restored forest")

(remnant.lnorm<-fitsad(x=remnant_abund,sad="lnorm"))
(remnant.lnorm.oc <- octavpred(remnant.lnorm))
(remnant.ls<-fitsad(x=remnant_abund,sad="ls"))
(remnant.ls.oc <- octavpred(remnant.ls))
(remnant.poilog<-fitsad(x=remnant_abund,sad="poilog"))
(remnant.poilog.oc<-octavpred(remnant.poilog))
(remnant.power<-fitsad(x=remnant_abund,sad="power"))
(remnant.power.oc<-octavpred(remnant.power))
remnant.SAD.AIC <- AICtab(remnant.lnorm,remnant.ls,remnant.poilog,remnant.power,
                          base=TRUE,weights=TRUE)
remnant.SAD.AIC
summary(remnant.poilog)

(restored.lnorm<-fitsad(x=restored_abund,sad="lnorm"))
(restored.lnorm.oc <- octavpred(restored.lnorm))
(restored.ls<-fitsad(x=restored_abund,sad="ls"))
(restored.ls.oc <- octavpred(restored.ls))
(restored.poilog<-fitsad(x=restored_abund,sad="poilog"))
(restored.poilog.oc<-octavpred(restored.poilog))
(restored.power<-fitsad(x=restored_abund,sad="power"))
(restored.power.oc<-octavpred(restored.power))

restored.SAD.AIC <- AICtab(restored.lnorm,restored.ls,restored.poilog,restored.power,
                          base=TRUE,weights=TRUE)
restored.SAD.AIC
summary(restored.poilog)


par(mar=c(5,5,2,2))
plot(remnant.oc)
title("A) Remnant forest",adj=1)
lines(remnant.lnorm.oc,col="#336B87",lwd=2)
lines(remnant.ls.oc,col="#598234",lwd=2)
lines(remnant.poilog.oc,col="#A43820",lwd=2)
lines(remnant.power.oc,col="#C99E10",lwd=2)

plot(restored.oc)
title("B) Restored forest",adj=1)
lines(restored.lnorm.oc,col="#336B87",lwd=2)
lines(restored.ls.oc,col="#598234",lwd=2)
lines(restored.poilog.oc,col="#A43820",lwd=2)
lines(restored.power.oc,col="#C99E10",lwd=2)

legend('right',
       c("Lognormal","Fisher's log-series","Poisson-lognormal","Power-law"),
       col=c("#336B87","#598234","#A43820","#C99E10"),
      lty=1,cex=1,lwd=6)

                                                      ############################   
                                                      #      Soil Chemistry      #
                                                      ############################ 
#Soil chemistry differences among habitat types and hosts
str(root_metadata)
OM_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$OM,
                                 subset(root_metadata,HabitatType == "Restored Forest")$OM,
                                 paired=FALSE, var.equal=FALSE)
OM_test

N_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$N,
                  subset(root_metadata,HabitatType == "Restored Forest")$N,
                  paired=FALSE, var.equal=FALSE)
N_test

P_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$P,
                  subset(root_metadata,HabitatType == "Restored Forest")$P,
                  paired=FALSE, var.equal=FALSE)
P_test

K_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$K,
                 subset(root_metadata,HabitatType == "Restored Forest")$K,
                 paired=FALSE, var.equal=FALSE)
K_test

Mg_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$Mg,
                  subset(root_metadata,HabitatType == "Restored Forest")$Mg,
                  paired=FALSE, var.equal=FALSE)
Mg_test

Ca_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$Ca,
                  subset(root_metadata,HabitatType == "Restored Forest")$Ca,
                  paired=FALSE, var.equal=FALSE)
Ca_test

Na_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$Na,
                  subset(root_metadata,HabitatType == "Restored Forest")$Na,
                  paired=FALSE, var.equal=FALSE)
Na_test

S_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$S,
                  subset(root_metadata,HabitatType == "Restored Forest")$S,
                  paired=FALSE, var.equal=FALSE)
S_test

pH_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$pH,
                  subset(root_metadata,HabitatType == "Restored Forest")$pH,
                  paired=FALSE, var.equal=FALSE)
pH_test

H_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$H,
                  subset(root_metadata,HabitatType == "Restored Forest")$H,
                  paired=FALSE, var.equal=FALSE)
H_test

CEC_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$CEC,
                  subset(root_metadata,HabitatType == "Restored Forest")$CEC,
                  paired=FALSE, var.equal=FALSE)
CEC_test

K.cation_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$K.cation,
                  subset(root_metadata,HabitatType == "Restored Forest")$K.cation,
                  paired=FALSE, var.equal=FALSE)
K.cation_test

Mg.cation_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$Mg.cation,
                  subset(root_metadata,HabitatType == "Restored Forest")$Mg.cation,
                  paired=FALSE, var.equal=FALSE)
Mg.cation_test

Ca.cation_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$Ca.cation,
                  subset(root_metadata,HabitatType == "Restored Forest")$Ca.cation,
                  paired=FALSE, var.equal=FALSE)
Ca.cation_test

H.cation_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$H.cation,
                  subset(root_metadata,HabitatType == "Restored Forest")$H.cation,
                  paired=FALSE, var.equal=FALSE)
H.cation_test

Na.cation_test <- t.test(subset(root_metadata,HabitatType == "Remnant Forest")$Na.cation,
                  subset(root_metadata,HabitatType == "Restored Forest")$Na.cation,
                  paired=FALSE, var.equal=FALSE)
Na.cation_test

                                                  #######################
                                                  ###Bipartite Network###
                                                  #######################
#load packages
library(sna);library(bipartite);library(Hmisc);library(ggpubr); library(Rmisc); library(emmeans);library(multcompView); 
library(multcomp); require(gridExtra); library(tidyverse)

#Subset vt.melt to include only root samples for bipartite analysis
haka_root_df<-vt.melt[which(vt.melt$SampleType=='roots'),]

#Convert a datatable into a network matrix for further use in bipartite.
#Create output as "list" type to calculate individual metrics and generate individual graphs for each plot.
#Prefereabble for different webs do not have to include all species names, i.e. they can be of different dimensions 
#(ragged). As such they are better suited for webs with non-comparable species sets (differing community rich/comps).
#Generates a weighted matrix for each plot
haka_bipartite <- frame2webs(haka_root_df,varnames=c("Species","Host","Plot","Abundance"),type.out="list")

haka_composite_bipartite<-frame2webs(haka_root_df,varnames=c("Species","Host","HabitatType","Abundance"),
                                     type.out="list")

####################################################
##  igraph observed composite bipartite networks  ##
####################################################
#Afforested koa nodes sized by degree centrality:
ak_bi <- graph.incidence(haka_composite_bipartite$`Afforested koa`, weighted = T)
V(ak_bi)$label<-V(ak_bi)$name
V(ak_bi)$frame.color = "black"
deg=centr_degree(ak_bi,mode="all")
V(ak_bi)$size = sqrt(deg$res)
ak_cust_layout<-layout.reingold.tilford(ak_bi,circular=T)
#Colour by Host. Need to assign colours to match previously made graphs
V(ak_bi)$Host=as.character(sample$Host[match(V(ak_bi)$name,sample$Host)])
V(ak_bi)$Host
V(ak_bi)$color=V(ak_bi)$Host
V(ak_bi)$color=gsub("Metrosideros polymorpha","#D73027",V(ak_bi)$color)
V(ak_bi)$color=gsub("Acacia koa","#FC8D59",V(ak_bi)$color)
V(ak_bi)$color=gsub("Coprosma rhynchocarpa","#FEE090",V(ak_bi)$color)
V(ak_bi)$color=gsub("Myrsine lessertiana","#008000",V(ak_bi)$color)
V(ak_bi)$color=gsub("Cheirodendron trigynum","#E0F3F8",V(ak_bi)$color)
V(ak_bi)$color=gsub("Rubus hawaiiensis","#008B8B",V(ak_bi)$color)
V(ak_bi)$color=gsub("Grass","#4575B4",V(ak_bi)$color)
V(ak_bi)$color=gsub(NA,"white",V(ak_bi)$color)

par(mar=c(0,0,0,0))
#With labels
set.seed(606)
plot(ak_bi, vertex.shape="circle",vertex.color=V(ak_bi)$color,
                     edge.width=1, edge.color="grey",layout=ak_cust_layout,
                     vertex.label=ifelse(degree(ak_bi)>10,V(ak_bi)$label,NA),vertex.label.font=2,
                     vertex.label.color="black",vertex.label.cex=1.5,vertex.label.dist=0)
#Legend to be used with networks coloured by AMF family
legend(title="Host",'topleft',
       c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
                    "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),pch=21,
       pt.bg=c("#D73027","#FC8D59","#FEE090","#008000","#E0F3F8","#008B8B","#4575B4"),pt.cex=3,
       box.lty=0,cex=0.75,text.font=4)

#Without labels
set.seed(606)
plot(ak_bi, vertex.shape="circle",vertex.color=col[as.numeric(V(ak_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=ak_cust_layout,
     vertex.label=NA,
     vertex.label.color="black",vertex.label.cex=1.5,vertex.label.dist=0)
#Legend to be used with networks coloured by AMF family
legend(title="Host",'topleft',
       c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
         "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),pch=21,
       pt.bg=c("#D73027","#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB","#4575B4"),pt.cex=3,
       box.lty=0,cex=0.75,text.font=4)

#save as .tiff file width = 1000

#Remnant ohia
ro_bi <- graph.incidence(haka_composite_bipartite$`Remnant ohia`, weighted = T)
V(ro_bi)$label<-V(ro_bi)$name
V(ro_bi)$frame.color = "black"
deg=centr_degree(ro_bi,mode="all")
V(ro_bi)$size = sqrt(deg$res)
ro_cust_layout<-layout.reingold.tilford(ro_bi,circular=T)
#Colour by Host. Need to assign colours to match previously made graphs
V(ro_bi)$Host=as.character(sample$Host[match(V(ro_bi)$name,sample$Host)])
V(ro_bi)$Host <- as.character(V(ro_bi)$Host)
V(ro_bi)$Host <- ifelse(is.na(V(ro_bi)$Host),'AMF',V(ro_bi)$Host)
V(ro_bi)$Host
V(ro_bi)$color=V(ro_bi)$Host
V(ro_bi)$color=gsub("Metrosideros polymorpha","#D73027",V(ro_bi)$color)
V(ro_bi)$color=gsub("Acacia koa","#FC8D59",V(ro_bi)$color)
V(ro_bi)$color=gsub("Coprosma rhynchocarpa","#FEE090",V(ro_bi)$color)
V(ro_bi)$color=gsub("Myrsine lessertiana","#008000",V(ro_bi)$color)
V(ro_bi)$color=gsub("Cheirodendron trigynum","#E0F3F8",V(ro_bi)$color)
V(ro_bi)$color=gsub("Rubus hawaiiensis","#008B8B",V(ro_bi)$color)
V(ro_bi)$color=gsub("Grass","#4575B4",V(ro_bi)$color)
V(ro_bi)$color=gsub("AMF","white",V(ro_bi)$color)

#With labels
set.seed(605)
RO_composite_bipartite <- plot(ro_bi, vertex.color=V(ro_bi)$color,
     edge.width=1, edge.color="grey",layout=ro_cust_layout,
     vertex.label=ifelse(degree(ro_bi)>10,V(ro_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#Legend to be used with networks coloured by AMF family
legend(title="Host",'topleft',
       c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
         "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),pch=21,
       pt.bg=c("#D73027","#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB","#4575B4"),pt.cex=3,
       box.lty=0,cex=0.75,text.font=4)

#Without labels
set.seed(605)
RO_composite_bipartite <- plot(ro_bi,vertex.color=col[as.numeric(V(ro_bi)$type)+1],
                               edge.width=1, edge.color="grey",layout=ro_cust_layout,
                               vertex.label=NA)
#Legend to be used with networks coloured by AMF family
legend(title="Host",'topleft',
       c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
         "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),pch=21,
       pt.bg=c("#D73027","#FC8D59","#FEE090","#FFFFBF","#E0F3F8","#91BFDB","#4575B4"),pt.cex=3,
       box.lty=0,cex=0.75,text.font=4)
#save as .tiff file width = 1000

###########################
##  Bipartite analyses   ##
###########################
#Generate null networks for each plot to be used for null model analyses. swap.web function used by Sepp et al., 2018
#Sepp et al. 2018 used this method in their paper examining plant AMF networks
#The "quasiswapcount" algorithm was chosen due to its property of maintaining marginal totals and network connectance. 
#The randomized matrices should not, therefore, exceed the realistic symbiont holding capacities of plant species.
# It is equivalent to the using the quasiswapcount algorithm and permatswap function in vegan
null_ak1 <- swap.web(haka_bipartite$AK1,N=1000)
null_ak2 <- swap.web(haka_bipartite$AK2,N=1000)
null_ak3 <- swap.web(haka_bipartite$AK3,N=1000)
null_ak4 <- swap.web(haka_bipartite$AK4,N=1000)
null_ak5 <- swap.web(haka_bipartite$AK5,N=1000)
null_ak6 <- swap.web(haka_bipartite$AK6,N=1000)
null_ro1 <- swap.web(haka_bipartite$RO1,N=1000)
null_ro2 <- swap.web(haka_bipartite$RO2,N=1000)
null_ro3 <- swap.web(haka_bipartite$RO3,N=1000)
null_ro4 <- swap.web(haka_bipartite$RO4,N=1000)
null_ro5 <- swap.web(haka_bipartite$RO5,N=1000)
null_ro6 <- swap.web(haka_bipartite$RO6,N=1000)
  
#Create labels for null  values (to be used for null dataframes created later)
null_hab_type_labs<-as.data.frame(rep(c("Afforested koa","Remnant ohia"),each=6000))
null_plot_labs<-as.data.frame(rep(c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6"),each=1000))


#######################
###   Connectance   ###
#######################
connectance <- lapply(haka_bipartite,networklevel, index='connectance')

con_dat <- data.frame(unlist(connectance))
rownames(con_dat) <- c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6")
con_dat <-cbind(con_dat,c("Afforested koa","Afforested koa","Afforested koa","Afforested koa",
                          "Afforested koa","Afforested koa","Remnant ohia","Remnant ohia","Remnant ohia",
                          "Remnant ohia","Remnant ohia","Remnant ohia"))
con_dat <-cbind(con_dat,c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6"))
colnames(con_dat) <- c("connectance","HabitatType","plot")

conn_p <- ggplot(con_dat,aes(x=HabitatType,y=connectance,fill=HabitatType)) +
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",width=0.1,size=2) +
  stat_summary(fun.y=mean,geom="point",size=6,shape=21,stroke=2) +
  ylab("Network Connectance") + xlab("Habitat Type") +
  ylim(0.2,0.6) +
  xlab("") +
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("A)") 

plot(conn_p)
ggsave("haka_connectance.tiff", 
       plot = conn_p)

#Observed welch t-test
connectance_welch <- t.test(subset(con_dat,HabitatType == "Afforested koa")$connectance,
                           subset(con_dat,HabitatType == "Remnant ohia")$connectance,
                           paired=FALSE, var.equal=FALSE)
connectance_welch

#Observed vs null network analysis - need to use a different null model than swap.web because swap.web constrains
#network connectance to be the same as the observed networks. Chose r2dtable as it is complete randomization.
#r2dtable maintains margnial sums but randomizes everything else (including connectacne)
null_r2d_ak1 <- nullmodel(haka_bipartite$AK1,method=1,N=1000)
null_r2d_ak2 <- nullmodel(haka_bipartite$AK2,method=1, N=1000)
null_r2d_ak3 <- nullmodel(haka_bipartite$AK3,method=1, N=1000)
null_r2d_ak4 <- nullmodel(haka_bipartite$AK4,method=1, N=1000)
null_r2d_ak5 <- nullmodel(haka_bipartite$AK5,method=1, N=1000)
null_r2d_ak6 <- nullmodel(haka_bipartite$AK6,method=1, N=1000)
null_r2d_ro1 <- nullmodel(haka_bipartite$RO1,method=1, N=1000)
null_r2d_ro2 <- nullmodel(haka_bipartite$RO2,method=1, N=1000)
null_r2d_ro3 <- nullmodel(haka_bipartite$RO3,method=1, N=1000)
null_r2d_ro4 <- nullmodel(haka_bipartite$RO4,method=1, N=1000)
null_r2d_ro5 <- nullmodel(haka_bipartite$RO5,method=1, N=1000)
null_r2d_ro6 <- nullmodel(haka_bipartite$RO6,method=1, N=1000)

connect_null_ak1<-unlist(sapply(null_r2d_ak1,networklevel,index='connectance'))
connect_null_ak2<-unlist(sapply(null_r2d_ak2,networklevel,index='connectance'))
connect_null_ak3<-unlist(sapply(null_r2d_ak3,networklevel,index='connectance'))
connect_null_ak4<-unlist(sapply(null_r2d_ak4,networklevel,index='connectance'))
connect_null_ak5<-unlist(sapply(null_r2d_ak5,networklevel,index='connectance'))
connect_null_ak6<-unlist(sapply(null_r2d_ak6,networklevel,index='connectance'))
connect_null_ro1<-unlist(sapply(null_r2d_ro1,networklevel,index='connectance'))
connect_null_ro2<-unlist(sapply(null_r2d_ro2,networklevel,index='connectance'))
connect_null_ro3<-unlist(sapply(null_r2d_ro3,networklevel,index='connectance'))
connect_null_ro4<-unlist(sapply(null_r2d_ro4,networklevel,index='connectance'))
connect_null_ro5<-unlist(sapply(null_r2d_ro5,networklevel,index='connectance'))
connect_null_ro6<-unlist(sapply(null_r2d_ro6,networklevel,index='connectance'))

connect_null<-list(c(AK1=connect_null_ak1,AK2=connect_null_ak2,AK3=connect_null_ak3,AK4=connect_null_ak4,
                     AK5=connect_null_ak5,AK6=connect_null_ak6,RO1=connect_null_ro1,RO2=connect_null_ro2,
                     RO3=connect_null_ro3,RO4=connect_null_ro4,RO5=connect_null_ro5,RO6=connect_null_ro6))

con_null_df <- as.data.frame(connect_null)
write.csv(as_tibble(con_null_df ),file="conn_null_df.csv")

con_null_df<-cbind(con_null_df,null_hab_type_labs,null_plot_labs)
colnames(con_null_df) <- c("connectance","HabitatType","plot")
connect_null_df<- summarySE(con_null_df, measurevar="connectance", 
                            groupvars=c("HabitatType"))

connect_null_plot <- ggplot(data=connect_null_df,aes(x=HabitatType,y=connectance)) +
  geom_errorbar(aes(ymin=connectance-sd,ymax=connectance+sd),width=0.1,colour="black",size=2) +
  geom_point(colour="black",pch=21,fill="white",size=5,stroke=2) +
  stat_summary(data=con_dat,aes(x=HabitatType,y=connectance,stat="identity"),
               fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(data=con_dat,aes(x=HabitatType,y=connectance,stat="identity"),
               fun.y=mean,geom="point",colour="black",size=5) +
    ylab("Network Connectance") + xlab("HabitatType") +
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,1)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) 

plot(connect_null_plot)

#t-tests comparing observed to null connectance
ak.conn.null.t <- t.test(subset(con_dat,HabitatType == "Afforested koa")$connectance,
                                subset(con_null_df,HabitatType == "Afforested koa")$connectance,
                                paired=FALSE, var.equal=FALSE)
ak.conn.null.t

ro.conn.null.t <- t.test(subset(con_dat,HabitatType == "Remnant ohia")$connectance,
                         subset(con_null_df,HabitatType == "Remnant ohia")$connectance,
                         paired=FALSE, var.equal=FALSE)
ro.conn.null.t


##############################
###   Network nestedness   ###
##############################
nestedness <-lapply(haka_bipartite,networklevel,index='NODF')

nest_dat <- data.frame(unlist(nestedness))
rownames(nest_dat) <- c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6")
nest_dat <-cbind(nest_dat,c("Afforested koa","Afforested koa","Afforested koa","Afforested koa",
                            "Afforested koa","Afforested koa","Remnant ohia","Remnant ohia","Remnant ohia",
                            "Remnant ohia","Remnant ohia","Remnant ohia"))
nest_dat <-cbind(nest_dat,c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6"))
colnames(nest_dat) <- c("nestedness","HabitatType","plot")

nestp <- ggplot(nest_dat,aes(x=HabitatType,y=nestedness)) +   
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(fun.y=mean,geom="point",color="black",size=6,shape=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Network Nestedness") + xlab("Habitat Type") +
  ylim(25,75) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("B)") 

plot(nestp)
ggsave("haka_nestedness.tiff", 
       plot = nestp)

#Observed welch t-test
nestedness_welch <- t.test(subset(nest_dat,HabitatType == "Afforested koa")$nestedness,
                            subset(nest_dat,HabitatType == "Remnant ohia")$nestedness,
                            paired=FALSE, var.equal=FALSE)
nestedness_welch

#Observed vs null network analysis
nest_null_ak1<-unlist(sapply(null_ak1,networklevel,index='NODF'))
nest_null_ak2<-unlist(sapply(null_ak2,networklevel,index='NODF'))
nest_null_ak3<-unlist(sapply(null_ak3,networklevel,index='NODF'))
nest_null_ak4<-unlist(sapply(null_ak4,networklevel,index='NODF'))
nest_null_ak5<-unlist(sapply(null_ak5,networklevel,index='NODF'))
nest_null_ak6<-unlist(sapply(null_ak6,networklevel,index='NODF'))
nest_null_ro1<-unlist(sapply(null_ro1,networklevel,index='NODF'))
nest_null_ro2<-unlist(sapply(null_ro2,networklevel,index='NODF'))
nest_null_ro3<-unlist(sapply(null_ro3,networklevel,index='NODF'))
nest_null_ro4<-unlist(sapply(null_ro4,networklevel,index='NODF'))
nest_null_ro5<-unlist(sapply(null_ro5,networklevel,index='NODF'))
nest_null_ro6<-unlist(sapply(null_ro6,networklevel,index='NODF'))

nest_null<-list(c(AK1=nest_null_ak1,AK2=nest_null_ak2,AK3=nest_null_ak3,AK4=nest_null_ak4,AK5=nest_null_ak5,
                  AK6=nest_null_ak6,RO1=nest_null_ro1,RO2=nest_null_ro2,RO3=nest_null_ro3,RO4=nest_null_ro4,
                  RO5=nest_null_ro5,RO6=nest_null_ro6))

#create a null datafrom from the list of null nested values above
nest_null_df <- as.data.frame(nest_null)
write.csv(as_tibble(nest_null_df ),file="nest_null_df.csv")

#Add metadata label columns to null nestedness dataframe
nest_null_df<-cbind(nest_null_df,null_hab_type_labs,null_plot_labs)
colnames(nest_null_df) <- c("nestedness","HabitatType","plot")
nest_null_sum<- summarySE(nest_null_df, measurevar="nestedness", 
                         groupvars=c("HabitatType"))

nest_null_plot <- ggplot(data=nest_null_sum,aes(x=HabitatType,y=nestedness)) +
  geom_errorbar(aes(ymin=nestedness-sd,ymax=nestedness+sd),width=0.1,colour="black",size=2) +
  geom_point(colour="black",pch=21,fill="white",size=5,stroke=2) +
  stat_summary(data=nest_dat,aes(x=HabitatType,y=nestedness,stat="identity"),
               fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(data=nest_dat,aes(x=HabitatType,y=nestedness,stat="identity"),
               fun.y=mean,geom="point",colour="black",size=5) +
  ylab("Network Nestedness") + xlab("HabitatType") +
  scale_y_continuous(breaks=seq(0,100,by=10),limits=c(40,70)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) 

plot(nest_null_plot)
ggsave("null_nestedness_97_otu_table_Jan7.tiff", 
       plot = nest_null_plot)

#t-tests comparing observed to null nestedness
ak.nest.null.t <- t.test(subset(nest_dat,HabitatType == "Afforested koa")$nestedness,
                         subset(nest_null_df,HabitatType == "Afforested koa")$nestedness,
                         paired=FALSE, var.equal=FALSE)
ak.nest.null.t

ro.nest.null.t <- t.test(subset(nest_dat,HabitatType == "Remnant ohia")$nestedness,
                         subset(nest_null_df,HabitatType == "Remnant ohia")$nestedness,
                         paired=FALSE, var.equal=FALSE)
ro.nest.null.t


######################
###   Modularity   ###
######################
mod_ak1 <- NOS(haka_bipartite$AK1) 
mod_ak2 <- NOS(haka_bipartite$AK2)
mod_ak3 <- NOS(haka_bipartite$AK3)
mod_ak4 <- NOS(haka_bipartite$AK4)
mod_ak5 <- NOS(haka_bipartite$AK5)
mod_ak6 <- NOS(haka_bipartite$AK6)
mod_ro1 <- NOS(haka_bipartite$RO1)
mod_ro2 <- NOS(haka_bipartite$RO2)
mod_ro3 <- NOS(haka_bipartite$RO3)
mod_ro4 <- NOS(haka_bipartite$RO4)
mod_ro5 <- NOS(haka_bipartite$RO5)
mod_ro6 <- NOS(haka_bipartite$RO6)

mod_list <-list(c(AK1=mod_ak1$mod,AK2= mod_ak2$mod,AK3=mod_ak3$mod,
                  AK4=mod_ak4$mod,AK5=mod_ak5$mod,AK6=mod_ak6$mod,
                  RO1=mod_ro1$mod,RO2=mod_ro2$mod,RO3=mod_ro3$mod,
                  RO4=mod_ro4$mod,RO5=mod_ro5$mod,RO6=mod_ro6$mod))
mod_dat <- as.data.frame(mod_list)

mod_dat <-cbind(mod_dat,c("Afforested koa","Afforested koa","Afforested koa","Afforested koa",
                        "Afforested koa","Afforested koa",
                        "Remnant ohia","Remnant ohia","Remnant ohia",
                        "Remnant ohia","Remnant ohia","Remnant ohia"))
mod_dat <-cbind(mod_dat,c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6"))
colnames(mod_dat) <- c("modularity","HabitatType","plot")

mod_p <- ggplot(mod_dat,aes(x=HabitatType,y=modularity)) +
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(fun.y=mean,geom="point",color="black",size=6,pch=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Network Modularity (Q)") + xlab("Habitat Type") +
  ylim(0.25,0.75) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("C)") 

plot(mod_p)
ggsave("haka_modularity.tiff", 
       plot = mod_p)

#Observed welch t-test
modularity_welch <- t.test(subset(mod_dat,HabitatType == "Afforested koa")$modularity,
                           subset(mod_dat,HabitatType == "Remnant ohia")$modularity,
                           paired=FALSE, var.equal=FALSE)
modularity_welch


#Observed vs null network analysis (takes a long time!!)
mod_null_ak1<-unlist(sapply(null_ak1,NOS))
mod_null_ak2<-unlist(sapply(null_ak2,NOS))
mod_null_ak3<-unlist(sapply(null_ak3,NOS))
mod_null_ak4<-unlist(sapply(null_ak4,NOS))
mod_null_ak5<-unlist(sapply(null_ak5,NOS))
mod_null_ak6<-unlist(sapply(null_ak6,NOS))
mod_null_ro1<-unlist(sapply(null_ro1,NOS))
mod_null_ro2<-unlist(sapply(null_ro2,NOS))
mod_null_ro3<-unlist(sapply(null_ro3,NOS))
mod_null_ro4<-unlist(sapply(null_ro4,NOS))
mod_null_ro5<-unlist(sapply(null_ro5,NOS))
mod_null_ro6<-unlist(sapply(null_ro6,NOS))

mod_null<-list(c(AK1=mod_null_ak1,AK2=mod_null_ak2,AK3=mod_null_ak3,AK4=mod_null_ak4,AK5=mod_null_ak5,
                 AK6=mod_null_ak6,RO1=mod_null_ro1,RO2=mod_null_ro2,RO3=mod_null_ro3,RO4=mod_null_ro4,
                 RO5=mod_null_ro5,RO6=mod_null_ro6))
mod_null_df <- as.data.frame(mod_null)
write.csv(as_tibble(mod_null_df ),file="mod_null_df.csv")

mod_null_df<-cbind(mod_null_df,null_hab_type_labs,null_plot_labs)
colnames(mod_null_df) <- c("modularity","HabitatType","plot")
mod_null_sum<- summarySE(mod_null_df, measurevar="modularity", 
                        groupvars=c("HabitatType"))

mod_null_plot <- ggplot(data=mod_null_sum,aes(x=HabitatType,y=modularity)) +
   geom_errorbar(aes(ymin=modularity-ci,ymax=modularity+ci),width=0.1,colour="black",size=2) +
  geom_point(colour="black",pch=21,fill="white",size=5,stroke=2) +
  stat_summary(data=mod_dat,aes(x=HabitatType,y=modularity,stat="identity"),
               fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(data=mod_dat,aes(x=HabitatType,y=modularity,stat="identity"),
               fun.y=mean,geom="point",color="black",size=5) +
  ylab("Network Modularity") + xlab("HabitatType") +
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0.1,0.6)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) 

plot(mod_null_plot)
ggsave("haka_null_modularity.tiff", 
       plot = mod_null_plot)

#t-tests comparing observed to null modularity
ak.mod.null.t <- t.test(subset(mod_dat,HabitatType == "Afforested koa")$modularity,
                         subset(mod_null_df,HabitatType == "Afforested koa")$modularity,
                         paired=FALSE, var.equal=FALSE)
ak.mod.null.t

ro.mod.null.t <- t.test(subset(mod_dat,HabitatType == "Remnant ohia")$modularity,
                         subset(mod_null_df,HabitatType == "Remnant ohia")$modularity,
                         paired=FALSE, var.equal=FALSE)
ro.mod.null.t



###########################
###   Linkage Density   ###
###########################
link_dens <-lapply(haka_bipartite, networklevel, index='linkage density')
# Mean number of interactions per species
link_dens_dat <- data.frame(unlist(link_dens))

link_dens_dat <- cbind(link_dens_dat,rep(c("Restored","Remnant"),each=6))
link_dens_dat <-cbind(link_dens_dat,c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6"))

colnames(link_dens_dat) <- c("LinkageDensity","HabitatType","plot")
link_dens_dat$HabitatType <- factor(link_dens_dat$HabitatType, levels=c("Restored","Remnant"))

link_dens_plot <- ggplot(link_dens_dat,aes(x=HabitatType,y=LinkageDensity)) +   
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(fun.y=mean,geom="point",color="black",size=6,pch=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Linkage Density") + xlab("Habitat Type") +
  ylim(6,10) +
  theme_classic() +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=20)) +
  theme(axis.text.y=element_text(colour="black",size=20)) +
  theme(axis.title.y=element_text(colour="black",size=22)) +
  theme(axis.title.x=element_text(colour="white")) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("B)") 

plot(link_dens_plot)
ggsave("Fig_1B_linkage_density.tiff", 
       plot = link_dens_plot)

#Observed welch t-test
link_dens_welch <- t.test(subset(link_dens_dat,HabitatType == "Restored")$LinkageDensity,
                           subset(link_dens_dat,HabitatType == "Remnant")$LinkageDensity,
                           paired=FALSE, var.equal=FALSE)
link_dens_welch

#Observed vs null analysis
link_dens_null_ak1<-unlist(sapply(null_ak1,networklevel,index='linkage density'))
link_dens_null_ak2<-unlist(sapply(null_ak2,networklevel,index='linkage density'))
link_dens_null_ak3<-unlist(sapply(null_ak3,networklevel,index='linkage density'))
link_dens_null_ak4<-unlist(sapply(null_ak4,networklevel,index='linkage density'))
link_dens_null_ak5<-unlist(sapply(null_ak5,networklevel,index='linkage density'))
link_dens_null_ak6<-unlist(sapply(null_ak6,networklevel,index='linkage density'))
link_dens_null_ro1<-unlist(sapply(null_ro1,networklevel,index='linkage density'))
link_dens_null_ro2<-unlist(sapply(null_ro2,networklevel,index='linkage density'))
link_dens_null_ro3<-unlist(sapply(null_ro3,networklevel,index='linkage density'))
link_dens_null_ro4<-unlist(sapply(null_ro4,networklevel,index='linkage density'))
link_dens_null_ro5<-unlist(sapply(null_ro5,networklevel,index='linkage density'))
link_dens_null_ro6<-unlist(sapply(null_ro6,networklevel,index='linkage density'))

link_dens_null<-list(c(AK1=link_dens_null_ak1,AK2=link_dens_null_ak2,AK3=link_dens_null_ak3,AK4=link_dens_null_ak4,
                       AK5=link_dens_null_ak5,AK6=link_dens_null_ak6,RO1=link_dens_null_ro1,RO2=link_dens_null_ro2,
                       RO3=link_dens_null_ro3,RO4=link_dens_null_ro4,RO5=link_dens_null_ro5,RO6=link_dens_null_ro6))

#create a null datafrom from the list of null nested values above
link_dens_null_df <- as.data.frame(link_dens_null)
write.csv(as_tibble(link_dens_null_df),file="link_dens_null_df.csv")

link_dens_null_df<-cbind(link_dens_null_df,null_hab_type_labs,null_plot_labs)
colnames(link_dens_null_df) <- c("LinkageDensity","HabitatType","plot")
link_dens_null_sum<- summarySE(link_dens_null_df, measurevar="LinkageDensity", 
                          groupvars=c("HabitatType"))

link_dens_null_plot <- ggplot(data=link_dens_null_sum,aes(x=HabitatType,y=LinkageDensity)) +
  geom_point(colour="black",pch=21,fill="white",size=5,stroke=2) +
  geom_errorbar(aes(ymin=LinkageDensity-sd,ymax=LinkageDensity+sd),width=0.1,colour="black",size=2) +
  stat_summary(data=link_dens_dat,aes(x=HabitatType,y=LinkageDensity,stat="identity"),
               fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(data=link_dens_dat,aes(x=HabitatType,y=LinkageDensity,stat="identity"),
               fun.y=mean,geom="point",colour="black",size=5) +
  ylab("Linkage Density") + xlab("Habitat Type") +
  scale_y_continuous(breaks=seq(0,20,by=2),limits=c(4,12)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) 

plot(link_dens_null_plot)
ggsave("null_nestedness_97_otu_table_Jan7.tiff", 
       plot = nest_null_plot)

#t-tests comparing observed to null Linkage Density
ak.link.dens.null.t <- t.test(subset(link_dens_dat,HabitatType == "Afforested koa")$LinkageDensity,
                         subset(link_dens_null_df,HabitatType == "Afforested koa")$LinkageDensity,
                         paired=FALSE, var.equal=FALSE)
ak.link.dens.null.t

ro.link.dens.null.t <- t.test(subset(link_dens_dat,HabitatType == "Remnant ohia")$LinkageDensity,
                         subset(link_dens_null_df,HabitatType == "Remnant ohia")$LinkageDensity,
                         paired=FALSE, var.equal=FALSE)
ro.link.dens.null.t

###################################
###   Links per Plant Species   ###
###################################
link <-lapply(haka_bipartite, grouplevel, level="higher", index='mean number of links')

link_dat <- data.frame(unlist(link))

link_dat <- cbind(link_dat,rep(c("Restored","Remnant"),each=6))
link_dat <-cbind(link_dat,rep(c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6")))
colnames(link_dat) <- c("Mean.Num.Links","HabitatType","plot")
link_dat$HabitatType <- factor(link_dat$HabitatType, levels=c("Restored","Remnant"))

plant_links <- ggplot(link_dat,aes(x=HabitatType,y=Mean.Num.Links)) +   
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(fun.y=mean,geom="point",color="black",size=6,pch=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Mean Number of Links") + xlab("Habitat Type") +
  scale_y_continuous(breaks=seq(40,54,by=2),limits=c(43,55)) +
  theme_classic() +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=20)) +
  theme(axis.text.y=element_text(colour="black",size=20)) +
  theme(axis.title.y=element_text(colour="black",size=22)) +
  theme(axis.title.x=element_text(colour="white")) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("C)") 

plot(plant_links)
ggsave("Fig_1C_haka_plant_links.tiff", 
       plot = plant_links)

#Observed welch t-test
link_welch <- t.test(subset(link_dat,HabitatType == "Restored")$Mean.Num.Links,
                          subset(link_dat,HabitatType == "Remnant")$Mean.Num.Links,
                          paired=FALSE, var.equal=FALSE)
link_welch

#Observed vs null analysis
links_null_ak1<-unlist(sapply(null_ak1,grouplevel, level="higher", index='mean number of links'))
links_null_ak2<-unlist(sapply(null_ak2,grouplevel, level="higher", index='mean number of links'))
links_null_ak3<-unlist(sapply(null_ak3,grouplevel, level="higher", index='mean number of links'))
links_null_ak4<-unlist(sapply(null_ak4,grouplevel, level="higher", index='mean number of links'))
links_null_ak5<-unlist(sapply(null_ak5,grouplevel, level="higher", index='mean number of links'))
links_null_ak6<-unlist(sapply(null_ak6,grouplevel, level="higher", index='mean number of links'))
links_null_ro1<-unlist(sapply(null_ro1,grouplevel, level="higher", index='mean number of links'))
links_null_ro2<-unlist(sapply(null_ro2,grouplevel, level="higher", index='mean number of links'))
links_null_ro3<-unlist(sapply(null_ro3,grouplevel, level="higher", index='mean number of links'))
links_null_ro4<-unlist(sapply(null_ro4,grouplevel, level="higher", index='mean number of links'))
links_null_ro5<-unlist(sapply(null_ro5,grouplevel, level="higher", index='mean number of links'))
links_null_ro6<-unlist(sapply(null_ro6,grouplevel, level="higher", index='mean number of links'))

links_null<-list(c(AK1=links_null_ak1,AK2=links_null_ak2,AK3=links_null_ak3,AK4=links_null_ak4,
                       AK5=links_null_ak5,AK6=links_null_ak6,RO1=links_null_ro1,RO2=links_null_ro2,
                       RO3=links_null_ro3,RO4=links_null_ro4,RO5=links_null_ro5,RO6=links_null_ro6))

#create a null datafrom from the list of null values above
links_null_df <- as.data.frame(links_null)
write.csv(as_tibble(links_null_df),file="links_null_df.csv")

links_null_df<-cbind(links_null_df,null_hab_type_labs,null_plot_labs)
colnames(links_null_df) <- c("Mean.Num.Links","HabitatType","plot")
links_null_sum<- summarySE(links_null_df, measurevar="Mean.Num.Links", 
                               groupvars=c("HabitatType"))

links_null_plot <- ggplot(data=links_null_sum,aes(x=HabitatType,y=Mean.Num.Links)) +
  geom_point(colour="black",pch=21,fill="white",size=5,stroke=2) +
  geom_errorbar(aes(ymin=Mean.Num.Links-sd,ymax=Mean.Num.Links+sd),width=0.1,colour="black",size=2) +
  stat_summary(data=link_dat,aes(x=HabitatType,y=Mean.Num.Links,stat="identity"),
               fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(data=link_dat,aes(x=HabitatType,y=Mean.Num.Links,stat="identity"),
               fun.y=mean,geom="point",colour="black",size=5) +
  ylab("Linkage Density") + xlab("Habitat Type") +
  scale_y_continuous(breaks=seq(62,100,by=4),limits=c(40,60)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) 

plot(links_null_plot)
ggsave("haka_plant_link_null.tiff", 
       plot = links_null_plot)

#t-tests comparing observed to null mean plant links
ak.link.null.t <- t.test(subset(link_dat,HabitatType == "Afforested koa")$Mean.Num.Links,
                              subset(links_null_df,HabitatType == "Afforested koa")$Mean.Num.Links,
                              paired=FALSE, var.equal=FALSE)
ak.link.null.t

ro.link.null.t <- t.test(subset(link_dat,HabitatType == "Remnant ohia")$Mean.Num.Links,
                              subset(links_null_df,HabitatType == "Remnant ohia")$Mean.Num.Links,
                              paired=FALSE, var.equal=FALSE)
ro.link.null.t

#############################################
###   Host robustness to AMF extinction   ###
#############################################
host_robustness <- lapply(haka_bipartite,grouplevel,index='robustness', level= "higher",nrep=1000)

robust_dat<-data.frame(unlist(host_robustness))
rownames(robust_dat) <- c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6")
robust_dat <-cbind(robust_dat,c("Afforested koa","Afforested koa","Afforested koa","Afforested koa",
                                          "Afforested koa","Afforested koa","Remnant ohia","Remnant ohia","Remnant ohia",
                                          "Remnant ohia","Remnant ohia","Remnant ohia"))
colnames(robust_dat) <- c("robustness","HabitatType")
robust_p <- ggplot(robust_dat,aes(x=HabitatType,y=robustness,fill=HabitatType)) +
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(fun.y=mean,geom="point",color="black",size=6,pch=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Host Robustness") + xlab("Habitat Type") +
  ylim(0.6,0.8) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("F)") 

plot(robust_p)
ggsave("haka_host_robustness.tiff", 
       plot = robust_p)

#Observed welch t-test
robust_welch <- t.test(subset(robust_dat,HabitatType == "Afforested koa")$robustness,
                     subset(robust_dat,HabitatType == "Remnant ohia")$robustness,
                     paired=FALSE, var.equal=FALSE)
robust_welch

#Observed vs null analysis
robust_null_ak1<-unlist(sapply(null_ak1,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ak2<-unlist(sapply(null_ak2,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ak3<-unlist(sapply(null_ak3,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ak4<-unlist(sapply(null_ak4,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ak5<-unlist(sapply(null_ak5,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ak6<-unlist(sapply(null_ak6,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ro1<-unlist(sapply(null_ro1,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ro2<-unlist(sapply(null_ro2,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ro3<-unlist(sapply(null_ro3,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ro4<-unlist(sapply(null_ro4,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ro5<-unlist(sapply(null_ro5,grouplevel, level="higher", index='robustness', nrep=1000))
robust_null_ro6<-unlist(sapply(null_ro6,grouplevel, level="higher", index='robustness', nrep=1000))

robust_null<-list(c(AK1=robust_null_ak1,AK2=robust_null_ak2,AK3=robust_null_ak3,AK4=robust_null_ak4,
                   AK5=robust_null_ak5,AK6=robust_null_ak6,RO1=robust_null_ro1,RO2=robust_null_ro2,
                   RO3=robust_null_ro3,RO4=robust_null_ro4,RO5=robust_null_ro5,RO6=robust_null_ro6))

#create a null dataframe from the list of null values above
robust_null_df <- as.data.frame(robust_null)
write.csv(as_tibble(robust_null_df),file="robust_null_df.csv")

robust_null_df<-cbind(robust_null_df,null_hab_type_labs,null_plot_labs)
colnames(robust_null_df) <- c("robustness","HabitatType","plot")
robust_null_sum<- summarySE(robust_null_df, measurevar="robustness", 
                           groupvars=c("HabitatType"))

robust_null_plot <- ggplot(data=robust_null_sum,aes(x=HabitatType,y=robustness)) +
  geom_point(colour="black",pch=21,fill="white",size=5,stroke=2) +
  geom_errorbar(aes(ymin=robustness-sd,ymax=robustness+sd),width=0.1,colour="black",size=2) +
  stat_summary(data=robust_dat,aes(x=HabitatType,y=robustness,stat="identity"),
               fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(data=robust_dat,aes(x=HabitatType,y=robustness,stat="identity"),
               fun.y=mean,geom="point",colour="black",size=5) +
  ylab("Host Robustness") + xlab("Habitat Type") +
  scale_y_continuous(breaks=seq(0.6,1,by=0.1),limits=c(0.6,1)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) 

plot(robust_null_plot)
ggsave("haka_host_robust_null.tiff", 
       plot = robust_null_plot)

#t-tests comparing observed to null mean plant links
ak.robust.null.t <- t.test(subset(robust_dat,HabitatType == "Afforested koa")$robustness,
                         subset(robust_null_df,HabitatType == "Afforested koa")$robustness,
                         paired=FALSE, var.equal=FALSE)
ak.robust.null.t

ro.robust.null.t <- t.test(subset(robust_dat,HabitatType == "Remnant ohia")$robustness,
                         subset(robust_null_df,HabitatType == "Remnant ohia")$robustness,
                         paired=FALSE, var.equal=FALSE)
ro.robust.null.t


##################################
###   Network specialization   ###
##################################
network_specialization <- lapply(haka_bipartite, networklevel, index= 'H2')

H2_dat <- data.frame(unlist(network_specialization))
rownames(H2_dat) <- c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6")
H2_dat <-cbind(H2_dat,c("Restored","Restored","Restored","Restored",
                                    "Restored","Restored","Remnant","Remnant","Remnant",
                                    "Remnant","Remnant","Remnant"))
H2_dat <-cbind(H2_dat,c("AK1","AK2","AK3","AK4","AK5","AK6","RO1","RO2","RO3","RO4","RO5","RO6"))
H2_dat$HabitatType <- factor(H2_dat$HabitatType, levels=c("Restored","Remnant"))
colnames(H2_dat) <- c("networkspecialization","HabitatType","plot")

H2_p <- ggplot(H2_dat,aes(x=HabitatType,y=networkspecialization,fill=HabitatType)) +
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(fun.y=mean,geom="point",color="black",size=6,pch=21,stroke=2,aes(fill=HabitatType)) +
  ylab("Network Specialization (H2')") + xlab("Habitat Type") +
  ylim(0.1,0.4) +
  theme_classic() +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=20)) +
  theme(axis.text.y=element_text(colour="black",size=20)) +
  theme(axis.title.y=element_text(colour="black",size=22)) +
  theme(axis.title.x=element_text(colour="white")) +
  theme(legend.position="none") +
  scale_fill_manual(values=c("#336B87","#88A550")) +
  stat_compare_means(size=10,label="p.signif",label.x=1.5) +
  ggtitle("A)") 

plot(H2_p)
ggsave("Fig_1A_haka_H2.tiff", 
       plot = H2_p)

#Observed welch t-test
H2_welch <- t.test(subset(H2_dat,HabitatType == "Restored")$networkspecialization,
                       subset(H2_dat,HabitatType == "Remnant")$networkspecialization,
                       paired=FALSE, var.equal=FALSE)
H2_welch

#Observed vs null network analysis
net_spec_null_ak1<-unlist(sapply(null_ak1,networklevel,index='H2'))
net_spec_null_ak2<-unlist(sapply(null_ak2,networklevel,index='H2'))
net_spec_null_ak3<-unlist(sapply(null_ak3,networklevel,index='H2'))
net_spec_null_ak4<-unlist(sapply(null_ak4,networklevel,index='H2'))
net_spec_null_ak5<-unlist(sapply(null_ak5,networklevel,index='H2'))
net_spec_null_ak6<-unlist(sapply(null_ak6,networklevel,index='H2'))
net_spec_null_ro1<-unlist(sapply(null_ro1,networklevel,index='H2'))
net_spec_null_ro2<-unlist(sapply(null_ro2,networklevel,index='H2'))
net_spec_null_ro3<-unlist(sapply(null_ro3,networklevel,index='H2'))
net_spec_null_ro4<-unlist(sapply(null_ro4,networklevel,index='H2'))
net_spec_null_ro5<-unlist(sapply(null_ro5,networklevel,index='H2'))
net_spec_null_ro6<-unlist(sapply(null_ro6,networklevel,index='H2'))

net_spec_null<-list(c(AK1=net_spec_null_ak1,AK2=net_spec_null_ak2,AK3=net_spec_null_ak3,AK4=net_spec_null_ak4,
                     AK5=net_spec_null_ak5,AK6=net_spec_null_ak6,RO1=net_spec_null_ro1,RO2=net_spec_null_ro2,
                     RO3=net_spec_null_ro3,RO4=net_spec_null_ro4,RO5=net_spec_null_ro5,RO6=net_spec_null_ro6))

net_spec_null_df <- as.data.frame(net_spec_null)
write.csv(as_tibble(net_spec_null_df),file="net_spec_null_df.csv")

net_spec_null_df<-cbind(net_spec_null_df,null_hab_type_labs,null_plot_labs)
colnames(net_spec_null_df) <- c("networkspecialization","HabitatType","plot")
net_spec_null_sum<- summarySE(net_spec_null_df, measurevar="networkspecialization", 
                            groupvars=c("HabitatType"))

H2_null_plot <- ggplot(data=net_spec_null_sum,aes(x=HabitatType,y=networkspecialization)) +
  geom_errorbar(aes(ymin=networkspecialization-sd,ymax=networkspecialization+sd),width=0.1,colour="black",size=2) +
  geom_point(pch=21, colour="black",size=5, fill="white",stroke=2) +
  stat_summary(data=H2_dat,aes(x=HabitatType,y=networkspecialization,stat="identity"),
               fun.data=mean_sdl, fun.args=list(mult=1),geom="errorbar",color="black",width=0.1,size=2) +
  stat_summary(data=H2_dat,aes(x=HabitatType,y=networkspecialization,stat="identity"),
               fun.y=mean,geom="point",colour="black",size=5) +
  ylab("Network Specialization (H2')") + xlab("Habitat Type") +
  scale_y_continuous(breaks=seq(0,1,by=0.1),limits=c(0,0.5)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) 

plot(H2_null_plot)

#t-tests comparing observed to null mean H2
ak.H2.null.t <- t.test(subset(H2_dat,HabitatType == "Afforested koa")$networkspecialization,
                           subset(net_spec_null_df,HabitatType == "Afforested koa")$networkspecialization,
                           paired=FALSE, var.equal=FALSE)
ak.H2.null.t

ro.H2.null.t <- t.test(subset(H2_dat,HabitatType == "Remnant ohia")$networkspecialization,
                           subset(net_spec_null_df,HabitatType == "Remnant ohia")$networkspecialization,
                           paired=FALSE, var.equal=FALSE)
ro.H2.null.t

##########################################
###   Host symbiont range ("degree")   ###
##########################################
host_degree <- lapply(haka_bipartite,specieslevel,index='degree', level= "higher")

host_degree_dat<-data.frame(unlist(host_degree))
host_ak_labs<-as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
                                  "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),times=6))
host_ro_labs_1<-as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                                    "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),times=5))
host_ro_labs_2<-as.data.frame(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
                                "Cheirodendron trigynum","Rubus hawaiiensis","Grass"))
names(host_ak_labs)<-names(host_ro_labs_1)
names(host_ro_labs_1)<-names(host_ro_labs_2)
host_labs <- rbind(host_ak_labs,host_ro_labs_1,host_ro_labs_2)

hab_ak_labs<-as.data.frame(rep("Afforested koa",times=42))
hab_ro_labs<-as.data.frame(rep("Remnant ohia",times=37))
names(hab_ak_labs)<-names(hab_ro_labs)
hab_labs<-rbind(hab_ak_labs,hab_ro_labs)

host_degree_df <-cbind(host_degree_dat,host_labs,hab_labs)
colnames(host_degree_df) <- c("degree","Host","HabitatType")
host_degree_df$Host<- factor(host_degree_df$Host,levels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                                                      "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                           labels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))


host_degree_p <- ggplot(host_degree_df,aes(x=HabitatType,y=degree,fill=Host)) + 
  geom_boxplot(size=0.75, outlier.alpha = 0) +
  geom_point(pch=21,position=position_jitterdodge(),size=3) +
  ylab("Host Symbiont Range (no. of links)") + xlab("Habitat Type") +
  scale_fill_manual(values=c("#D73027","#FC8D59","#FEE090","#008000","#E0F3F8","#008B8B","#4575B4")) +
  scale_y_continuous(breaks=seq(20,80,by=10),limits=c(20,80)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=30,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.text = element_text(face="italic",size=15)) +
  ggtitle("A)") 

plot(host_degree_p)
ggsave("Fig_6A_host_symbiont_range.tiff", 
       plot = host_degree_p,width=13, height=8)

#Degree among hosts and habitat type Stats
degree_glm1 <- glm(degree ~ Host*HabitatType, data=host_degree_df, family = quasipoisson(link="log")) 
par(mfrow = c(2, 2))
plot(degree_glm1)
degree_chi.anova<-anova(degree_glm1,test="Chisq")
degree_chi.anova
#Do pairwise contrasts using estimated margnial means (EMM)
#Least-squares means are discussed, and the term ``estimated marginal means''
#is suggested, in Searle, Speed, and Milliken (1980) Population marginal means
#in the linear model: An alternative to least squares means, The American
#Statistician 34(4), 216-221
library(emmeans);library(multcompView)
#Create an interaction column in the dataframe that accounts for each unique host x habitat type combination
host_degree_df$interaction<- interaction(host_degree_df$HabitatType,host_degree_df$Host)
#Run glm again on the interaction term
degree_glm2 <- glm(degree ~ interaction, data=host_degree_df, family = quasipoisson(link="log"))
#Calculate EMM on interactions
degree_emmeans <- emmeans(degree_glm2, specs="interaction")
degree_emmeans
degree_posthoc.pairs = pairs(degree_emmeans)
degree_posthoc.pairs
#Create letters for interaction differences
degree_mc_letters<-cld(degree_emmeans,Letters="abcdefg")
degree_mc_letters

#Generate null values for host specialization and save as a list
degree_null_AK1<-unlist(sapply(null_ak1,specieslevel,level='higher',index='degree'))
degree_null_AK2<-unlist(sapply(null_ak2,specieslevel,level='higher',index='degree'))
degree_null_AK3<-unlist(sapply(null_ak3,specieslevel,level='higher',index='degree'))
degree_null_AK4<-unlist(sapply(null_ak4,specieslevel,level='higher',index='degree'))
degree_null_AK5<-unlist(sapply(null_ak5,specieslevel,level='higher',index='degree'))
degree_null_AK6<-unlist(sapply(null_ak6,specieslevel,level='higher',index='degree'))
degree_null_RO1<-unlist(sapply(null_ro1,specieslevel,level='higher',index='degree'))
degree_null_RO2<-unlist(sapply(null_ro2,specieslevel,level='higher',index='degree'))
degree_null_RO3<-unlist(sapply(null_ro3,specieslevel,level='higher',index='degree'))
degree_null_RO4<-unlist(sapply(null_ro4,specieslevel,level='higher',index='degree'))
degree_null_RO5<-unlist(sapply(null_ro5,specieslevel,level='higher',index='degree'))
degree_null_RO6<-unlist(sapply(null_ro6,specieslevel,level='higher',index='degree'))


degree_null <- list(c(degree_null_AK1, degree_null_AK2, degree_null_AK3, degree_null_AK4, degree_null_AK5, degree_null_AK6,
                    degree_null_RO1, degree_null_RO2, degree_null_RO3, degree_null_RO4, degree_null_RO5, degree_null_RO6))

#Create null network specialization dataframe
degree_null_df<- as.data.frame(degree_null)
write.csv(as_tibble(degree_null_df),file="degree_null_df.csv")
#Create df of habitat labels
ak_labs <- as.data.frame(rep("Afforested koa",times=42000))
ro_labs <-as.data.frame(rep("Remnant ohia",times=37000))
names(ak_labs) <- names (ro_labs)
degree_hab_labs<-rbind(ak_labs, ro_labs)
#Create df of host labels
ak_host_labs <- as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa",
                                    "Myrsine lessertiana","Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                                  each=1000,times=6))
ro_host_labs_1 <- as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                                      "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),each=1000,times=5))
ro_host_labs_2 <- as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa",
                                      "Myrsine lessertiana","Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                                    each=1000))
names(ak_host_labs)<-names(ro_host_labs_1)
names(ro_host_labs_1)<-names(ro_host_labs_2)
degree_host_labs<-rbind(ak_host_labs, ro_host_labs_1, ro_host_labs_2)
#Create df of plot labels
ak_plot_labs <- as.data.frame(rep(c("AK1","AK2","AK3","AK4","AK5","AK6"),each=7000))
ro_plot_labs_1<- as.data.frame(rep(c("RO1","RO2","RO3","RO4","RO5"),each=6000))
ro_plot_labs_2<- as.data.frame(rep("RO6",times=7000))
names(ak_plot_labs)<-names(ro_plot_labs_1)
names(ro_plot_labs_1)<-names(ro_plot_labs_2)
degree_plot_labs<-rbind(ak_plot_labs,ro_plot_labs_1,ro_plot_labs_2)

#Add metadata label columns to null dataframe
degree_null_df<-cbind(degree_null_df,degree_hab_labs,degree_host_labs,degree_plot_labs)
colnames(degree_null_df) <- c("degree","HabitatType","Host","Plot")
degree_null_df$Host<- factor(degree_null_df$Host,levels=c("Metrosideros polymorpha", "Acacia koa", 
                                                          "Coprosma rhynchocarpa", "Myrsine lessertiana", 
                                                          "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                           labels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Create summary of null values for the plot
degree_null_sum<- summarySE(degree_null_df, measurevar="degree", 
                          groupvars=c("HabitatType","Host"))

#Get colour codes for the 7 hosts used in the observed plot
degree_palette = get_palette("RdYlBu",7)
degree_palette

#Null plot
degree_plot <- ggplot(degree_null_sum,aes(x=HabitatType, y=degree, fill=Host),inherit.aes=FALSE) +
  geom_point(size=6,position=position_jitterdodge()) +
  geom_errorbar(aes(ymin=degree-sd,ymax=degree+sd),width=0.1,size=1, position=position_jitterdodge()) +
  geom_boxplot(data=host_degree_df, aes(x=HabitatType, y=degree, fill=Host),size=0.75,outlier.alpha=0) +
  geom_point(data=host_degree_df, aes(x=HabitatType, y=degree, fill=Host),pch=21,position=position_jitterdodge(),size=3) +
  scale_fill_brewer(palette='RdYlBu') +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=30,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.text = element_text(face="italic")) 

plot(degree_plot)

#t-tests comparing observed to null mean
#Afforested koa
ak.degree.Metrosideros.t.test <- t.test(subset(host_degree_df,HabitatType == "Afforested koa" & 
                                          Host=="Metrosideros polymorpha")$degree,
                                 subset(degree_null_df,HabitatType == "Afforested koa" & 
                                          Host=="Metrosideros polymorpha")$degree,
                                 paired=FALSE, var.equal = FALSE)
ak.degree.Metrosideros.t.test


ak.degree.Acacia.t.test <- t.test(subset(host_degree_df,HabitatType == "Afforested koa" & 
                                    Host=="Acacia koa")$degree,
                           subset(degree_null_df,HabitatType == "Afforested koa" & 
                                    Host=="Acacia koa")$degree,
                           paired=FALSE, var.equal = FALSE)
ak.degree.Acacia.t.test

ak.degree.Coprosma.t.test <- t.test(subset(host_degree_df,HabitatType == "Afforested koa" & 
                                      Host=="Coprosma rhynchocarpa")$degree,
                             subset(degree_null_df,HabitatType == "Afforested koa" & 
                                      Host=="Coprosma rhynchocarpa")$degree,
                             paired=FALSE,var.equal = FALSE)
ak.degree.Coprosma.t.test

ak.degree.Myrsine.t.test <- t.test(subset(host_degree_df,HabitatType == "Afforested koa" & 
                                     Host=="Myrsine lessertiana")$degree,
                            subset(degree_null_df,HabitatType == "Afforested koa" & 
                                     Host=="Myrsine lessertiana")$degree,
                            paired=FALSE,var.equal = FALSE)
ak.degree.Myrsine.t.test

ak.degree.Cheirodendron.t.test <- t.test(subset(host_degree_df,HabitatType == "Afforested koa" & 
                                           Host=="Cheirodendron trigynum")$degree,
                                  subset(degree_null_df,HabitatType == "Afforested koa" & 
                                           Host=="Cheirodendron trigynum")$degree,
                                  paired=FALSE,var.equal = FALSE)
ak.degree.Cheirodendron.t.test

ak.degree.Rubus.t.test <- t.test(subset(host_degree_df,HabitatType == "Afforested koa" & 
                                   Host=="Rubus hawaiiensis")$degree,
                          subset(degree_null_df,HabitatType == "Afforested koa" & 
                                   Host=="Rubus hawaiiensis")$degree,
                          paired=FALSE,var.equal=FALSE)
ak.degree.Rubus.t.test

ak.degree.Grass.t.test <- t.test(subset(host_degree_df,HabitatType == "Afforested koa" & 
                                   Host=="Grass")$degree,
                          subset(degree_null_df,HabitatType == "Afforested koa" & 
                                   Host=="Grass")$degree,
                          paired=FALSE,var.equal = FALSE)
ak.degree.Grass.t.test

#Remnant ohia
ro.degree.Metrosideros.t.test <- t.test(subset(host_degree_df,HabitatType == "Remnant ohia" & 
                                          Host=="Metrosideros polymorpha")$degree,
                                 subset(degree_null_df,HabitatType == "Remnant ohia" & 
                                          Host=="Metrosideros polymorpha")$degree,
                                 paired=FALSE,var.equal = FALSE)
ro.degree.Metrosideros.t.test

ro.degree.Acacia.t.test <- t.test(subset(host_degree_df,HabitatType == "Remnant ohia" & 
                                    Host=="Acacia koa")$degree,
                           subset(degree_null_df,HabitatType == "Remnant ohia" & 
                                    Host=="Acacia koa")$degree,
                           paired=FALSE,var.equal = FALSE)
ro.degree.Acacia.t.test

#*** Only one Coprosma found in our RO plots ****#
ro.degree.Coprosma.t.test <- t.test(subset(degree_null_df,HabitatType == "Remnant ohia" & 
                                      Host=="Coprosma rhynchocarpa")$degree,
                             mu=subset(host_degree_df,HabitatType == "Remnant ohia" & 
                                         Host=="Coprosma rhynchocarpa")$degree)
ro.degree.Coprosma.t.test
#*** Only one Coprosma found in our RO plots ****#

ro.degree.Myrsine.t.test <- t.test(subset(host_degree_df,HabitatType == "Remnant ohia" & 
                                     Host=="Myrsine lessertiana")$degree,
                            subset(degree_null_df,HabitatType == "Remnant ohia" & 
                                     Host=="Myrsine lessertiana")$degree,
                            paired=FALSE,var.equal = FALSE)
ro.degree.Myrsine.t.test

ro.degree.Cheirodendron.t.test <- t.test(subset(host_degree_df,HabitatType == "Remnant ohia" & 
                                           Host=="Cheirodendron trigynum")$degree,
                                  subset(degree_null_df,HabitatType == "Remnant ohia" & 
                                           Host=="Cheirodendron trigynum")$degree,
                                  paired=FALSE, var.equal = FALSE)
ro.degree.Cheirodendron.t.test

ro.degree.Rubus.t.test <- t.test(subset(host_degree_df,HabitatType == "Remnant ohia" & 
                                   Host=="Rubus hawaiiensis")$degree,
                          subset(degree_null_df,HabitatType == "Remnant ohia" & 
                                   Host=="Rubus hawaiiensis")$degree,
                          paired=FALSE,var.equal = FALSE)
ro.degree.Rubus.t.test

ro.degree.Grass.t.test <- t.test(subset(host_degree_df,HabitatType == "Remnant ohia" & 
                                   Host=="Grass")$degree,
                          subset(degree_null_df,HabitatType == "Remnant ohia" & 
                                   Host=="Grass")$degree,
                          paired=FALSE,var.equal = FALSE)
ro.degree.Grass.t.test
#########################################
###   Host specialization on AMF (d') ###
#########################################
host_spec_ak1 <- dfun(t(haka_bipartite$AK1))
host_spec_ak2 <- dfun(t(haka_bipartite$AK2))
host_spec_ak3 <- dfun(t(haka_bipartite$AK3))
host_spec_ak4 <- dfun(t(haka_bipartite$AK4))
host_spec_ak5 <- dfun(t(haka_bipartite$AK5))
host_spec_ak6 <- dfun(t(haka_bipartite$AK6))
host_spec_ro1 <- dfun(t(haka_bipartite$RO1))
host_spec_ro2 <- dfun(t(haka_bipartite$RO2))
host_spec_ro3 <- dfun(t(haka_bipartite$RO3))
host_spec_ro4 <- dfun(t(haka_bipartite$RO4))
host_spec_ro5 <- dfun(t(haka_bipartite$RO5))
host_spec_ro6 <- dfun(t(haka_bipartite$RO6))

host_spec_list <-list(c(AK1=host_spec_ak1$dprime,AK2= host_spec_ak2$dprime,AK3=host_spec_ak3$dprime,
                        AK4=host_spec_ak4$dprime,AK5=host_spec_ak5$dprime,AK6=host_spec_ak6$dprime,
                        RO1=host_spec_ro1$dprime,RO2=host_spec_ro2$dprime,RO3=host_spec_ro3$dprime,
                        RO4=host_spec_ro4$dprime,RO5=host_spec_ro5$dprime,RO6=host_spec_ro6$dprime))
host_spec_df <- as.data.frame(host_spec_list)


host_ak_labs<-as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
                                  "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),times=6))
host_ro_labs_1<-as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                                    "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),times=5))
host_ro_labs_2<-as.data.frame(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa","Myrsine lessertiana",
                              "Cheirodendron trigynum","Rubus hawaiiensis","Grass"))

names(host_ak_labs)<-names(host_ro_labs_1)
names(host_ak_labs)<-names(host_ro_labs_2)

host_labs <- rbind(host_ak_labs,host_ro_labs_1)
names(host_labs)<-names(host_ro_labs_2)
host_labs<-rbind(host_labs,host_ro_labs_2)

hab_ak_labs<-as.data.frame(rep("Afforested koa",times=42))
hab_ro_labs<-as.data.frame(rep("Remnant ohia",times=37))
names(hab_ak_labs)<-names(hab_ro_labs)
hab_labs<-rbind(hab_ak_labs,hab_ro_labs)

host_spec_df <-cbind(host_spec_df,host_labs,hab_labs)
colnames(host_spec_df) <- c("dprime","Host","HabitatType")
host_spec_df$Host<- factor(host_spec_df$Host,levels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                                          "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"),
                     labels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                              "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))


host_dfun_p <- ggplot(host_spec_df,aes(x=HabitatType,y=dprime,fill=Host)) + 
  geom_boxplot(size=0.75, outlier.alpha = 0) +
  geom_point(pch=21,position=position_jitterdodge(),size=3) +
  ylab("Host Specialization on AM fungi (d')") + xlab("Habitat Type") +
  ylim(0,0.7) +
  scale_fill_manual(values=c("#D73027","#FC8D59","#FEE090","#008000","#E0F3F8","#008B8B","#4575B4")) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=30,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.text = element_text(face="italic",size=15)) +
  ggtitle("B)") 

plot(host_dfun_p)
ggsave("Fig_6B_host_specialization.tiff", 
       plot = host_dfun_p, width=13, height=8)

#Specialization among hosts and habitat type Stats
host_spec_glm1 <- glm(dprime ~ HabitatType*Host, data=host_spec_df, family = quasipoisson(link="log")) 
par(mfrow = c(2, 2))
plot(host_spec_glm1)
host_spec_chi.anova<-anova(host_spec_glm1,test="Chisq")
host_spec_chi.anova
#Posthoc analysis (same as for Host symbiont range)
host_spec_df$interaction<- interaction(host_spec_df$HabitatType,host_spec_df$Host)
dprime_glm2<-glm(dprime ~ interaction, data=host_spec_df, family = quasipoisson(link="log")) 
dprime_emmeans <- emmeans(dprime_glm2, specs="interaction")
dprime_emmeans
dprime_posthoc.pairs = pairs(dprime_emmeans)
dprime_posthoc.pairs
#Create letters for interaction differences
dprime_mc_letters<-cld(dprime_emmeans,Letters="abcdefg")
dprime_mc_letters

### d' null model analysis

#Generate null values for host specialization and save as a list
dfun_null_AK1<-unlist(sapply(null_ak1,specieslevel,level='higher',index='d'))
dfun_null_AK2<-unlist(sapply(null_ak2,specieslevel,level='higher',index='d'))
dfun_null_AK3<-unlist(sapply(null_ak3,specieslevel,level='higher',index='d'))
dfun_null_AK4<-unlist(sapply(null_ak4,specieslevel,level='higher',index='d'))
dfun_null_AK5<-unlist(sapply(null_ak5,specieslevel,level='higher',index='d'))
dfun_null_AK6<-unlist(sapply(null_ak6,specieslevel,level='higher',index='d'))
dfun_null_RO1<-unlist(sapply(null_ro1,specieslevel,level='higher',index='d'))
dfun_null_RO2<-unlist(sapply(null_ro2,specieslevel,level='higher',index='d'))
dfun_null_RO3<-unlist(sapply(null_ro3,specieslevel,level='higher',index='d'))
dfun_null_RO4<-unlist(sapply(null_ro4,specieslevel,level='higher',index='d'))
dfun_null_RO5<-unlist(sapply(null_ro5,specieslevel,level='higher',index='d'))
dfun_null_RO6<-unlist(sapply(null_ro6,specieslevel,level='higher',index='d'))


dfun_null <- list(c(dfun_null_AK1, dfun_null_AK2, dfun_null_AK3, dfun_null_AK4, dfun_null_AK5, dfun_null_AK6,
                    dfun_null_RO1, dfun_null_RO2, dfun_null_RO3, dfun_null_RO4, dfun_null_RO5, dfun_null_RO6))

#Create null network specialization dataframe
dfun_null_df<- as.data.frame(dfun_null)
write.csv(as_tibble(dfun_null_df),file="dfun_null_df.csv")
#Create df of habitat labels
ak_labs <- as.data.frame(rep("Afforested koa",times=42000))
ro_labs <-as.data.frame(rep("Remnant ohia",times=37000))
names(ak_labs) <- names (ro_labs)
dfun_hab_labs<-rbind(ak_labs, ro_labs)
#Create df of host labels
ak_host_labs <- as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa",
                                    "Myrsine lessertiana","Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                                  each=1000,times=6))
ro_host_labs_1 <- as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Myrsine lessertiana",
                                      "Cheirodendron trigynum","Rubus hawaiiensis","Grass"),each=1000,times=5))
ro_host_labs_2 <- as.data.frame(rep(c("Metrosideros polymorpha","Acacia koa","Coprosma rhynchocarpa",
                                      "Myrsine lessertiana","Cheirodendron trigynum","Rubus hawaiiensis","Grass"),
                                    each=1000))
names(ak_host_labs)<-names(ro_host_labs_1)
names(ak_host_labs)<-names(ro_host_labs_2)
dfun_host_labs<-rbind(ak_host_labs, ro_host_labs_1, ro_host_labs_2)
#Create df of plot labels
ak_plot_labs <- as.data.frame(rep(c("AK1","AK2","AK3","AK4","AK5","AK6"),each=7000))
ro_plot_labs_1<- as.data.frame(rep(c("RO1","RO2","RO3","RO4","RO5"),each=6000))
ro_plot_labs_2<- as.data.frame(rep("RO6",times=7000))
names(ak_plot_labs)<-names(ro_plot_labs_1)
names(ro_plot_labs_1)<-names(ro_plot_labs_2)
dfun_plot_labs<-rbind(ak_plot_labs,ro_plot_labs_1,ro_plot_labs_2)

#Add metadata label columns to null dataframe
dfun_null_df<-cbind(dfun_null_df,dfun_hab_labs,dfun_host_labs,dfun_plot_labs)
colnames(dfun_null_df) <- c("dprime","HabitatType","Host","Plot")
dfun_null_df$Host<- factor(dfun_null_df$Host,levels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                                                      "Myrsine lessertiana", "Cheirodendron trigynum", 
                                                      "Rubus hawaiiensis", "Grass"),
                           labels=c("Metrosideros polymorpha", "Acacia koa", "Coprosma rhynchocarpa", 
                                    "Myrsine lessertiana", "Cheirodendron trigynum", "Rubus hawaiiensis", "Grass"))

#Create summary of null values for the plot
dfun_null_sum<- summarySE(dfun_null_df, measurevar="dprime", 
                          groupvars=c("HabitatType","Host"))

#Get colour codes for the 7 hosts used in the observed plot
dfun_palette = get_palette("RdYlBu",7)
dfun_palette

#Null plot
dfun_plot <- ggplot(dfun_null_sum,aes(x=HabitatType, y=dprime, fill=Host),inherit.aes=FALSE) +
  geom_point(size=6,position=position_jitterdodge()) +
  geom_errorbar(aes(ymin=dprime-sd,ymax=dprime+sd),width=0.1,size=1, position=position_jitterdodge()) +
  geom_boxplot(data=host_spec_df, aes(x=HabitatType, y=dprime, fill=Host),size=0.75,outlier.alpha=0) +
  geom_point(data=host_spec_df, aes(x=HabitatType, y=dprime, fill=Host),pch=21,position=position_jitterdodge(),size=3) +
  scale_fill_brewer(palette='RdYlBu') +
  scale_y_continuous(breaks=seq(0,0.8,by=0.2),limits=c(0,0.8)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(angle=30,hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  theme(legend.text = element_text(face="italic")) 

plot(dfun_plot)

#d' t-tests comparing observed to null mean
#Afforested koa
ak.Metrosideros.t.test <- t.test(subset(host_spec_df,HabitatType == "Afforested koa" & 
                                               Host=="Metrosideros polymorpha")$dprime,
                                 subset(dfun_null_df,HabitatType == "Afforested koa" & 
                                                  Host=="Metrosideros polymorpha")$dprime,
                                 paired=FALSE, var.equal=FALSE)
ak.Metrosideros.t.test


ak.Acacia.t.test <- t.test(subset(host_spec_df,HabitatType == "Afforested koa" & 
                                          Host=="Acacia koa")$dprime,
                                 subset(dfun_null_df,HabitatType == "Afforested koa" & 
                                                  Host=="Acacia koa")$dprime,
                                 paired=FALSE, var.equal=FALSE)
ak.Acacia.t.test

ak.Coprosma.t.test <- t.test(subset(host_spec_df,HabitatType == "Afforested koa" & 
                                          Host=="Coprosma rhynchocarpa")$dprime,
                                 subset(dfun_null_df,HabitatType == "Afforested koa" & 
                                                  Host=="Coprosma rhynchocarpa")$dprime,
                                 paired=FALSE,var.equal=FALSE)
ak.Coprosma.t.test

ak.Myrsine.t.test <- t.test(subset(host_spec_df,HabitatType == "Afforested koa" & 
                                          Host=="Myrsine lessertiana")$dprime,
                                 subset(dfun_null_df,HabitatType == "Afforested koa" & 
                                                  Host=="Myrsine lessertiana")$dprime,
                                 paired=FALSE,var.equal=FALSE)
ak.Myrsine.t.test

ak.Cheirodendron.t.test <- t.test(subset(host_spec_df,HabitatType == "Afforested koa" & 
                                          Host=="Cheirodendron trigynum")$dprime,
                                 subset(dfun_null_df,HabitatType == "Afforested koa" & 
                                                  Host=="Cheirodendron trigynum")$dprime,
                                 paired=FALSE,var.equal=FALSE)
ak.Cheirodendron.t.test

ak.Rubus.t.test <- t.test(subset(host_spec_df,HabitatType == "Afforested koa" & 
                                          Host=="Rubus hawaiiensis")$dprime,
                                 subset(dfun_null_df,HabitatType == "Afforested koa" & 
                                                  Host=="Rubus hawaiiensis")$dprime,
                                 paired=FALSE,var.equal=FALSE)
ak.Rubus.t.test

ak.Grass.t.test <- t.test(subset(host_spec_df,HabitatType == "Afforested koa" & 
                                          Host=="Grass")$dprime,
                                 subset(dfun_null_df,HabitatType == "Afforested koa" & 
                                                  Host=="Grass")$dprime,
                                 paired=FALSE,var.equal=FALSE)
ak.Grass.t.test

#Remnant ohia
ro.Metrosideros.t.test <- t.test(subset(host_spec_df,HabitatType == "Remnant ohia" & 
                                          Host=="Metrosideros polymorpha")$dprime,
                                 subset(dfun_null_df,HabitatType == "Remnant ohia" & 
                                                  Host=="Metrosideros polymorpha")$dprime,
                                 paired=FALSE,var.equal=FALSE)
ro.Metrosideros.t.test

ro.Acacia.t.test <- t.test(subset(host_spec_df,HabitatType == "Remnant ohia" & 
                                    Host=="Acacia koa")$dprime,
                           subset(dfun_null_df,HabitatType == "Remnant ohia" & 
                                            Host=="Acacia koa")$dprime,
                           paired=FALSE,var.equal=FALSE)
ro.Acacia.t.test

                                #*** Only one Coprosma found in our RO plots ****#
ro.Coprosma.t.test <- t.test(subset(dfun_null_df,HabitatType == "Remnant ohia" & 
                                              Host=="Coprosma rhynchocarpa")$dprime,
                             mu=subset(host_spec_df,HabitatType == "Remnant ohia" & 
                                         Host=="Coprosma rhynchocarpa")$dprime)
ro.Coprosma.t.test
                                .#*** Only one Coprosma found in our RO plots ****#

ro.Myrsine.t.test <- t.test(subset(host_spec_df,HabitatType == "Remnant ohia" & 
                                     Host=="Myrsine lessertiana")$dprime,
                            subset(dfun_null_df,HabitatType == "Remnant ohia" & 
                                             Host=="Myrsine lessertiana")$dprime,
                            paired=FALSE,var.equal=FALSE)
ro.Myrsine.t.test

ro.Cheirodendron.t.test <- t.test(subset(host_spec_df,HabitatType == "Remnant ohia" & 
                                           Host=="Cheirodendron trigynum")$dprime,
                                  subset(dfun_null_df,HabitatType == "Remnant ohia" & 
                                                   Host=="Cheirodendron trigynum")$dprime,
                                  paired=FALSE, var.equal=FALSE)
ro.Cheirodendron.t.test

ro.Rubus.t.test <- t.test(subset(host_spec_df,HabitatType == "Remnant ohia" & 
                                   Host=="Rubus hawaiiensis")$dprime,
                          subset(dfun_null_df,HabitatType == "Remnant ohia" & 
                                           Host=="Rubus hawaiiensis")$dprime,
                          paired=FALSE,var.equal=FALSE)
ro.Rubus.t.test

ro.Grass.t.test <- t.test(subset(host_spec_df,HabitatType == "Remnant ohia" & 
                                   Host=="Grass")$dprime,
                          subset(dfun_null_df,HabitatType == "Remnant ohia" & 
                                           Host=="Grass")$dprime,
                          paired=FALSE,var.equal=FALSE)
ro.Grass.t.test

#########################################
###   Observed bipartite metric plots ###
#########################################
require(gridExtra)
Fig_1<- grid.arrange(H2_p, link_dens_plot, plant_links,
                              ncol=3,nrow=1)
ggsave("Fig_1_observed_bipartite_metrics.tiff", 
       plot = Fig_1, width = 15, height = 10)

obs_plot2<-grid.arrange(host_degree_p, host_dfun_p, ncol=1, nrow=2)

ggsave("Host_degree_and_dprime.tiff",
       plot= obs_plot2, width = 12, height = 20)


#####################################  
#### bipartite networks by plot #####
#####################################
library(igraph);library(qgraph)

                            ###Afforested koa sites###
###Example using AK1
#create igraph object from incidence matrix
ak1_bi <- graph.incidence(haka_bipartite$AK1, weighted = T)
#Identify that only want hosts labelled
V(ak1_bi)$label<-V(ak1_bi)$name
#Change border colour of vertices to white to make graph better looking
V(ak1_bi)$frame.color = "black"
#Color by node type. First colour is AMF, second is plant
col<-c("#F1F3CE","#265C00")
#Change edge widths to reflect the weight of interaction
E(ak1_bi)$width<-E(ak1_bi)$weight/500
#Change the size of vertices to reflect tings like degree of centrality
deg = centr_degree(ak1_bi, mode="all")
V(ak1_bi)$size = sqrt(deg$res)
#Create custom layout to better visualize nodes associating with only one host species
ak1_cust_layout<-layout.reingold.tilford(ak1_bi,circular=T)
#plot
set.seed(605)
plot(ak1_bi, vertex.size=6, vertex.color=col[as.numeric(V(ak1_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=ak1_cust_layout,
     vertex.label=ifelse(degree(ak1_bi)>40,V(ak1_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#Colour by Family. First need to assign colours to match previously made graphs
net_tax=read.csv("haka_tax.csv", header = TRUE,row.names = 1)
V(ak1_bi)$Family=as.character(net_tax$Family[match(V(ak1_bi)$name,net_tax$Species)])
V(ak1_bi)$Family
V(ak1_bi)$color=V(ak1_bi)$Family
show_col(hue_pal()(8))
V(ak1_bi)$color=gsub("Acaulosporaceae","#F8766D",V(ak1_bi)$color)
V(ak1_bi)$color=gsub("Archaeosporaceae","#CD9600",V(ak1_bi)$color)
V(ak1_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(ak1_bi)$color)
V(ak1_bi)$color=gsub("Diversisporaceae","#00BE67",V(ak1_bi)$color)
V(ak1_bi)$color=gsub("Gigasporaceae","#00BFC4",V(ak1_bi)$color)
V(ak1_bi)$color=gsub("Glomeraceae","#00A9FF",V(ak1_bi)$color)
V(ak1_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(ak1_bi)$color)
V(ak1_bi)$color=gsub("unassigned","#FF61CC",V(ak1_bi)$color)
set.seed(605)
plot(ak1_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=ak1_cust_layout,
     vertex.label=ifelse(degree(ak1_bi)>40,V(ak1_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#Legend to be used with networks coloured by AMF family
legend('topright',c("Acaulosporaceae","Archaeosporaceae","Claroideoglomeraceae","Diversisporaceae","Gigasporaceae",
         "Glomeraceae","Paraglomeraceae","unassigned"),pch=21,
       pt.bg=c("#F8766D","#CD9600","#7CAE00","#00BE67","#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
       pt.cex=4,cex=2)

###AK2
AK2_bi <- graph.incidence(haka_bipartite$AK2, weighted = T)
V(AK2_bi)$label<-V(AK2_bi)$name
V(AK2_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(AK2_bi)$width<-E(AK2_bi)$weight/500
deg = centr_degree(AK2_bi, mode="all")
AK2_cust_layout<-layout.reingold.tilford(AK2_bi,circular=T)
set.seed(605)
plot(AK2_bi, vertex.size=6, vertex.color=col[as.numeric(V(AK2_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=AK2_cust_layout,
     vertex.label=ifelse(degree(AK2_bi)>50,V(AK2_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#AK2 coloured by Family. 
V(AK2_bi)$Family=as.character(net_tax$Family[match(V(AK2_bi)$name,net_tax$Species)])
V(AK2_bi)$Family
V(AK2_bi)$color=V(AK2_bi)$Family
V(AK2_bi)$color=gsub("Acaulosporaceae","#F8766D",V(AK2_bi)$color)
V(AK2_bi)$color=gsub("Archaeosporaceae","#CD9600",V(AK2_bi)$color)
V(AK2_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(AK2_bi)$color)
V(AK2_bi)$color=gsub("Diversisporaceae","#00BE67",V(AK2_bi)$color)
V(AK2_bi)$color=gsub("Gigasporaceae","#00BFC4",V(AK2_bi)$color)
V(AK2_bi)$color=gsub("Glomeraceae","#00A9FF",V(AK2_bi)$color)
V(AK2_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(AK2_bi)$color)
V(AK2_bi)$color=gsub("unassigned","#FF61CC",V(AK2_bi)$color)
set.seed(605)
plot(AK2_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=AK2_cust_layout,
     vertex.label=ifelse(degree(AK2_bi)>50,V(AK2_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###AK3
AK3_bi <- graph.incidence(haka_bipartite$AK3, weighted = T)
V(AK3_bi)$label<-V(AK3_bi)$name
V(AK3_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(AK3_bi)$width<-E(AK3_bi)$weight/500
deg = centr_degree(AK3_bi, mode="all")
AK3_cust_layout<-layout.reingold.tilford(AK3_bi,circular=T)
set.seed(605)
plot(AK3_bi, vertex.size=6, vertex.color=col[as.numeric(V(AK3_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=AK3_cust_layout,
     vertex.label=ifelse(degree(AK3_bi)>50,V(AK3_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#AK3 coloured by Family. 
V(AK3_bi)$Family=as.character(net_tax$Family[match(V(AK3_bi)$name,net_tax$Species)])
V(AK3_bi)$Family
V(AK3_bi)$color=V(AK3_bi)$Family
V(AK3_bi)$color=gsub("Acaulosporaceae","#F8766D",V(AK3_bi)$color)
V(AK3_bi)$color=gsub("Archaeosporaceae","#CD9600",V(AK3_bi)$color)
V(AK3_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(AK3_bi)$color)
V(AK3_bi)$color=gsub("Diversisporaceae","#00BE67",V(AK3_bi)$color)
V(AK3_bi)$color=gsub("Gigasporaceae","#00BFC4",V(AK3_bi)$color)
V(AK3_bi)$color=gsub("Glomeraceae","#00A9FF",V(AK3_bi)$color)
V(AK3_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(AK3_bi)$color)
V(AK3_bi)$color=gsub("unassigned","#FF61CC",V(AK3_bi)$color)
set.seed(605)
plot(AK3_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=AK3_cust_layout,
     vertex.label=ifelse(degree(AK3_bi)>50,V(AK3_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###AK4
AK4_bi <- graph.incidence(haka_bipartite$AK4, weighted = T)
V(AK4_bi)$label<-V(AK4_bi)$name
V(AK4_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(AK4_bi)$width<-E(AK4_bi)$weight/500
deg = centr_degree(AK4_bi, mode="all")
AK4_cust_layout<-layout.reingold.tilford(AK4_bi,circular=T)
set.seed(605)
plot(AK4_bi, vertex.size=6, vertex.color=col[as.numeric(V(AK4_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=AK4_cust_layout,
     vertex.label=ifelse(degree(AK4_bi)>50,V(AK4_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#AK4 coloured by Family. 
V(AK4_bi)$Family=as.character(net_tax$Family[match(V(AK4_bi)$name,net_tax$Species)])
V(AK4_bi)$Family
V(AK4_bi)$color=V(AK4_bi)$Family
V(AK4_bi)$color=gsub("Acaulosporaceae","#F8766D",V(AK4_bi)$color)
V(AK4_bi)$color=gsub("Archaeosporaceae","#CD9600",V(AK4_bi)$color)
V(AK4_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(AK4_bi)$color)
V(AK4_bi)$color=gsub("Diversisporaceae","#00BE67",V(AK4_bi)$color)
V(AK4_bi)$color=gsub("Gigasporaceae","#00BFC4",V(AK4_bi)$color)
V(AK4_bi)$color=gsub("Glomeraceae","#00A9FF",V(AK4_bi)$color)
V(AK4_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(AK4_bi)$color)
V(AK4_bi)$color=gsub("unassigned","#FF61CC",V(AK4_bi)$color)
set.seed(605)
plot(AK4_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=AK4_cust_layout,
     vertex.label=ifelse(degree(AK4_bi)>50,V(AK4_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###AK5
AK5_bi <- graph.incidence(haka_bipartite$AK5, weighted = T)
V(AK5_bi)$label<-V(AK5_bi)$name
V(AK5_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(AK5_bi)$width<-E(AK5_bi)$weight/500
deg = centr_degree(AK5_bi, mode="all")
AK5_cust_layout<-layout.reingold.tilford(AK5_bi,circular=T)
set.seed(605)
plot(AK5_bi, vertex.size=6, vertex.color=col[as.numeric(V(AK5_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=AK5_cust_layout,
     vertex.label=ifelse(degree(AK5_bi)>50,V(AK5_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#AK5 coloured by Family. 
V(AK5_bi)$Family=as.character(net_tax$Family[match(V(AK5_bi)$name,net_tax$Species)])
V(AK5_bi)$Family
V(AK5_bi)$color=V(AK5_bi)$Family
V(AK5_bi)$color=gsub("Acaulosporaceae","#F8766D",V(AK5_bi)$color)
V(AK5_bi)$color=gsub("Archaeosporaceae","#CD9600",V(AK5_bi)$color)
V(AK5_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(AK5_bi)$color)
V(AK5_bi)$color=gsub("Diversisporaceae","#00BE67",V(AK5_bi)$color)
V(AK5_bi)$color=gsub("Gigasporaceae","#00BFC4",V(AK5_bi)$color)
V(AK5_bi)$color=gsub("Glomeraceae","#00A9FF",V(AK5_bi)$color)
V(AK5_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(AK5_bi)$color)
V(AK5_bi)$color=gsub("unassigned","#FF61CC",V(AK5_bi)$color)
set.seed(605)
plot(AK5_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=AK5_cust_layout,
     vertex.label=ifelse(degree(AK5_bi)>50,V(AK5_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###AK6
AK6_bi <- graph.incidence(haka_bipartite$AK6, weighted = T)
V(AK6_bi)$label<-V(AK6_bi)$name
V(AK6_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(AK6_bi)$width<-E(AK6_bi)$weight/500
deg = centr_degree(AK6_bi, mode="all")
AK6_cust_layout<-layout.reingold.tilford(AK6_bi,circular=T)
set.seed(605)
plot(AK6_bi, vertex.size=6, vertex.color=col[as.numeric(V(AK6_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=AK6_cust_layout,
     vertex.label=ifelse(degree(AK6_bi)>50,V(AK6_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#AK6 coloured by Family. 
V(AK6_bi)$Family=as.character(net_tax$Family[match(V(AK6_bi)$name,net_tax$Species)])
V(AK6_bi)$Family
V(AK6_bi)$color=V(AK6_bi)$Family
V(AK6_bi)$color=gsub("Acaulosporaceae","#F8766D",V(AK6_bi)$color)
V(AK6_bi)$color=gsub("Archaeosporaceae","#CD9600",V(AK6_bi)$color)
V(AK6_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(AK6_bi)$color)
V(AK6_bi)$color=gsub("Diversisporaceae","#00BE67",V(AK6_bi)$color)
V(AK6_bi)$color=gsub("Gigasporaceae","#00BFC4",V(AK6_bi)$color)
V(AK6_bi)$color=gsub("Glomeraceae","#00A9FF",V(AK6_bi)$color)
V(AK6_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(AK6_bi)$color)
V(AK6_bi)$color=gsub("unassigned","#FF61CC",V(AK6_bi)$color)
set.seed(605)
plot(AK6_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=AK6_cust_layout,
     vertex.label=ifelse(degree(AK6_bi)>50,V(AK6_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

          ###Remnant ohia sites###
###RO1
ro1_bi <- graph.incidence(haka_bipartite$RO1, weighted = T)
V(ro1_bi)$label<-V(ro1_bi)$name
V(ro1_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(ro1_bi)$width<-E(ro1_bi)$weight/500
deg = centr_degree(ro1_bi, mode="all")
ro1_cust_layout<-layout.reingold.tilford(ro1_bi,circular=T)
set.seed(605)
plot(ro1_bi, vertex.size=6, vertex.color=col[as.numeric(V(ro1_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=ro1_cust_layout,
     vertex.label=ifelse(degree(ro1_bi)>50,V(ro1_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#RO1 coloured by Family. 
V(ro1_bi)$Family=as.character(net_tax$Family[match(V(ro1_bi)$name,net_tax$Species)])
V(ro1_bi)$color=V(ro1_bi)$Family
V(ro1_bi)$color=gsub("Acaulosporaceae","#F8766D",V(ro1_bi)$color)
V(ro1_bi)$color=gsub("Archaeosporaceae","#CD9600",V(ro1_bi)$color)
V(ro1_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(ro1_bi)$color)
V(ro1_bi)$color=gsub("Diversisporaceae","#00BE67",V(ro1_bi)$color)
V(ro1_bi)$color=gsub("Gigasporaceae","#00BFC4",V(ro1_bi)$color)
V(ro1_bi)$color=gsub("Glomeraceae","#00A9FF",V(ro1_bi)$color)
V(ro1_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(ro1_bi)$color)
V(ro1_bi)$color=gsub("unassigned","#FF61CC",V(ro1_bi)$color)
set.seed(605)
plot(ro1_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=ro1_cust_layout,
     vertex.label=ifelse(degree(ro1_bi)>50,V(ro1_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###RO2
RO2_bi <- graph.incidence(haka_bipartite$RO2, weighted = T)
V(RO2_bi)$label<-V(RO2_bi)$name
V(RO2_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(RO2_bi)$width<-E(RO2_bi)$weight/500
deg = centr_degree(RO2_bi, mode="all")
RO2_cust_layout<-layout.reingold.tilford(RO2_bi,circular=T)
set.seed(605)
plot(RO2_bi, vertex.size=6, vertex.color=col[as.numeric(V(RO2_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=RO2_cust_layout,
     vertex.label=ifelse(degree(RO2_bi)>50,V(RO2_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#RO2 coloured by Family. 
V(RO2_bi)$Family=as.character(net_tax$Family[match(V(RO2_bi)$name,net_tax$Species)])
V(RO2_bi)$color=V(RO2_bi)$Family
V(RO2_bi)$color=gsub("Acaulosporaceae","#F8766D",V(RO2_bi)$color)
V(RO2_bi)$color=gsub("Archaeosporaceae","#CD9600",V(RO2_bi)$color)
V(RO2_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(RO2_bi)$color)
V(RO2_bi)$color=gsub("Diversisporaceae","#00BE67",V(RO2_bi)$color)
V(RO2_bi)$color=gsub("Gigasporaceae","#00BFC4",V(RO2_bi)$color)
V(RO2_bi)$color=gsub("Glomeraceae","#00A9FF",V(RO2_bi)$color)
V(RO2_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(RO2_bi)$color)
V(RO2_bi)$color=gsub("unassigned","#FF61CC",V(RO2_bi)$color)
set.seed(605)
plot(RO2_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=RO2_cust_layout,
     vertex.label=ifelse(degree(RO2_bi)>50,V(RO2_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###RO3
RO3_bi <- graph.incidence(haka_bipartite$RO3, weighted = T)
V(RO3_bi)$label<-V(RO3_bi)$name
V(RO3_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(RO3_bi)$width<-E(RO3_bi)$weight/500
deg = centr_degree(RO3_bi, mode="all")
RO3_cust_layout<-layout.reingold.tilford(RO3_bi,circular=T)
set.seed(605)
plot(RO3_bi, vertex.size=6, vertex.color=col[as.numeric(V(RO3_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=RO3_cust_layout,
     vertex.label=ifelse(degree(RO3_bi)>50,V(RO3_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#RO3 coloured by Family. 
V(RO3_bi)$Family=as.character(net_tax$Family[match(V(RO3_bi)$name,net_tax$Species)])
V(RO3_bi)$color=V(RO3_bi)$Family
V(RO3_bi)$color=gsub("Acaulosporaceae","#F8766D",V(RO3_bi)$color)
V(RO3_bi)$color=gsub("Archaeosporaceae","#CD9600",V(RO3_bi)$color)
V(RO3_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(RO3_bi)$color)
V(RO3_bi)$color=gsub("Diversisporaceae","#00BE67",V(RO3_bi)$color)
V(RO3_bi)$color=gsub("Gigasporaceae","#00BFC4",V(RO3_bi)$color)
V(RO3_bi)$color=gsub("Glomeraceae","#00A9FF",V(RO3_bi)$color)
V(RO3_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(RO3_bi)$color)
V(RO3_bi)$color=gsub("unassigned","#FF61CC",V(RO3_bi)$color)
set.seed(605)
plot(RO3_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=RO3_cust_layout,
     vertex.label=ifelse(degree(RO3_bi)>50,V(RO3_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###RO4
RO4_bi <- graph.incidence(haka_bipartite$RO4, weighted = T)
V(RO4_bi)$label<-V(RO4_bi)$name
V(RO4_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(RO4_bi)$width<-E(RO4_bi)$weight/500
deg = centr_degree(RO4_bi, mode="all")
RO4_cust_layout<-layout.reingold.tilford(RO4_bi,circular=T)
set.seed(605)
plot(RO4_bi, vertex.size=6, vertex.color=col[as.numeric(V(RO4_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=RO4_cust_layout,
     vertex.label=ifelse(degree(RO4_bi)>50,V(RO4_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#RO4 coloured by Family. 
V(RO4_bi)$Family=as.character(net_tax$Family[match(V(RO4_bi)$name,net_tax$Species)])
V(RO4_bi)$Family
V(RO4_bi)$color=V(RO4_bi)$Family
V(RO4_bi)$color=gsub("Acaulosporaceae","#F8766D",V(RO4_bi)$color)
V(RO4_bi)$color=gsub("Archaeosporaceae","#CD9600",V(RO4_bi)$color)
V(RO4_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(RO4_bi)$color)
V(RO4_bi)$color=gsub("Diversisporaceae","#00BE67",V(RO4_bi)$color)
V(RO4_bi)$color=gsub("Gigasporaceae","#00BFC4",V(RO4_bi)$color)
V(RO4_bi)$color=gsub("Glomeraceae","#00A9FF",V(RO4_bi)$color)
V(RO4_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(RO4_bi)$color)
V(RO4_bi)$color=gsub("unassigned","#FF61CC",V(RO4_bi)$color)
set.seed(605)
plot(RO4_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=RO4_cust_layout,
     vertex.label=ifelse(degree(RO4_bi)>50,V(RO4_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###RO5
RO5_bi <- graph.incidence(haka_bipartite$RO5, weighted = T)
V(RO5_bi)$label<-V(RO5_bi)$name
V(RO5_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(RO5_bi)$width<-E(RO5_bi)$weight/500
deg = centr_degree(RO5_bi, mode="all")
RO5_cust_layout<-layout.reingold.tilford(RO5_bi,circular=T)
set.seed(605)
plot(RO5_bi, vertex.size=6, vertex.color=col[as.numeric(V(RO5_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=RO5_cust_layout,
     vertex.label=ifelse(degree(RO5_bi)>50,V(RO5_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#RO5 coloured by Family. 
V(RO5_bi)$Family=as.character(net_tax$Family[match(V(RO5_bi)$name,net_tax$Species)])
V(RO5_bi)$Family
V(RO5_bi)$color=V(RO5_bi)$Family
V(RO5_bi)$color=gsub("Acaulosporaceae","#F8766D",V(RO5_bi)$color)
V(RO5_bi)$color=gsub("Archaeosporaceae","#CD9600",V(RO5_bi)$color)
V(RO5_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(RO5_bi)$color)
V(RO5_bi)$color=gsub("Diversisporaceae","#00BE67",V(RO5_bi)$color)
V(RO5_bi)$color=gsub("Gigasporaceae","#00BFC4",V(RO5_bi)$color)
V(RO5_bi)$color=gsub("Glomeraceae","#00A9FF",V(RO5_bi)$color)
V(RO5_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(RO5_bi)$color)
V(RO5_bi)$color=gsub("unassigned","#FF61CC",V(RO5_bi)$color)
set.seed(605)
plot(RO5_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=RO5_cust_layout,
     vertex.label=ifelse(degree(RO5_bi)>50,V(RO5_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)

###RO6
RO6_bi <- graph.incidence(haka_bipartite$RO6, weighted = T)
V(RO6_bi)$label<-V(RO6_bi)$name
V(RO6_bi)$frame.color = "black"
col<-c("#F1F3CE","#265C00")
E(RO6_bi)$width<-E(RO6_bi)$weight/500
deg = centr_degree(RO6_bi, mode="all")
RO6_cust_layout<-layout.reingold.tilford(RO6_bi,circular=T)
set.seed(605)
plot(RO6_bi, vertex.size=6, vertex.color=col[as.numeric(V(RO6_bi)$type)+1],
     edge.width=1, edge.color="grey",layout=RO6_cust_layout,
     vertex.label=ifelse(degree(RO6_bi)>50,V(RO6_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)
#RO6 coloured by Family. 
V(RO6_bi)$Family=as.character(net_tax$Family[match(V(RO6_bi)$name,net_tax$Species)])
V(RO6_bi)$Family
V(RO6_bi)$color=V(RO6_bi)$Family
V(RO6_bi)$color=gsub("Acaulosporaceae","#F8766D",V(RO6_bi)$color)
V(RO6_bi)$color=gsub("Archaeosporaceae","#CD9600",V(RO6_bi)$color)
V(RO6_bi)$color=gsub("Claroideoglomeraceae","#7CAE00",V(RO6_bi)$color)
V(RO6_bi)$color=gsub("Diversisporaceae","#00BE67",V(RO6_bi)$color)
V(RO6_bi)$color=gsub("Gigasporaceae","#00BFC4",V(RO6_bi)$color)
V(RO6_bi)$color=gsub("Glomeraceae","#00A9FF",V(RO6_bi)$color)
V(RO6_bi)$color=gsub("Paraglomeraceae","#C77CFF",V(RO6_bi)$color)
V(RO6_bi)$color=gsub("unassigned","#FF61CC",V(RO6_bi)$color)
set.seed(605)
plot(RO6_bi, vertex.size=6,
     edge.width=1, edge.color="grey",layout=RO6_cust_layout,
     vertex.label=ifelse(degree(RO6_bi)>50,V(RO6_bi)$label,NA),
     vertex.label.color="black",vertex.label.font=4,vertex.label.cex=1.5,vertex.label.dist=1.5)



                                                 #############################
                                                 ### Fungal Fungal Network ###
                                                 #############################
library(devtools)
install_github("zdk123/SpiecEasi", force=TRUE)
library(SpiecEasi); library(ramify); source("taxa_summary.R", local = TRUE)
library(ggnetwork);library(intergraph);library(Matrix); library(igraph)

############################################################################################################
### Composite co-occurrence networks (habitat type x sample type) for visual and identifying keystone sp ###
############################################################################################################
### SpiecEasi ###

### Subset phyloseq object into different hab types and then by roots and soil
ro.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant Forest")
ro.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant Forest" & SampleType=="soil")
ak.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Restored Forest")
ak.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Restored Forest" & SampleType=="soil")

### Construct networks using Spiec Easi (save in case R crashes as an RDS file)
se.ro.r<- spiec.easi(ro.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro.r,"spiec_easi_files/se.ro.r.RData")
se.ro.r <- readRDS("spiec_easi_files/se.ro.r.RData")

se.ak.r<- spiec.easi(ak.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak.r,"spiec_easi_files/se.ak.r.RData")
se.ak.r <- readRDS("spiec_easi_files/se.ak.r.RData")

se.ro.s<- spiec.easi(ro.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro.s,"spiec_easi_files/se.ro.s.RData")
se.ro.s <- readRDS("spiec_easi_files/se.ro.s.RData")

se.ak.s<- spiec.easi(ak.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak.s,"spiec_easi_files/se.ak.s.RData")
se.ak.s <- readRDS("spiec_easi_files/se.ak.s.RData")


### Convert to igraph objects for plotting and computing graph statistics.
#Bring in edge weights
ig.ro.r<-adj2igraph(se.ro.r$refit$stars,vertex.attr=list(name=taxa_names(ro.r.physeq)))
ig.ak.r<-adj2igraph(se.ak.r$refit$stars,vertex.attr=list(name=taxa_names(ak.r.physeq)))
ig.ro.s<-adj2igraph(se.ro.s$refit$stars,vertex.attr=list(name=taxa_names(ro.s.physeq)))
ig.ak.s<-adj2igraph(se.ak.s$refit$stars,vertex.attr=list(name=taxa_names(ak.s.physeq)))

##RO roots plot
V(ig.ro.r)$label<-V(ig.ro.r)$name
V(ig.ro.r)$frame.color = "black"

#Calculate number of samples each AM fungal species (node) is detected, and add as vertex attribute
ro.r_dat =fast_melt(ro.r.physeq)
ro.r_no.samples = ro.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro.r_no.samples = as.data.frame(ro.r_no.samples)
V(ig.ro.r)$no.samples=as.numeric(ro.r_no.samples$no.samples[match(V(ig.ro.r)$name,ro.r_no.samples$TaxaID)])
ig.ro.r = delete_vertices(ig.ro.r,V(ig.ro.r)$no.samples == 0)

#Calculate relative abundance of each AM fungal species (node), and add as vertex attribute
ro.r.rel_abund <- summarize_taxa(ro.r.physeq,"Species")
V(ig.ro.r)$rel.abund =as.numeric(ro.r.rel_abund$meanRA[match(V(ig.ro.r)$name,ro.r_no.samples$TaxaID)])

#Important nodes
#First cacluate betweenness (nodes with high betweenness centrality represent key connector species in network)
V(ig.ro.r)$betweenness = igraph::betweenness(ig.ro.r,weights=E(ig.ro.r),normalized=FALSE)
#Then calculate node degree (nodes with high degree represent network hubs)
V(ig.ro.r)$degree = igraph::degree(ig.ro.r,mode="all",normalized=FALSE)
#Size nodes according to the two metrics
V(ig.ro.r)$size = (V(ig.ro.r)$betweenness*V(ig.ro.r)$degree)/100

#Colour by Family. First need to assign colours to match previously made graphs
brewer.pal(9,"Paired")
net_tax=read.csv("haka_feb_large_tax.csv", header = TRUE,row.names = 1)
V(ig.ro.r)$Family=as.character(net_tax$Family[match(V(ig.ro.r)$name,net_tax$Name)])
V(ig.ro.r)$color=V(ig.ro.r)$Family
V(ig.ro.r)$color=gsub("Acaulosporaceae","#8D230F",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Ambisporaceae","#4CB5F5",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Archaeosporaceae","#B7B8B6",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Claroideoglomeraceae","#3F681C",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Diversisporaceae","#FFBB00",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Geosiphonaceae","#00293C",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Gigasporaceae","#1E656D",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Glomeraceae","#F62A00",V(ig.ro.r)$color)
V(ig.ro.r)$color=gsub("Paraglomeraceae","#F1F3CE",V(ig.ro.r)$color)

par(mfrow=c(1,1))
par(mar=c(0,0,0,0))
set.seed(606)
plot(ig.ro.r, vertex.size=7, vertex.shape="circle",vertex.color=V(ig.ro.r)$color, 
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA)


##RO soil network plot
V(ig.ro.s)$label<-V(ig.ro.s)$name
V(ig.ro.s)$frame.color = "black"
ro.s_dat =fast_melt(ro.s.physeq)
ro.s_no.samples = ro.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro.s_no.samples = as.data.frame(ro.s_no.samples)
V(ig.ro.s)$no.samples=as.character(ro.s_no.samples$no.samples[match(V(ig.ro.s)$name,ro.s_no.samples$TaxaID)])
ig.ro.s = delete_vertices(ig.ro.s,V(ig.ro.s)$no.samples == 0)
ro.s.rel_abund <- summarize_taxa(ro.s.physeq,"Species")
V(ig.ro.s)$rel.abund =as.numeric(ro.s.rel_abund$meanRA[match(V(ig.ro.s)$name,ro.s_no.samples$TaxaID)])
V(ig.ro.s)$betweenness = igraph::betweenness(ig.ro.s,weights=E(ig.ro.s),normalized=FALSE)
V(ig.ro.s)$degree = igraph::degree(ig.ro.s,mode="all",normalized=FALSE)
V(ig.ro.s)$size = (V(ig.ro.s)$betweenness*V(ig.ro.s)$degree)/100
V(ig.ro.s)$Family=as.character(net_tax$Family[match(V(ig.ro.s)$name,net_tax$Name)])
V(ig.ro.s)$color=V(ig.ro.s)$Family
V(ig.ro.s)$color=gsub("Acaulosporaceae","#8D230F",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Ambisporaceae","#4CB5F5",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Archaeosporaceae","#B7B8B6",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Claroideoglomeraceae","#3F681C",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Diversisporaceae","#FFBB00",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Geosiphonaceae","#00293C",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Gigasporaceae","#1E656D",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Glomeraceae","#F62A00",V(ig.ro.s)$color)
V(ig.ro.s)$color=gsub("Paraglomeraceae","#F1F3CE",V(ig.ro.s)$color)

set.seed(606)
plot(ig.ro.s, vertex.shape="circle",vertex.color=V(ig.ro.s)$color, 
     vertex.size= 7, edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA)

##AK root network plot
V(ig.ak.r)$label<-V(ig.ak.r)$name
V(ig.ak.r)$frame.color = "black"
ak.r_dat =fast_melt(ak.r.physeq)
ak.r_no.samples = ak.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak.r_no.samples = as.data.frame(ak.r_no.samples)
V(ig.ak.r)$no.samples=as.character(ak.r_no.samples$no.samples[match(V(ig.ak.r)$name,ak.r_no.samples$TaxaID)])
ig.ak.r = delete_vertices(ig.ak.r,V(ig.ak.r)$no.samples == 0)
ak.r.rel_abund <- summarize_taxa(ak.r.physeq,"Species")
V(ig.ak.r)$rel.abund =as.numeric(ak.r.rel_abund$meanRA[match(V(ig.ak.r)$name,ak.r_no.samples$TaxaID)])
V(ig.ak.r)$betweenness = igraph::betweenness(ig.ak.r,weights=E(ig.ak.r),normalized=FALSE)
V(ig.ak.r)$degree = igraph::degree(ig.ak.r,mode="all",normalized=FALSE)
V(ig.ak.r)$size = (V(ig.ak.r)$betweenness*V(ig.ak.r)$degree)/100
V(ig.ak.r)$Family=as.character(net_tax$Family[match(V(ig.ak.r)$name,net_tax$Name)])
V(ig.ak.r)$color=V(ig.ak.r)$Family
V(ig.ak.r)$color=gsub("Acaulosporaceae","#8D230F",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Ambisporaceae","#4CB5F5",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Archaeosporaceae","#B7B8B6",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Claroideoglomeraceae","#3F681C",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Diversisporaceae","#FFBB00",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Geosiphonaceae","#00293C",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Gigasporaceae","#1E656D",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Glomeraceae","#F62A00",V(ig.ak.r)$color)
V(ig.ak.r)$color=gsub("Paraglomeraceae","#F1F3CE",V(ig.ak.r)$color)

set.seed(606)
plot(ig.ak.r, vertex.shape="circle",vertex.color=V(ig.ak.r)$color, 
     vertex.size=7, edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA)

##AK soil network plot
V(ig.ak.s)$label<-V(ig.ak.s)$name
V(ig.ak.s)$frame.color = "black"
ak.s_dat =fast_melt(ak.s.physeq)
ak.s_no.samples = ak.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak.s_no.samples = as.data.frame(ak.s_no.samples)
V(ig.ak.s)$no.samples=as.character(ak.s_no.samples$no.samples[match(V(ig.ak.s)$name,ak.s_no.samples$TaxaID)])
ig.ak.s = delete_vertices(ig.ak.s,V(ig.ak.s)$no.samples == 0)
ak.s.rel_abund <- summarize_taxa(ak.s.physeq,"Species")
V(ig.ak.s)$rel.abund =as.numeric(ak.s.rel_abund$meanRA[match(V(ig.ak.s)$name,ak.s_no.samples$TaxaID)])
V(ig.ak.s)$betweenness = igraph::betweenness(ig.ak.s,weights=E(ig.ak.s),normalized=FALSE)
V(ig.ak.s)$degree = igraph::degree(ig.ak.s,mode="all",normalized=FALSE)
V(ig.ak.s)$size = (V(ig.ak.s)$betweenness*V(ig.ak.s)$degree)/100
V(ig.ak.s)$Family=as.character(net_tax$Family[match(V(ig.ak.s)$name,net_tax$Name)])
V(ig.ak.s)$color=V(ig.ak.s)$Family
V(ig.ak.s)$color=gsub("Acaulosporaceae","#8D230F",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Ambisporaceae","#4CB5F5",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Archaeosporaceae","#B7B8B6",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Claroideoglomeraceae","#3F681C",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Diversisporaceae","#FFBB00",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Geosiphonaceae","#00293C",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Gigasporaceae","#1E656D",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Glomeraceae","#F62A00",V(ig.ak.s)$color)
V(ig.ak.s)$color=gsub("Paraglomeraceae","#F1F3CE",V(ig.ak.s)$color)

set.seed(606)
plot(ig.ak.s, vertex.shape="circle",vertex.color=V(ig.ak.s)$color, 
     vertex.size = 7, edge.width=3, edge.color="dark grey",layout=ir.ak.s_cust_layout,
     vertex.label=NA)

#Manuscript plot
par(mar=c(0.75,0.75,0.75,0.75))
par(mfrow=c(2,1))
set.seed(606)
plot(ig.ro.r, vertex.size=7, vertex.shape="circle",vertex.color=V(ig.ro.r)$color, 
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Remnant Forest")
set.seed(606)
plot(ig.ak.r, vertex.shape="circle",vertex.color=V(ig.ak.r)$color, 
     vertex.size= 7, edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Restored Forest")

legend(title="AM fungal Family",'right',
       c("Acaulosporaceae","Ambisporaceae","Archaeosporaceae","Claroideoglomeraceae",
         "Diversisporaceae","Geosiphonaceae","Gigasporaceae","Glomeraceae","Paraglomeraceae"),pch=21,
       pt.bg=c("#8D230F","#4CB5F5","#B7B8B6","#3F681C","#FFBB00","#00293C","#1E656D","#F62A00","#F1F3CE"),pt.cex=4,
       box.lty=0,cex=1,text.font=2)

##############################
#####  CHARACTERISTICS   #####
##############################
#Betweenness centrality
ig.ro.r_between <- igraph::betweenness(ig.ro.r,weights=E(ig.ro.r),normalized=TRUE)
ig.ro.s_between <- igraph::betweenness(ig.ro.s,weights=E(ig.ro.s),normalized=TRUE)
ig.ak.r_between <- igraph::betweenness(ig.ak.r,weights=E(ig.ak.r),normalized=TRUE)
ig.ak.s_between <- igraph::betweenness(ig.ak.s,weights=E(ig.ak.s),normalized=TRUE)

#Connectedness (degree)
ig.ro1.r_degree <- igraph::degree(ig.ro1.r, mode="all", normalized=TRUE)
ig.ro2.r_degree <- igraph::degree(ig.ro2.r, mode="all", normalized=TRUE)
ig.ro3.r_degree <- igraph::degree(ig.ro3.r, mode="all", normalized=TRUE)
ig.ro4.r_degree <- igraph::degree(ig.ro4.r, mode="all", normalized=TRUE)


## Welch t-tests
#Centrality
centrality_between_hab_roots <- t.test(ig.ro.r_between,ig.ak.r_between,
                                          paired=FALSE, var.equal=FALSE)
centrality_between_hab_roots

centrality_between_hab_soil <- t.test(ig.ro.s_between,ig.ak.s_between,
                                         paired=FALSE, var.equal=FALSE)
centrality_between_hab_soil

centrality_within_RO <- t.test(ig.ro.r_between, ig.ro.s_between,
                                  paired=FALSE, var.equal=FALSE)
centrality_within_RO

connectedness_within_AK <- t.test(subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "roots")$Connectedness,
                                  subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "soil")$Connectedness,
                                  paired=FALSE, var.equal=FALSE)
connectedness_within_AK

#Density
Density_between_hab_roots <- t.test(subset(fungal_networks,
                                           HabitatType == "Afforested koa" & SampleType == "roots")$Density,
                                    subset(fungal_networks,
                                           HabitatType == "Remnant ohia" & SampleType == "roots")$Density,
                                    paired=FALSE, var.equal=FALSE)
Density_between_hab_roots

Density_between_hab_soil <- t.test(subset(fungal_networks,
                                          HabitatType == "Afforested koa" & SampleType == "soil")$Density,
                                   subset(fungal_networks,
                                          HabitatType == "Remnant ohia" & SampleType == "soil")$Density,
                                   paired=FALSE, var.equal=FALSE)
Density_between_hab_soil

Density_within_RO <- t.test(subset(fungal_networks,
                                   HabitatType == "Remnant ohia" & SampleType == "roots")$Density,
                            subset(fungal_networks,
                                   HabitatType == "Remnant ohia" & SampleType == "soil")$Density,
                            paired=FALSE, var.equal=FALSE)
Density_within_RO

Density_within_AK <- t.test(subset(fungal_networks,
                                   HabitatType == "Afforested koa" & SampleType == "roots")$Density,
                            subset(fungal_networks,
                                   HabitatType == "Afforested koa" & SampleType == "soil")$Density,
                            paired=FALSE, var.equal=FALSE)
Density_within_AK


##############################
#####   KEYSTONE SPECIES   ###
##############################
#Remnant ohia roots
ro.r.between <- as.data.frame(V(ig.ro.r)$betweenness,normalized=TRUE)
ro.r.degree <- as.data.frame(V(ig.ro.r)$degree,mode="all",normalized=TRUE)
ro.r.rel.abund <- as.data.frame(V(ig.ro.r)$rel.abund)
ro.r.fam <- as.data.frame(V(ig.ro.r)$Family)
ro.r.species <- as.data.frame(net_tax$Species[match(V(ig.ro.r)$name,net_tax$Name)])
ro.r.no.samples <- as.data.frame(V(ig.ro.r)$no.samples)
sample_data(ro.r.physeq)
ro.r.total.samples <- as.data.frame(rep(237),times=162)
ro.r.keystone <- cbind(ro.r.between,ro.r.degree,ro.r.rel.abund,ro.r.fam,ro.r.species,ro.r.no.samples,ro.r.total.samples)
rownames(ro.r.keystone) <- V(ig.ro.r)$name
colnames(ro.r.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family","Species","No.Samples","TotalSamples")
ro.r.keystone$No.Samples <- as.numeric(ro.r.keystone$No.Samples)
ro.r.keystone$perc.samples <- ro.r.keystone$No.Samples / ro.r.keystone$TotalSamples
ro.r.keystone$Prevalence <- ro.r.keystone$perc.samples * ro.r.keystone$RelativeAbundance


ro.r.keystone.plot <- ggplot(ro.r.keystone,aes(x=Degree,y=Betweenness)) +
  geom_point(aes(size=Prevalence,colour=Family),position="jitter") +
  ylab("Betweenness Centrality (normalized)") + 
  xlab("Node degree") +
  scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x)10^x),
                     labels=trans_format("log10",math_format(10^.x)),limits=c(10^0,10^4)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  scale_colour_manual(values=c("#34888C","#8D230F","#B7B8B6","#3F681C",
                               "#FFBB00","#00293C","#1E656D","#F62A00","#F1F3CE")) +
  guides(colour=guide_legend(override.aes=list(size=5))) +
  ggtitle("A)")

ro.r.keystone.plot

#Afforested koa roots
ak.r.between <- as.data.frame(V(ig.ak.r)$betweenness,normalized=TRUE)
ak.r.degree <- as.data.frame(V(ig.ak.r)$degree,mode="all",normalized=TRUE)
ak.r.rel.abund <- as.data.frame(V(ig.ak.r)$rel.abund)
ak.r.fam <- as.data.frame(V(ig.ak.r)$Family)
ak.r.species <- as.data.frame(net_tax$Species[match(V(ig.ak.r)$name,net_tax$Name)])
ak.r.no.samples <- as.data.frame(V(ig.ak.r)$no.samples)
sample_data(ak.r.physeq)
ak.r.total.samples <- as.data.frame(rep(281),times=167)
ak.r.keystone <- cbind(ak.r.between,ak.r.degree,ak.r.rel.abund,ak.r.fam,ak.r.species,ak.r.no.samples,ak.r.total.samples)
rownames(ak.r.keystone) <- V(ig.ak.r)$name
colnames(ak.r.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family","Species","No.Samples","TotalSamples")
ak.r.keystone$No.Samples <- as.numeric(ak.r.keystone$No.Samples)
ak.r.keystone$perc.samples <- ak.r.keystone$No.Samples / ak.r.keystone$TotalSamples
ak.r.keystone$Prevalence <- ak.r.keystone$perc.samples * ak.r.keystone$RelativeAbundance

ak.r.keystone.plot <- ggplot(ak.r.keystone,aes(x=Degree,y=Betweenness)) +
  geom_point(aes(size=Prevalence,colour=Family),position="jitter") +
  ylab("Betweenness Centrality (normalized)") + xlab("Node degree") +
  scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x)10^x),
                     labels=trans_format("log10",math_format(10^.x)),limits=c(10^0,10^4)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  scale_colour_manual(values=c("#34888C","#8D230F","#B7B8B6","#3F681C",
                               "#FFBB00","#00293C","#1E656D","#F62A00","#F1F3CE")) +
  guides(colour=guide_legend(override.aes=list(size=5))) +
  ggtitle("B)")

ak.r.keystone.plot

#Remnant ohia soil
ro.s.between <- as.data.frame(V(ig.ro.s)$betweenness)
ro.s.degree <- as.data.frame(V(ig.ro.s)$degree)
ro.s.rel.abund <- as.data.frame(V(ig.ro.s)$rel.abund)
ro.s.fam <- as.data.frame(V(ig.ro.s)$Family)
ro.s.species <- as.data.frame(net_tax$Species[match(V(ig.ro.s)$name,net_tax$Name)])
ro.s.no.samples <- as.data.frame(V(ig.ro.s)$no.samples)
sample_data(ro.s.physeq)
ro.s.total.samples <- as.data.frame(rep(212),times=138)
ro.s.keystone <- cbind(ro.s.between,ro.s.degree,ro.s.rel.abund,ro.s.fam,ro.s.species,ro.s.no.samples,ro.s.total.samples)
rownames(ro.s.keystone) <- V(ig.ro.s)$name
colnames(ro.s.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family","Species","No.Samples","TotalSamples")
ro.s.keystone$No.Samples <- as.numeric(ro.s.keystone$No.Samples)
ro.s.keystone$perc.samples <- ro.s.keystone$No.Samples / ro.s.keystone$TotalSamples
ro.s.keystone$Prevalence <- ro.s.keystone$perc.samples * ro.s.keystone$RelativeAbundance

ro.s.keystone.plot <- ggplot(ro.s.keystone,aes(x=Degree,y=Betweenness)) +
  geom_point(aes(size=Prevalence,colour=Family),position="jitter") +
  ylab("Betweenness Centrality") + xlab("Node degree") +
  scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x)10^x),
                     labels=trans_format("log10",math_format(10^.x)),limits=c(imits=c(10^0,10^4))) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  scale_colour_manual(values=c("#34888C","#8D230F","#B7B8B6","#3F681C",
                               "#FFBB00","#00293C","#1E656D","#F62A00","#F1F3CE")) +
  guides(colour=guide_legend(override.aes=list(size=5))) +
  ggtitle("B)")

ro.s.keystone.plot


#Afforested koa soil
ak.s.between <- as.data.frame(V(ig.ak.s)$betweenness)
ak.s.degree <- as.data.frame(V(ig.ak.s)$degree)
ak.s.rel.abund <- as.data.frame(V(ig.ak.s)$rel.abund)
ak.s.fam <- as.data.frame(V(ig.ak.s)$Family)
ak.s.species <- as.data.frame(net_tax$Species[match(V(ig.ak.s)$name,net_tax$Name)])
ak.s.no.samples <- as.data.frame(V(ig.ak.s)$no.samples)
sample_data(ak.s.physeq)
ak.s.total.samples <- as.data.frame(rep(212),times=138)
ak.s.keystone <- cbind(ak.s.between,ak.s.degree,ak.s.rel.abund,ak.s.fam,ak.s.species,ak.s.no.samples,ak.s.total.samples)
rownames(ak.s.keystone) <- V(ig.ak.s)$name
colnames(ak.s.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family","Species","No.Samples","TotalSamples")
ak.s.keystone$No.Samples <- as.numeric(ak.s.keystone$No.Samples)
ak.s.keystone$perc.samples <- ak.s.keystone$No.Samples / ak.s.keystone$TotalSamples
ak.s.keystone$Prevalence <- ak.s.keystone$perc.samples * ak.s.keystone$RelativeAbundance

ak.s.keystone.plot <- ggplot(ak.s.keystone,aes(x=Degree,y=Betweenness)) +
  geom_point(aes(size=Prevalence,colour=Family),position="jitter") +
  ylab("Betweenness Centrality") + xlab("Node degree") +
  scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x)10^x),
                     labels=trans_format("log10",math_format(10^.x)),limits=c(10^0,10^4)) +
  theme(text=element_text(colour="black",size=15)) + 
  theme(axis.text.x=element_text(hjust=1,colour="black",size=15)) +
  theme(axis.text.y=element_text(colour="black",size=15)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank()) +
  scale_colour_manual(values=c("#34888C","#8D230F","#B7B8B6","#3F681C",
                               "#FFBB00","#00293C","#1E656D","#F62A00","#F1F3CE")) +
  guides(colour=guide_legend(override.aes=list(size=5))) +
  ggtitle("D)")

ak.s.keystone.plot

#Keystone species plot
keystone_plot<- grid.arrange(ro.r.keystone.plot + theme(legend.position = "none"), 
                             ak.r.keystone.plot + theme(legend.position = "none"), 
                             ncol=2,nrow=1)
ggsave("Fig_4_keystone_species.tiff", 
       plot = keystone_plot, width = 12, height = 5)

#Export keystone tables
write.table(ro.r.keystone,"keystone_tables/ro.r.keystone.txt",sep="\t")
write.table(ro.s.keystone,"keystone_tables/ro.s.keystone.txt",sep="\t")
write.table(ak.r.keystone,"keystone_tables/ak.r.keystone.txt",sep="\t")
write.table(ak.s.keystone,"keystone_tables/ak.s.keystone.txt",sep="\t")

##############################################################################################################
### Co-occurrence networks by plot (habitat type x sample type x plot) for overall network characteristics ###
##############################################################################################################
#### Subset phyloseq objects ####
ro1.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="roots" & Plot =="RO1")
ro2.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="roots" & Plot =="RO2")
ro3.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="roots" & Plot =="RO3")
ro4.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="roots" & Plot =="RO4")
ro5.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="roots" & Plot =="RO5")
ro6.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="roots" & Plot =="RO6")
ro1.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="soil" & Plot =="RO1")
ro2.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="soil" & Plot =="RO2")
ro3.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="soil" & Plot =="RO3")
ro4.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="soil" & Plot =="RO4")
ro5.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="soil" & Plot =="RO5")
ro6.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Remnant ohia" & SampleType=="soil" & Plot =="RO6")
ak1.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="roots" & Plot =="AK1")
ak2.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="roots" & Plot =="AK2")
ak3.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="roots" & Plot =="AK3")
ak4.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="roots" & Plot =="AK4")
ak5.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="roots" & Plot =="AK5")
ak6.r.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="roots" & Plot =="AK6")
ak1.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="soil" & Plot == "AK1")
ak2.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="soil" & Plot == "AK2")
ak3.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="soil" & Plot == "AK3")
ak4.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="soil" & Plot == "AK4")
ak5.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="soil" & Plot == "AK5")
ak6.s.physeq = subset_samples(vt_rel_abund,HabitatType == "Afforested koa" & SampleType=="soil" & Plot == "AK6")

#### Make Spiec-Easi objects and SAVE ####
se.ro1.r<- spiec.easi(ro1.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro1.r,"spiec_easi_files/se.ro1.r.RData")
se.ro2.r<- spiec.easi(ro2.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro2.r,"spiec_easi_files/se.ro2.r.RData")
se.ro3.r<- spiec.easi(ro3.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro3.r,"spiec_easi_files/se.ro3.r.RData")
se.ro4.r<- spiec.easi(ro4.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro4.r,"spiec_easi_files/se.ro4.r.RData")
se.ro5.r<- spiec.easi(ro5.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro5.r,"spiec_easi_files/se.ro5.r.RData")
se.ro6.r<- spiec.easi(ro6.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro6.r,"spiec_easi_files/se.ro6.r.RData")
se.ro1.s<- spiec.easi(ro1.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro1.s,"spiec_easi_files/se.ro1.s.RData")
se.ro2.s<- spiec.easi(ro2.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro2.s,"spiec_easi_files/se.ro2.s.RData")
se.ro3.s<- spiec.easi(ro3.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro3.s,"spiec_easi_files/se.ro3.s.RData")
se.ro4.s<- spiec.easi(ro4.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro4.s,"spiec_easi_files/se.ro4.s.RData")
se.ro5.s<- spiec.easi(ro5.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro5.s,"spiec_easi_files/se.ro5.s.RData")
se.ro6.s<- spiec.easi(ro6.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ro6.s,"spiec_easi_files/se.ro6.s.RData")
se.ak1.r<- spiec.easi(ak1.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak1.r,"spiec_easi_files/se.ak1.r.RData")
se.ak2.r<- spiec.easi(ak2.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak2.r,"spiec_easi_files/se.ak2.r.RData")
se.ak3.r<- spiec.easi(ak3.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak3.r,"spiec_easi_files/se.ak3.r.RData")
se.ak4.r<- spiec.easi(ak4.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak4.r,"spiec_easi_files/se.ak4.r.RData")
se.ak5.r<- spiec.easi(ak5.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak5.r,"spiec_easi_files/se.ak5.r.RData")
se.ak6.r<- spiec.easi(ak6.r.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak6.r,"spiec_easi_files/se.ak6.r.RData")
se.ak1.s<- spiec.easi(ak1.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak1.s,"spiec_easi_files/se.ak1.s.RData")
se.ak2.s<- spiec.easi(ak2.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak2.s,"spiec_easi_files/se.ak2.s.RData")
se.ak3.s<- spiec.easi(ak3.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak3.s,"spiec_easi_files/se.ak3.s.RData")
se.ak4.s<- spiec.easi(ak4.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak4.s,"spiec_easi_files/se.ak4.s.RData")
se.ak5.s<- spiec.easi(ak5.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak5.s,"spiec_easi_files/se.ak5.s.RData")
se.ak6.s<- spiec.easi(ak6.s.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,icov.select.params=list(rep.num=100))
saveRDS(se.ak6.s,"spiec_easi_files/se.ak6.s.RData")

#### Bring networks back into R (in case R crashes like a mother fucker)    ####
se.ro1.r <- readRDS("spiec_easi_files/se.ro1.r.RData")
se.ro2.r <- readRDS("spiec_easi_files/se.ro2.r.RData")
se.ro3.r <- readRDS("spiec_easi_files/se.ro3.r.RData")
se.ro4.r <- readRDS("spiec_easi_files/se.ro4.r.RData")
se.ro5.r <- readRDS("spiec_easi_files/se.ro5.r.RData")
se.ro6.r <- readRDS("spiec_easi_files/se.ro6.r.RData")
se.ro1.s <- readRDS("spiec_easi_files/se.ro1.s.RData")
se.ro2.s <- readRDS("spiec_easi_files/se.ro2.s.RData")
se.ro3.s <- readRDS("spiec_easi_files/se.ro3.s.RData")
se.ro4.s <- readRDS("spiec_easi_files/se.ro4.s.RData")
se.ro5.s <- readRDS("spiec_easi_files/se.ro5.s.RData")
se.ro6.s <- readRDS("spiec_easi_files/se.ro6.s.RData")
se.ak1.r <- readRDS("spiec_easi_files/se.ak1.r.RData")
se.ak2.r <- readRDS("spiec_easi_files/se.ak2.r.RData")
se.ak3.r <- readRDS("spiec_easi_files/se.ak3.r.RData")
se.ak4.r <- readRDS("spiec_easi_files/se.ak4.r.RData")
se.ak5.r <- readRDS("spiec_easi_files/se.ak5.r.RData")
se.ak6.r <- readRDS("spiec_easi_files/se.ak6.r.RData")
se.ak1.s <- readRDS("spiec_easi_files/se.ak1.s.RData")
se.ak2.s <- readRDS("spiec_easi_files/se.ak2.s.RData")
se.ak3.s <- readRDS("spiec_easi_files/se.ak3.s.RData")
se.ak4.s <- readRDS("spiec_easi_files/se.ak4.s.RData")
se.ak5.s <- readRDS("spiec_easi_files/se.ak5.s.RData")
se.ak6.s <- readRDS("spiec_easi_files/se.ak6.s.RData")

#### Convert to igraph networks ####
ig.ro1.r<-adj2igraph(se.ro1.r$refit$stars,vertex.attr=list(name=taxa_names(ro1.r.physeq)))
ro1.r_dat =fast_melt(ro1.r.physeq)
ro1.r_no.samples = ro1.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro1.r_no.samples = as.data.frame(ro1.r_no.samples)
V(ig.ro1.r)$no.samples=as.numeric(ro1.r_no.samples$no.samples[match(V(ig.ro1.r)$name,ro1.r_no.samples$TaxaID)])
ig.ro1.r = delete_vertices(ig.ro1.r,V(ig.ro1.r)$no.samples == 0)

ig.ro2.r<-adj2igraph(se.ro2.r$refit$stars,vertex.attr=list(name=taxa_names(ro2.r.physeq)))
ro2.r_dat =fast_melt(ro2.r.physeq)
ro2.r_no.samples = ro2.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro2.r_no.samples = as.data.frame(ro2.r_no.samples)
V(ig.ro2.r)$no.samples=as.numeric(ro2.r_no.samples$no.samples[match(V(ig.ro2.r)$name,ro2.r_no.samples$TaxaID)])
ig.ro2.r = delete_vertices(ig.ro2.r,V(ig.ro2.r)$no.samples == 0)

ig.ro3.r<-adj2igraph(se.ro3.r$refit$stars,vertex.attr=list(name=taxa_names(ro3.r.physeq)))
ro3.r_dat =fast_melt(ro3.r.physeq)
ro3.r_no.samples = ro3.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro3.r_no.samples = as.data.frame(ro3.r_no.samples)
V(ig.ro3.r)$no.samples=as.numeric(ro3.r_no.samples$no.samples[match(V(ig.ro3.r)$name,ro3.r_no.samples$TaxaID)])
ig.ro3.r = delete_vertices(ig.ro3.r,V(ig.ro3.r)$no.samples == 0)

ig.ro4.r<-adj2igraph(se.ro4.r$refit$stars,vertex.attr=list(name=taxa_names(ro4.r.physeq)))
ro4.r_dat =fast_melt(ro4.r.physeq)
ro4.r_no.samples = ro4.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro4.r_no.samples = as.data.frame(ro4.r_no.samples)
V(ig.ro4.r)$no.samples=as.numeric(ro4.r_no.samples$no.samples[match(V(ig.ro4.r)$name,ro4.r_no.samples$TaxaID)])
ig.ro4.r = delete_vertices(ig.ro4.r,V(ig.ro4.r)$no.samples == 0)

ig.ro5.r<-adj2igraph(se.ro5.r$refit$stars,vertex.attr=list(name=taxa_names(ro5.r.physeq)))
ro5.r_dat =fast_melt(ro5.r.physeq)
ro5.r_no.samples = ro5.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro5.r_no.samples = as.data.frame(ro5.r_no.samples)
V(ig.ro5.r)$no.samples=as.numeric(ro5.r_no.samples$no.samples[match(V(ig.ro5.r)$name,ro5.r_no.samples$TaxaID)])
ig.ro5.r = delete_vertices(ig.ro5.r,V(ig.ro5.r)$no.samples == 0)

ig.ro6.r<-adj2igraph(se.ro6.r$refit$stars,vertex.attr=list(name=taxa_names(ro6.r.physeq)))
ro6.r_dat =fast_melt(ro6.r.physeq)
ro6.r_no.samples = ro6.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro6.r_no.samples = as.data.frame(ro6.r_no.samples)
V(ig.ro6.r)$no.samples=as.numeric(ro6.r_no.samples$no.samples[match(V(ig.ro6.r)$name,ro6.r_no.samples$TaxaID)])
ig.ro6.r = delete_vertices(ig.ro6.r,V(ig.ro6.r)$no.samples == 0)

ig.ro1.s<-adj2igraph(se.ro1.s$refit$stars,vertex.attr=list(name=taxa_names(ro1.s.physeq)))
ro1.s_dat =fast_melt(ro1.s.physeq)
ro1.s_no.samples = ro1.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro1.s_no.samples = as.data.frame(ro1.s_no.samples)
V(ig.ro1.s)$no.samples=as.numeric(ro1.s_no.samples$no.samples[match(V(ig.ro1.s)$name,ro1.s_no.samples$TaxaID)])
ig.ro1.s = delete_vertices(ig.ro1.s,V(ig.ro1.s)$no.samples == 0)

ig.ro2.s<-adj2igraph(se.ro2.s$refit$stars,vertex.attr=list(name=taxa_names(ro2.s.physeq)))
ro2.s_dat =fast_melt(ro2.s.physeq)
ro2.s_no.samples = ro2.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro2.s_no.samples = as.data.frame(ro2.s_no.samples)
V(ig.ro2.s)$no.samples=as.numeric(ro2.s_no.samples$no.samples[match(V(ig.ro2.s)$name,ro2.s_no.samples$TaxaID)])
ig.ro2.s = delete_vertices(ig.ro2.s,V(ig.ro2.s)$no.samples == 0)

ig.ro3.s<-adj2igraph(se.ro3.s$refit$stars,vertex.attr=list(name=taxa_names(ro3.s.physeq)))
ro3.s_dat =fast_melt(ro3.s.physeq)
ro3.s_no.samples = ro3.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro3.s_no.samples = as.data.frame(ro3.s_no.samples)
V(ig.ro3.s)$no.samples=as.numeric(ro3.s_no.samples$no.samples[match(V(ig.ro3.s)$name,ro3.s_no.samples$TaxaID)])
ig.ro3.s = delete_vertices(ig.ro3.s,V(ig.ro3.s)$no.samples == 0)

ig.ro4.s<-adj2igraph(se.ro4.s$refit$stars,vertex.attr=list(name=taxa_names(ro4.s.physeq)))
ro4.s_dat =fast_melt(ro4.s.physeq)
ro4.s_no.samples = ro4.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro4.s_no.samples = as.data.frame(ro4.s_no.samples)
V(ig.ro4.s)$no.samples=as.numeric(ro4.s_no.samples$no.samples[match(V(ig.ro4.s)$name,ro4.s_no.samples$TaxaID)])
ig.ro4.s = delete_vertices(ig.ro4.s,V(ig.ro4.s)$no.samples == 0)

ig.ro5.s<-adj2igraph(se.ro5.s$refit$stars,vertex.attr=list(name=taxa_names(ro5.s.physeq)))
ro5.s_dat =fast_melt(ro5.s.physeq)
ro5.s_no.samples = ro5.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro5.s_no.samples = as.data.frame(ro5.s_no.samples)
V(ig.ro5.s)$no.samples=as.numeric(ro5.s_no.samples$no.samples[match(V(ig.ro5.s)$name,ro5.s_no.samples$TaxaID)])
ig.ro5.s = delete_vertices(ig.ro5.s,V(ig.ro5.s)$no.samples == 0)

ig.ro6.s<-adj2igraph(se.ro6.s$refit$stars,vertex.attr=list(name=taxa_names(ro6.s.physeq)))
ro6.s_dat =fast_melt(ro6.s.physeq)
ro6.s_no.samples = ro6.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ro6.s_no.samples = as.data.frame(ro6.s_no.samples)
V(ig.ro6.s)$no.samples=as.numeric(ro6.s_no.samples$no.samples[match(V(ig.ro6.s)$name,ro6.s_no.samples$TaxaID)])
ig.ro6.s = delete_vertices(ig.ro6.s,V(ig.ro6.s)$no.samples == 0)

ig.ak1.r<-adj2igraph(se.ak1.r$refit$stars,vertex.attr=list(name=taxa_names(ak1.r.physeq)))
ak1.r_dat =fast_melt(ak1.r.physeq)
ak1.r_no.samples = ak1.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak1.r_no.samples = as.data.frame(ak1.r_no.samples)
V(ig.ak1.r)$no.samples=as.numeric(ak1.r_no.samples$no.samples[match(V(ig.ak1.r)$name,ak1.r_no.samples$TaxaID)])
ig.ak1.r = delete_vertices(ig.ak1.r,V(ig.ak1.r)$no.samples == 0)

ig.ak2.r<-adj2igraph(se.ak2.r$refit$stars,vertex.attr=list(name=taxa_names(ak2.r.physeq)))
ak2.r_dat =fast_melt(ak2.r.physeq)
ak2.r_no.samples = ak2.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak2.r_no.samples = as.data.frame(ak2.r_no.samples)
V(ig.ak2.r)$no.samples=as.numeric(ak2.r_no.samples$no.samples[match(V(ig.ak2.r)$name,ak2.r_no.samples$TaxaID)])
ig.ak2.r = delete_vertices(ig.ak2.r,V(ig.ak2.r)$no.samples == 0)

ig.ak3.r<-adj2igraph(se.ak3.r$refit$stars,vertex.attr=list(name=taxa_names(ak3.r.physeq)))
ak3.r_dat =fast_melt(ak3.r.physeq)
ak3.r_no.samples = ak3.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak3.r_no.samples = as.data.frame(ak3.r_no.samples)
V(ig.ak3.r)$no.samples=as.numeric(ak3.r_no.samples$no.samples[match(V(ig.ak3.r)$name,ak3.r_no.samples$TaxaID)])
ig.ak3.r = delete_vertices(ig.ak3.r,V(ig.ak3.r)$no.samples == 0)

ig.ak4.r<-adj2igraph(se.ak4.r$refit$stars,vertex.attr=list(name=taxa_names(ak4.r.physeq)))
ak4.r_dat =fast_melt(ak4.r.physeq)
ak4.r_no.samples = ak4.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak4.r_no.samples = as.data.frame(ak4.r_no.samples)
V(ig.ak4.r)$no.samples=as.numeric(ak4.r_no.samples$no.samples[match(V(ig.ak4.r)$name,ak4.r_no.samples$TaxaID)])
ig.ak4.r = delete_vertices(ig.ak4.r,V(ig.ak4.r)$no.samples == 0)

ig.ak5.r<-adj2igraph(se.ak5.r$refit$stars,vertex.attr=list(name=taxa_names(ak5.r.physeq)))
ak5.r_dat =fast_melt(ak5.r.physeq)
ak5.r_no.samples = ak5.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak5.r_no.samples = as.data.frame(ak5.r_no.samples)
V(ig.ak5.r)$no.samples=as.numeric(ak5.r_no.samples$no.samples[match(V(ig.ak5.r)$name,ak5.r_no.samples$TaxaID)])
ig.ak5.r = delete_vertices(ig.ak5.r,V(ig.ak5.r)$no.samples == 0)

ig.ak6.r<-adj2igraph(se.ak6.r$refit$stars,vertex.attr=list(name=taxa_names(ak6.r.physeq)))
ak6.r_dat =fast_melt(ak6.r.physeq)
ak6.r_no.samples = ak6.r_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak6.r_no.samples = as.data.frame(ak6.r_no.samples)
V(ig.ak6.r)$no.samples=as.numeric(ak6.r_no.samples$no.samples[match(V(ig.ak6.r)$name,ak6.r_no.samples$TaxaID)])
ig.ak6.r = delete_vertices(ig.ak6.r,V(ig.ak6.r)$no.samples == 0)

ig.ak1.s<-adj2igraph(se.ak1.s$refit$stars,vertex.attr=list(name=taxa_names(ak1.s.physeq)))
ak1.s_dat =fast_melt(ak1.s.physeq)
ak1.s_no.samples = ak1.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak1.s_no.samples = as.data.frame(ak1.s_no.samples)
V(ig.ak1.s)$no.samples=as.numeric(ak1.s_no.samples$no.samples[match(V(ig.ak1.s)$name,ak1.s_no.samples$TaxaID)])
ig.ak1.s = delete_vertices(ig.ak1.s,V(ig.ak1.s)$no.samples == 0)

ig.ak2.s<-adj2igraph(se.ak2.s$refit$stars,vertex.attr=list(name=taxa_names(ak2.s.physeq)))
ak2.s_dat =fast_melt(ak2.s.physeq)
ak2.s_no.samples = ak2.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak2.s_no.samples = as.data.frame(ak2.s_no.samples)
V(ig.ak2.s)$no.samples=as.numeric(ak2.s_no.samples$no.samples[match(V(ig.ak2.s)$name,ak2.s_no.samples$TaxaID)])
ig.ak2.s = delete_vertices(ig.ak2.s,V(ig.ak2.s)$no.samples == 0)

ig.ak3.s<-adj2igraph(se.ak3.s$refit$stars,vertex.attr=list(name=taxa_names(ak3.s.physeq)))
ak3.s_dat =fast_melt(ak3.s.physeq)
ak3.s_no.samples = ak3.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak3.s_no.samples = as.data.frame(ak3.s_no.samples)
V(ig.ak3.s)$no.samples=as.numeric(ak3.s_no.samples$no.samples[match(V(ig.ak3.s)$name,ak3.s_no.samples$TaxaID)])
ig.ak3.s = delete_vertices(ig.ak3.s,V(ig.ak3.s)$no.samples == 0)

ig.ak4.s<-adj2igraph(se.ak4.s$refit$stars,vertex.attr=list(name=taxa_names(ak4.s.physeq)))
ak4.s_dat =fast_melt(ak4.s.physeq)
ak4.s_no.samples = ak4.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak4.s_no.samples = as.data.frame(ak4.s_no.samples)
V(ig.ak4.s)$no.samples=as.numeric(ak4.s_no.samples$no.samples[match(V(ig.ak4.s)$name,ak4.s_no.samples$TaxaID)])
ig.ak4.s = delete_vertices(ig.ak4.s,V(ig.ak4.s)$no.samples == 0)

ig.ak5.s<-adj2igraph(se.ak5.s$refit$stars,vertex.attr=list(name=taxa_names(ak5.s.physeq)))
ak5.s_dat =fast_melt(ak5.s.physeq)
ak5.s_no.samples = ak5.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak5.s_no.samples = as.data.frame(ak5.s_no.samples)
V(ig.ak5.s)$no.samples=as.numeric(ak5.s_no.samples$no.samples[match(V(ig.ak5.s)$name,ak5.s_no.samples$TaxaID)])
ig.ak5.s = delete_vertices(ig.ak5.s,V(ig.ak5.s)$no.samples == 0)

ig.ak6.s<-adj2igraph(se.ak6.s$refit$stars,vertex.attr=list(name=taxa_names(ak6.s.physeq)))
ak6.s_dat =fast_melt(ak6.s.physeq)
ak6.s_no.samples = ak6.s_dat[, list(no.samples=sum(count>0)),by=TaxaID]
ak6.s_no.samples = as.data.frame(ak6.s_no.samples)
V(ig.ak6.s)$no.samples=as.numeric(ak6.s_no.samples$no.samples[match(V(ig.ak6.s)$name,ak6.s_no.samples$TaxaID)])
ig.ak6.s = delete_vertices(ig.ak6.s,V(ig.ak6.s)$no.samples == 0)

obs_fungal_networks <- as.list(ig.ro1.r,ig.ro2.r,ig.ro3.r,ig.ro4.r,
                               ig.ro5.r,ig.ro6.r,ig.ro1.s,ig.ro2.s,ig.ro3.s,ig.ro4.s,ig.ro5.s,ig.ro6.s,ig.ak1.r,ig.ak2.r,
                               ig.ak3.r,ig.ak4.r,ig.ak5.r,ig.ak6.r,ig.ak1.s,ig.ak2.s,ig.ak3.s,ig.ak4.s,ig.ak5.s,ig.ak6.s)

#### Calculate network metrics and build dataframe of observed characteristics ####

#Centrality (betweenness)
ig.ro1.r_between <- igraph::betweenness(ig.ro1.r,weights=E(ig.ro1.r),normalized=FALSE)
ig.ro2.r_between <- igraph::betweenness(ig.ro2.r,weights=E(ig.ro2.r),normalized=FALSE)
ig.ro3.r_between <- igraph::betweenness(ig.ro3.r,weights=E(ig.ro3.r),normalized=FALSE)
ig.ro4.r_between <- igraph::betweenness(ig.ro4.r,weights=E(ig.ro4.r),normalized=FALSE)
ig.ro5.r_between <- igraph::betweenness(ig.ro5.r,weights=E(ig.ro5.r),normalized=FALSE)
ig.ro6.r_between <- igraph::betweenness(ig.ro6.r,weights=E(ig.ro6.r),normalized=FALSE)
ig.ro1.s_between <- igraph::betweenness(ig.ro1.s,weights=E(ig.ro1.s),normalized=FALSE)
ig.ro2.s_between <- igraph::betweenness(ig.ro2.s,weights=E(ig.ro2.s),normalized=FALSE)
ig.ro3.s_between <- igraph::betweenness(ig.ro3.s,weights=E(ig.ro3.s),normalized=FALSE)
ig.ro4.s_between <- igraph::betweenness(ig.ro4.s,weights=E(ig.ro4.s),normalized=FALSE)
ig.ro5.s_between <- igraph::betweenness(ig.ro5.s,weights=E(ig.ro5.s),normalized=FALSE)
ig.ro6.s_between <- igraph::betweenness(ig.ro6.s,weights=E(ig.ro6.s),normalized=FALSE)
ig.ak1.r_between <- igraph::betweenness(ig.ak1.r,weights=E(ig.ak1.r),normalized=FALSE)
ig.ak2.r_between <- igraph::betweenness(ig.ak2.r,weights=E(ig.ak2.r),normalized=FALSE)
ig.ak3.r_between <- igraph::betweenness(ig.ak3.r,weights=E(ig.ak3.r),normalized=FALSE)
ig.ak4.r_between <- igraph::betweenness(ig.ak4.r,weights=E(ig.ak4.r),normalized=FALSE)
ig.ak5.r_between <- igraph::betweenness(ig.ak5.r,weights=E(ig.ak5.r),normalized=FALSE)
ig.ak6.r_between <- igraph::betweenness(ig.ak6.r,weights=E(ig.ak6.r),normalized=FALSE)
ig.ak1.s_between <- igraph::betweenness(ig.ak1.s,weights=E(ig.ak1.s),normalized=FALSE)
ig.ak2.s_between <- igraph::betweenness(ig.ak2.s,weights=E(ig.ak2.s),normalized=FALSE)
ig.ak3.s_between <- igraph::betweenness(ig.ak3.s,weights=E(ig.ak3.s),normalized=FALSE)
ig.ak4.s_between <- igraph::betweenness(ig.ak4.s,weights=E(ig.ak4.s),normalized=FALSE)
ig.ak5.s_between <- igraph::betweenness(ig.ak5.s,weights=E(ig.ak5.s),normalized=FALSE)
ig.ak6.s_between <- igraph::betweenness(ig.ak6.s,weights=E(ig.ak6.s),normalized=FALSE)

haka_centrality <- c(mean(ig.ro1.r_between),mean(ig.ro2.r_between),mean(ig.ro3.r_between),
                        mean(ig.ro4.r_between),mean(ig.ro5.r_between),mean(ig.ro6.r_between),
                        mean(ig.ro1.s_between),mean(ig.ro2.s_between),mean(ig.ro3.s_between),
                        mean(ig.ro4.s_between),mean(ig.ro5.s_between),mean(ig.ro6.s_between),
                        mean(ig.ak1.r_between),mean(ig.ak2.r_between),mean(ig.ak3.r_between),
                        mean(ig.ak4.r_between),mean(ig.ak5.r_between),mean(ig.ak6.r_between),
                        mean(ig.ak1.s_between),mean(ig.ak2.s_between),mean(ig.ak3.s_between),
                        mean(ig.ak4.s_between),mean(ig.ak5.s_between),mean(ig.ak6.s_between))
haka_centrality <- as.data.frame(haka_centrality)

#Connectedness (Degree)
ig.ro1.r_degree <- igraph::degree(ig.ro1.r, mode="all", normalized=TRUE)
ig.ro2.r_degree <- igraph::degree(ig.ro2.r, mode="all", normalized=TRUE)
ig.ro3.r_degree <- igraph::degree(ig.ro3.r, mode="all", normalized=TRUE)
ig.ro4.r_degree <- igraph::degree(ig.ro4.r, mode="all", normalized=TRUE)
ig.ro5.r_degree <- igraph::degree(ig.ro5.r, mode="all", normalized=TRUE)
ig.ro6.r_degree <- igraph::degree(ig.ro6.r, mode="all", normalized=TRUE)
ig.ro1.s_degree <- igraph::degree(ig.ro1.s, mode="all", normalized=TRUE)
ig.ro2.s_degree <- igraph::degree(ig.ro2.s, mode="all", normalized=TRUE)
ig.ro3.s_degree <- igraph::degree(ig.ro3.s, mode="all", normalized=TRUE)
ig.ro4.s_degree <- igraph::degree(ig.ro4.s, mode="all", normalized=TRUE)
ig.ro5.s_degree <- igraph::degree(ig.ro5.s, mode="all", normalized=TRUE)
ig.ro6.s_degree <- igraph::degree(ig.ro6.s, mode="all", normalized=TRUE)
ig.ak1.r_degree <- igraph::degree(ig.ak1.r, mode="all", normalized=TRUE)
ig.ak2.r_degree <- igraph::degree(ig.ak2.r, mode="all", normalized=TRUE)
ig.ak3.r_degree <- igraph::degree(ig.ak3.r, mode="all", normalized=TRUE)
ig.ak4.r_degree <- igraph::degree(ig.ak4.r, mode="all", normalized=TRUE)
ig.ak5.r_degree <- igraph::degree(ig.ak5.r, mode="all", normalized=TRUE)
ig.ak6.r_degree <- igraph::degree(ig.ak6.r, mode="all", normalized=TRUE)
ig.ak1.s_degree <- igraph::degree(ig.ak1.s, mode="all", normalized=TRUE)
ig.ak2.s_degree <- igraph::degree(ig.ak2.s, mode="all", normalized=TRUE)
ig.ak3.s_degree <- igraph::degree(ig.ak3.s, mode="all", normalized=TRUE)
ig.ak4.s_degree <- igraph::degree(ig.ak4.s, mode="all", normalized=TRUE)
ig.ak5.s_degree <- igraph::degree(ig.ak5.s, mode="all", normalized=TRUE)
ig.ak6.s_degree <- igraph::degree(ig.ak6.s, mode="all", normalized=TRUE)
ig.ak1.r_degree <- igraph::degree(ig.ak1.r, mode="all", normalized=TRUE)
ig.ak2.r_degree <- igraph::degree(ig.ak2.r, mode="all", normalized=TRUE)
ig.ak3.r_degree <- igraph::degree(ig.ak3.r, mode="all", normalized=TRUE)
ig.ak4.r_degree <- igraph::degree(ig.ak4.r, mode="all", normalized=TRUE)
ig.ak5.r_degree <- igraph::degree(ig.ak5.r, mode="all", normalized=TRUE)
ig.ak6.r_degree <- igraph::degree(ig.ak6.r, mode="all", normalized=TRUE)
ig.ak1.s_degree <- igraph::degree(ig.ak1.s, mode="all", normalized=TRUE)
ig.ak2.s_degree <- igraph::degree(ig.ak2.s, mode="all", normalized=TRUE)
ig.ak3.s_degree <- igraph::degree(ig.ak3.s, mode="all", normalized=TRUE)
ig.ak4.s_degree <- igraph::degree(ig.ak4.s, mode="all", normalized=TRUE)
ig.ak5.s_degree <- igraph::degree(ig.ak5.s, mode="all", normalized=TRUE)
ig.ak6.s_degree <- igraph::degree(ig.ak6.s, mode="all", normalized=TRUE)

haka_connectedness <- c(mean(ig.ro1.r_degree),mean(ig.ro2.r_degree),mean(ig.ro3.r_degree),
                     mean(ig.ro4.r_degree),mean(ig.ro5.r_degree),mean(ig.ro6.r_degree),
                     mean(ig.ro1.s_degree),mean(ig.ro2.s_degree),mean(ig.ro3.s_degree),
                     mean(ig.ro4.s_degree),mean(ig.ro5.s_degree),mean(ig.ro6.s_degree),
                     mean(ig.ak1.r_degree),mean(ig.ak2.r_degree),mean(ig.ak3.r_degree),
                     mean(ig.ak4.r_degree),mean(ig.ak5.r_degree),mean(ig.ak6.r_degree),
                     mean(ig.ak1.s_degree),mean(ig.ak2.s_degree),mean(ig.ak3.s_degree),
                     mean(ig.ak4.s_degree),mean(ig.ak5.s_degree),mean(ig.ak6.s_degree))
haka_connectedness <- as.data.frame(haka_connectedness)


#Density
haka_density <- c(edge_density(ig.ro1.r),edge_density(ig.ro2.r),edge_density(ig.ro3.r),
             edge_density(ig.ro4.r),edge_density(ig.ro5.r),edge_density(ig.ro6.r),
             edge_density(ig.ro1.s),edge_density(ig.ro2.s),edge_density(ig.ro3.s),
             edge_density(ig.ro4.s),edge_density(ig.ro5.s),edge_density(ig.ro6.s),
             edge_density(ig.ak1.r),edge_density(ig.ak2.r),edge_density(ig.ak3.r),
             edge_density(ig.ak4.r),edge_density(ig.ak5.r),edge_density(ig.ak6.r),
             edge_density(ig.ak1.s),edge_density(ig.ak2.s),edge_density(ig.ak3.s),
             edge_density(ig.ak4.s),edge_density(ig.ak5.s),edge_density(ig.ak6.s))
haka_density <- as.data.frame(haka_density)

#Dataframe building
plot<-as.data.frame(c("RO1","RO2","RO3","RO4","RO5","RO6","RO1","RO2","RO3","RO4","RO5","RO6",
                      "AK1","AK2","AK3","AK4","AK5","AK6","AK1","AK2","AK3","AK4","AK5","AK6"))
hab_type <- as.data.frame(rep(c("Remnant ohia","Afforested koa"),each=12))
sample_type<- as.data.frame(rep(c("roots","soil"),each=6,times=2))

fungal_networks <- cbind(plot,hab_type,sample_type,haka_centrality, haka_connectedness, haka_density)
colnames(fungal_networks) <- c("Plot","HabitatType","SampleType","Centrality","Connectedness","Density")

## Welch t-tests
#Connectedness (degree)
connectedness_between_hab_roots <- t.test(subset(fungal_networks,
                                    HabitatType == "Afforested koa" & SampleType == "roots")$Connectedness,
                     subset(fungal_networks,HabitatType == "Remnant ohia" & SampleType == "roots")$Connectedness,
                     paired=FALSE, var.equal=FALSE)
connectedness_between_hab_roots

connectedness_between_hab_soil <- t.test(subset(fungal_networks,
                                                 HabitatType == "Afforested koa" & SampleType == "soil")$Connectedness,
                                          subset(fungal_networks,
                                                 HabitatType == "Remnant ohia" & SampleType == "soil")$Connectedness,
                                          paired=FALSE, var.equal=FALSE)
connectedness_between_hab_soil

connectedness_within_RO <- t.test(subset(fungal_networks,
                                            HabitatType == "Remnant ohia" & SampleType == "roots")$Connectedness,
                                     subset(fungal_networks,
                                            HabitatType == "Remnant ohia" & SampleType == "soil")$Connectedness,
                                     paired=FALSE, var.equal=FALSE)
connectedness_within_RO

connectedness_within_AK <- t.test(subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "roots")$Connectedness,
                                  subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "soil")$Connectedness,
                                  paired=FALSE, var.equal=FALSE)
connectedness_within_AK

#Centrality (betweenness)
Centrality_between_hab_roots <- t.test(subset(fungal_networks,
                                                 HabitatType == "Afforested koa" & SampleType == "roots")$Centrality,
                                          subset(fungal_networks,HabitatType == "Remnant ohia" & SampleType == "roots")$Centrality,
                                          paired=FALSE, var.equal=FALSE)
Centrality_between_hab_roots

Centrality_between_hab_soil <- t.test(subset(fungal_networks,
                                                HabitatType == "Afforested koa" & SampleType == "soil")$Centrality,
                                         subset(fungal_networks,
                                                HabitatType == "Remnant ohia" & SampleType == "soil")$Centrality,
                                         paired=FALSE, var.equal=FALSE)
Centrality_between_hab_soil

Centrality_within_RO <- t.test(subset(fungal_networks,
                                         HabitatType == "Remnant ohia" & SampleType == "roots")$Centrality,
                                  subset(fungal_networks,
                                         HabitatType == "Remnant ohia" & SampleType == "soil")$Centrality,
                                  paired=FALSE, var.equal=FALSE)
Centrality_within_RO

Centrality_within_AK <- t.test(subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "roots")$Centrality,
                                  subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "soil")$Centrality,
                                  paired=FALSE, var.equal=FALSE)
Centrality_within_AK

#Density
Density_between_hab_roots <- t.test(subset(fungal_networks,
                                                 HabitatType == "Afforested koa" & SampleType == "roots")$Density,
                                          subset(fungal_networks,
                                                 HabitatType == "Remnant ohia" & SampleType == "roots")$Density,
                                          paired=FALSE, var.equal=FALSE)
Density_between_hab_roots

Density_between_hab_soil <- t.test(subset(fungal_networks,
                                                HabitatType == "Afforested koa" & SampleType == "soil")$Density,
                                         subset(fungal_networks,
                                                HabitatType == "Remnant ohia" & SampleType == "soil")$Density,
                                         paired=FALSE, var.equal=FALSE)
Density_between_hab_soil

Density_within_RO <- t.test(subset(fungal_networks,
                                         HabitatType == "Remnant ohia" & SampleType == "roots")$Density,
                                  subset(fungal_networks,
                                         HabitatType == "Remnant ohia" & SampleType == "soil")$Density,
                                  paired=FALSE, var.equal=FALSE)
Density_within_RO

Density_within_AK <- t.test(subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "roots")$Density,
                                  subset(fungal_networks,
                                         HabitatType == "Afforested koa" & SampleType == "soil")$Density,
                                  paired=FALSE, var.equal=FALSE)
Density_within_AK

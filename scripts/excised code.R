# removed code


hakalau_map_zoom <-get_map(location=c(-155.320,19.83),zoom=14,maptype="satellite",color="color")
haka_map <- ggmap(hakalau_map_zoom) +
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.text=element_text(colour="black",size=8)) +
  theme(axis.title=element_text(colour="black",size=8))
plot(haka_map)


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


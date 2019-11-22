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


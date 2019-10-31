############################################################################################################
### Composite co-occurrence networks (habitat type x sample type) for visual and identifying keystone sp ###
############################################################################################################

### SpiecEasi ###
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)
library(car)
library(huge)
library(MASS)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(scales)


############## Using phyloseq

### Subset phyloseq object into different hab types and then by soil--using virtual taxa
RO.physeq = subset_samples(haka_VT_soil_physeq,HabitatType == "Remnant Forest")
AK.physeq = subset_samples(haka_VT_soil_physeq,HabitatType == "Restored Forest")

### Construct networks using Spiec Easi (save in case R crashes as an RDS file)
se.RO<- spiec.easi(RO.physeq, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
saveRDS(se.RO,"output/spiec_easi_files/se.RO.RData")
se.RO<- readRDS("output/spiec_easi_files/se.RO.RData")

se.AK<- spiec.easi(AK.physeq, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
saveRDS(se.AK,"output/spiec_easi_files/se.AK.RData")
se.AK <- readRDS("output/spiec_easi_files/se.AK.RData")

### Convert to igraph objects for plotting and computing graph statistics.
#Bring in edge weights
ig.RO<-adj2igraph(se.RO$refit$stars, vertex.attr=list(name=taxa_names(RO.physeq)))
ig.AK<-adj2igraph(se.AK$refit$stars, vertex.attr=list(name=taxa_names(AK.physeq)))

V(ig.RO)$label<-V(ig.RO)$name
V(ig.RO)$frame.color = "black"

V(ig.AK)$label<-V(ig.AK)$name
V(ig.AK)$frame.color = "black"

# now both are ready for igraph and plotting


#########  RO network in igrpah
#Calculate number of samples each AM fungal species (node) is detected, and add as vertex attribute
RO_dat =fast_melt(RO.physeq)
RO_no.samples = RO_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO_no.samples = as.data.frame(RO_no.samples)
V(ig.RO)$no.samples = as.numeric(RO_no.samples$no.samples[match(V(ig.RO)$Name, RO_no.samples$TaxaID)])
ig.RO = delete_vertices(ig.RO, V(ig.RO)$no.samples == 0)

#Calculate relative abundance of each AM fungal species (node), and add as vertex attribute
RO.rel_abund <- summarize_taxa(RO.physeq,"Species")
V(ig.RO)$rel.abund = as.numeric(RO.rel_abund$meanRA[match(V(ig.RO)$name, RO_no.samples$TaxaID)])

#Important nodes
#First calcuate betweenness (nodes w/ high betweenness centrality represent key connector species in network)
V(ig.RO)$betweenness = igraph::betweenness(ig.RO,weights=E(ig.RO),normalized=FALSE)
#Then calculate node degree (nodes with high degree represent network hubs)
V(ig.RO)$degree = igraph::degree(ig.RO,mode="all",normalized=FALSE)
#Size nodes according to the two metrics
V(ig.RO)$size = (V(ig.RO)$betweenness*V(ig.RO)$degree)/100

#Colour by Family. First need to assign colours to match previously made graphs
net_tax=read.csv("data/haka_soil_taxonomy.csv", header = TRUE,row.names = 1)
net_tax$OTU<-row.names(net_tax)
V(ig.RO)$Family=as.character(net_tax$Family[match(V(ig.RO)$name, net_tax$OTU)])

# if want to set colors manually with names....
pal<-brewer.pal(10,"Paired")

V(ig.RO)$color=V(ig.RO)$Family
V(ig.RO)$color=gsub("Acaulosporaceae",pal[1] ,V(ig.RO)$color)
V(ig.RO)$color=gsub("Ambisporaceae",pal[2],V(ig.RO)$color)
V(ig.RO)$color=gsub("Archaeosporaceae",pal[3],V(ig.RO)$color)
V(ig.RO)$color=gsub("Archaeosporales_uncultured",pal[4],V(ig.RO)$color)
V(ig.RO)$color=gsub("Claroideoglomeraceae",pal[5],V(ig.RO)$color)
V(ig.RO)$color=gsub("Diversisporaceae",pal[6],V(ig.RO)$color)
V(ig.RO)$color=gsub("Geosiphonaceae",pal[7],V(ig.RO)$color)
V(ig.RO)$color=gsub("Gigasporaceae",pal[8],V(ig.RO)$color)
V(ig.RO)$color=gsub("Glomeraceae",pal[9],V(ig.RO)$color)
V(ig.RO)$color=gsub("Paraglomeraceae",pal[10],V(ig.RO)$color)

######### RO network plot
GroupRO <- levels(as.factor(V(ig.RO)$Family))
ColorRO <- levels(as.factor(V(ig.RO)$color))

par(mfrow=c(1,1))
par(mar=c(0,2,2,0))
set.seed(606)
plot(ig.RO, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.RO)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr, vertex.label=V(ig.RO)$Family)
legend('topright',legend=GroupRO, pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.6, inset=c(0.6, 0.025))

######### 


#########  AK network in igrpah
#Calculate number of samples each AM fungal species (node) is detected, and add as vertex attribute
AK_dat =fast_melt(AK.physeq)
AK_no.samples = AK_dat[, list(no.samples=sum(count>0)),by=TaxaID]
AK_no.samples = as.data.frame(AK_no.samples)
V(ig.AK)$no.samples = as.numeric(AK_no.samples$no.samples[match(V(ig.AK)$name, AK_no.samples$TaxaID)])
ig.AK = delete_vertices(ig.AK, V(ig.AK)$no.samples == 0)

#Calculate relative abundance of each AM fungal species (node), and add as vertex attribute
AK.rel_abund <- summarize_taxa(AK.physeq,"Species")
V(ig.AK)$rel.abund = as.numeric(AK.rel_abund$meanRA[match(V(ig.AK)$name, AK_no.samples$TaxaID)])

#Important nodes
#First calcuate betweenness (nodes w/ high betweenness centrality represent key connector species in network)
V(ig.AK)$betweenness = igraph::betweenness(ig.AK,weights=E(ig.AK),normalized=FALSE)
#Then calculate node degree (nodes with high degree represent network hubs)
V(ig.AK)$degree = igraph::degree(ig.AK,mode="all",normalized=FALSE)
#Size nodes according to the two metrics
V(ig.AK)$size = (V(ig.AK)$betweenness*V(ig.AK)$degree)/100

#Colour by Family
V(ig.AK)$Family=as.character(net_tax$Family[match(V(ig.AK)$name, net_tax$OTU)])

# if want to set colors manually with names....
V(ig.AK)$color=V(ig.AK)$Family
V(ig.AK)$color=gsub("Acaulosporaceae",pal[1] ,V(ig.AK)$color)
V(ig.AK)$color=gsub("Ambisporaceae",pal[2],V(ig.AK)$color)
V(ig.AK)$color=gsub("Archaeosporaceae",pal[3],V(ig.AK)$color)
V(ig.AK)$color=gsub("Archaeosporales_uncultured",pal[4],V(ig.AK)$color)
V(ig.AK)$color=gsub("Claroideoglomeraceae",pal[5],V(ig.AK)$color)
V(ig.AK)$color=gsub("Diversisporaceae",pal[6],V(ig.AK)$color)
V(ig.AK)$color=gsub("Geosiphonaceae",pal[7],V(ig.AK)$color)
V(ig.AK)$color=gsub("Gigasporaceae",pal[8],V(ig.AK)$color)
V(ig.AK)$color=gsub("Glomeraceae",pal[9],V(ig.AK)$color)
V(ig.AK)$color=gsub("Paraglomeraceae",pal[10],V(ig.AK)$color)

#########  AK network plot
GroupAK <- levels(as.factor(V(ig.AK)$Family))
ColorAK <- levels(as.factor(V(ig.AK)$color))

par(mfrow=c(1,1))
par(mar=c(0,2,2,0))
set.seed(606)
plot(ig.AK, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.AK)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=V(ig.AK)$Family)
legend('topright',legend=GroupAK, pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.6, inset=c(0.6, 0.025))

######### 


#Manuscript plot
par(mar=c(0.75,0.75,0.75,0.75))
par(mfrow=c(1,2))
set.seed(606)
plot(ig.RO, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.RO)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Remnant Forest")
set.seed(606)
plot(ig.AK, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.AK)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Restored Forest")
legend('topright', title="AM fungal Family",'right', legend=GroupAK, pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.5, inset=c(0.3, 0.08))

dev.copy(pdf, "figures/Network.pdf", height=6, width=8)
dev.off() 


######### Degree distribution plots
dd.RO     <- degree.distribution(ig.RO)
dd.AK     <- degree.distribution(ig.AK)

par(mfrow=c(1,1), mar=c(4,4,2,2))
# RO degree distribution
plot(0:(length(dd.RO)-1), dd.RO, col="#88A550" , type='b', ylab="Frequency", xlab="Degree", 
    main="Degree Distributions", ylim=c(0,0.5), xlim=c(0,12), lwd=1.5)
points(0:(length(dd.AK)-1), dd.AK, ylim=c(0,0.5), xlim=c(0,12), col="#336B87", type='b', lwd=1.5)
legend("topright", c("RO","AK"), col=c("#88A550","#336B87"), pch=1, lty=1)
dev.copy(pdf, "figures/dd.mb.habitat.pdf", height=4, width=5)
dev.off() 
######### 

##############################################################################################################
### Co-occurrence networks by plot (habitat type x sample type x plot) for overall network characteristics ###
##############################################################################################################

RO1.physeq = subset_samples(RO.physeq, HabitatType == "Remnant Forest" & Plot =="RO1")
RO2.physeq = subset_samples(RO.physeq, HabitatType == "Remnant Forest" & Plot =="RO2")
RO3.physeq = subset_samples(RO.physeq, HabitatType == "Remnant Forest" & Plot =="RO3")
RO4.physeq = subset_samples(RO.physeq, HabitatType == "Remnant Forest" & Plot =="RO4")
RO5.physeq = subset_samples(RO.physeq, HabitatType == "Remnant Forest" & Plot =="RO5")
RO6.physeq = subset_samples(RO.physeq, HabitatType == "Remnant Forest" & Plot =="RO6")

AK1.physeq = subset_samples(AK.physeq, HabitatType == "Restored Forest"  & Plot == "AK1")
AK2.physeq = subset_samples(AK.physeq, HabitatType == "Restored Forest"  & Plot == "AK2")
AK3.physeq = subset_samples(AK.physeq, HabitatType == "Restored Forest"  & Plot == "AK3")
AK4.physeq = subset_samples(AK.physeq, HabitatType == "Restored Forest"  & Plot == "AK4")
AK5.physeq = subset_samples(AK.physeq, HabitatType == "Restored Forest"  & Plot == "AK5")
AK6.physeq = subset_samples(AK.physeq, HabitatType == "Restored Forest"  & Plot == "AK6")


#### Make Spiec-Easi objects and SAVE ####
se.RO1<- spiec.easi(RO1.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO1,"output/spiec_easi_files/se.RO1.RData")
se.RO2<- spiec.easi(RO2.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO2,"output/spiec_easi_files/se.RO2.RData")
se.RO3<- spiec.easi(RO3.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO3,"output/spiec_easi_files/se.RO3.RData")
se.RO4<- spiec.easi(RO4.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO4,"output/spiec_easi_files/se.RO4.RData")
se.RO5<- spiec.easi(RO5.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO5,"output/spiec_easi_files/se.RO5.RData")
se.RO6<- spiec.easi(RO6.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO6,"output/spiec_easi_files/se.RO6.RData")

se.AK1<- spiec.easi(AK1.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK1,"output/spiec_easi_files/se.AK1.RData")
se.AK2<- spiec.easi(AK2.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK2,"output/spiec_easi_files/se.AK2.RData")
se.AK3<- spiec.easi(AK3.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK3,"output/spiec_easi_files/se.AK3.RData")
se.AK4<- spiec.easi(AK4.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK4,"output/spiec_easi_files/se.AK4.RData")
se.AK5<- spiec.easi(AK5.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK5,"output/spiec_easi_files/se.AK5.RData")
se.AK6<- spiec.easi(AK6.physeq,method='mb',lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK6,"output/spiec_easi_files/se.AK6.RData")

#### Bring networks back into R (in case R crashes like a mother fucker)    ####
se.RO1<- readRDS("output/spiec_easi_files/se.RO1.RData")
se.RO2<- readRDS("output/spiec_easi_files/se.RO2.RData")
se.RO3<- readRDS("output/spiec_easi_files/se.RO3.RData")
se.RO4<- readRDS("output/spiec_easi_files/se.RO4.RData")
se.RO5<- readRDS("output/spiec_easi_files/se.RO5.RData")
se.RO6<- readRDS("output/spiec_easi_files/se.RO6.RData")

se.AK1<- readRDS("output/spiec_easi_files/se.AK1.RData")
se.AK2<- readRDS("output/spiec_easi_files/se.AK2.RData")
se.AK3<- readRDS("output/spiec_easi_files/se.AK3.RData")
se.AK4<- readRDS("output/spiec_easi_files/se.AK4.RData")
se.AK5<- readRDS("output/spiec_easi_files/se.AK5.RData")
se.AK6<- readRDS("output/spiec_easi_files/se.AK6.RData")

#### Convert to igraph networks ####

ig.RO1<-adj2igraph(se.RO1$refit$stars,vertex.attr=list(name=taxa_names(RO1.physeq)))
RO1_dat =fast_melt(RO1.physeq)
RO1_no.samples = RO1_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO1_no.samples = as.data.frame(RO1_no.samples)
V(ig.RO1)$no.samples=as.numeric(RO1_no.samples$no.samples[match(V(ig.RO1)$name, RO1_no.samples$TaxaID)])
ig.RO1= delete_vertices(ig.RO1,V(ig.RO1)$no.samples == 0)

ig.RO2<-adj2igraph(se.RO2$refit$stars,vertex.attr=list(name=taxa_names(RO2.physeq)))
RO2_dat =fast_melt(RO2.physeq)
RO2_no.samples = RO2_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO2_no.samples = as.data.frame(RO2_no.samples)
V(ig.RO2)$no.samples=as.numeric(RO2_no.samples$no.samples[match(V(ig.RO2)$name, RO2_no.samples$TaxaID)])
ig.RO2= delete_vertices(ig.RO2,V(ig.RO2)$no.samples == 0)

ig.RO3<-adj2igraph(se.RO3$refit$stars,vertex.attr=list(name=taxa_names(RO3.physeq)))
RO3_dat =fast_melt(RO3.physeq)
RO3_no.samples = RO3_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO3_no.samples = as.data.frame(RO3_no.samples)
V(ig.RO3)$no.samples=as.numeric(RO3_no.samples$no.samples[match(V(ig.RO3)$name, RO3_no.samples$TaxaID)])
ig.RO3= delete_vertices(ig.RO3,V(ig.RO3)$no.samples == 0)

ig.RO4<-adj2igraph(se.RO4$refit$stars,vertex.attr=list(name=taxa_names(RO4.physeq)))
RO4_dat =fast_melt(RO4.physeq)
RO4_no.samples = RO4_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO4_no.samples = as.data.frame(RO4_no.samples)
V(ig.RO4)$no.samples=as.numeric(RO4_no.samples$no.samples[match(V(ig.RO4)$name, RO4_no.samples$TaxaID)])
ig.RO4= delete_vertices(ig.RO4,V(ig.RO4)$no.samples == 0)

ig.RO5<-adj2igraph(se.RO5$refit$stars,vertex.attr=list(name=taxa_names(RO5.physeq)))
RO5_dat =fast_melt(RO5.physeq)
RO5_no.samples = RO5_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO5_no.samples = as.data.frame(RO5_no.samples)
V(ig.RO5)$no.samples=as.numeric(RO5_no.samples$no.samples[match(V(ig.RO5)$name, RO5_no.samples$TaxaID)])
ig.RO5= delete_vertices(ig.RO5,V(ig.RO5)$no.samples == 0)

ig.RO6<-adj2igraph(se.RO6$refit$stars,vertex.attr=list(name=taxa_names(RO6.physeq)))
RO6_dat =fast_melt(RO6.physeq)
RO6_no.samples = RO6_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO6_no.samples = as.data.frame(RO6_no.samples)
V(ig.RO6)$no.samples=as.numeric(RO6_no.samples$no.samples[match(V(ig.RO6)$name, RO6_no.samples$TaxaID)])
ig.RO6= delete_vertices(ig.RO6,V(ig.RO6)$no.samples == 0)



### AK igraph
ig.AK1<-adj2igraph(se.AK1$refit$stars,vertex.attr=list(name=taxa_names(AK1.physeq)))
AK1_dat =fast_melt(AK1.physeq)
AK1_no.samples = AK1_dat[, list(no.samples=sum(count>0)),by=TaxaID]
AK1_no.samples = as.data.frame(AK1_no.samples)
V(ig.AK1)$no.samples=as.numeric(AK1_no.samples$no.samples[match(V(ig.AK1)$name, AK1_no.samples$TaxaID)])
ig.AK1= delete_vertices(ig.AK1,V(ig.AK1)$no.samples == 0)

ig.AK2<-adj2igraph(se.AK2$refit$stars,vertex.attr=list(name=taxa_names(AK2.physeq)))
AK2_dat =fast_melt(AK2.physeq)
AK2_no.samples = AK2_dat[, list(no.samples=sum(count>0)),by=TaxaID]
AK2_no.samples = as.data.frame(AK2_no.samples)
V(ig.AK2)$no.samples=as.numeric(AK2_no.samples$no.samples[match(V(ig.AK2)$name, AK2_no.samples$TaxaID)])
ig.AK2= delete_vertices(ig.AK2,V(ig.AK2)$no.samples == 0)

ig.AK3<-adj2igraph(se.AK3$refit$stars,vertex.attr=list(name=taxa_names(AK3.physeq)))
AK3_dat =fast_melt(AK3.physeq)
AK3_no.samples = AK3_dat[, list(no.samples=sum(count>0)),by=TaxaID]
AK3_no.samples = as.data.frame(AK3_no.samples)
V(ig.AK3)$no.samples=as.numeric(AK3_no.samples$no.samples[match(V(ig.AK3)$name, AK3_no.samples$TaxaID)])
ig.AK3= delete_vertices(ig.AK3,V(ig.AK3)$no.samples == 0)

ig.AK4<-adj2igraph(se.AK4$refit$stars,vertex.attr=list(name=taxa_names(AK4.physeq)))
AK4_dat =fast_melt(AK4.physeq)
AK4_no.samples = AK4_dat[, list(no.samples=sum(count>0)),by=TaxaID]
AK4_no.samples = as.data.frame(AK4_no.samples)
V(ig.AK4)$no.samples=as.numeric(AK4_no.samples$no.samples[match(V(ig.AK4)$name, AK4_no.samples$TaxaID)])
ig.AK4= delete_vertices(ig.AK4,V(ig.AK4)$no.samples == 0)

ig.AK5<-adj2igraph(se.AK5$refit$stars,vertex.attr=list(name=taxa_names(AK5.physeq)))
AK5_dat =fast_melt(AK5.physeq)
AK5_no.samples = AK5_dat[, list(no.samples=sum(count>0)),by=TaxaID]
AK5_no.samples = as.data.frame(AK5_no.samples)
V(ig.AK5)$no.samples=as.numeric(AK5_no.samples$no.samples[match(V(ig.AK5)$name, AK5_no.samples$TaxaID)])
ig.AK5= delete_vertices(ig.AK5,V(ig.AK5)$no.samples == 0)

ig.AK6<-adj2igraph(se.AK6$refit$stars,vertex.attr=list(name=taxa_names(AK6.physeq)))
AK6_dat =fast_melt(AK6.physeq)
AK6_no.samples = AK6_dat[, list(no.samples=sum(count>0)),by=TaxaID]
AK6_no.samples = as.data.frame(AK6_no.samples)
V(ig.AK6)$no.samples=as.numeric(AK6_no.samples$no.samples[match(V(ig.AK6)$name, AK6_no.samples$TaxaID)])
ig.AK6= delete_vertices(ig.AK6,V(ig.AK6)$no.samples == 0)


obs_fungal_networks <- as.list(ig.RO1,ig.RO2,ig.RO3,ig.RO4,ig.RO5,ig.RO6,ig.AK1,ig.AK2,ig.AK3,ig.AK4,ig.AK5,ig.AK6)



##############################
#####  CHARACTERISTICS   #####
##############################
#Betweenness centrality
ig.RO_between <- igraph::betweenness(ig.RO, weights=E(ig.RO),normalized=TRUE)
ig.AK_between <- igraph::betweenness(ig.AK,weights=E(ig.AK),normalized=TRUE)

#Connectedness (degree)
ig.RO1_degree <- igraph::degree(ig.RO1, mode="all", normalized=TRUE)
ig.RO2_degree <- igraph::degree(ig.RO2, mode="all", normalized=TRUE)
ig.RO3_degree <- igraph::degree(ig.RO3, mode="all", normalized=TRUE)
ig.RO4_degree <- igraph::degree(ig.RO4, mode="all", normalized=TRUE)


## Welch t-tests
#Centrality
centrality_between_hab <- t.test(ig.RO_between,ig.AK_between,paired=FALSE, var.equal=FALSE)
centrality_between_hab

centrality_within_RO <- t.test(ig.RO_between, ig.RO_between,paired=FALSE, var.equal=FALSE)
centrality_within_RO

centrality_within_AK <- t.test(ig.AK_between, ig.AK_between,paired=FALSE, var.equal=FALSE)
centrality_within_AK



##############################
#####   KEYSTONE SPECIES   ###
##############################

# Remnant Forest
RO.between <- as.data.frame(V(ig.RO)$betweenness,normalized=TRUE)
RO.degree <- as.data.frame(V(ig.RO)$degree,mode="all",normalized=TRUE)
RO.rel.abund <- as.data.frame(V(ig.RO)$rel.abund)
RO.fam <- as.data.frame(V(ig.RO)$Family)
RO.species <- as.data.frame(net_tax$Species[match(V(ig.RO)$name,net_tax$OTU)])
RO.no.samples <- as.data.frame(V(ig.RO)$no.samples)
sample_data(RO.physeq)
RO.total.samples <- as.data.frame(rep(237),times=162)
RO.keystone <- cbind(RO.between,RO.degree,RO.rel.abund,RO.fam,RO.species,RO.no.samples,RO.total.samples)
rownames(RO.keystone) <- V(ig.RO)$name
colnames(RO.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family","Species","No.Samples","TotalSamples")
RO.keystone$No.Samples <- as.numeric(RO.keystone$No.Samples)
RO.keystone$perc.samples <- RO.keystone$No.Samples / RO.keystone$TotalSamples
RO.keystone$Prevalence <- RO.keystone$perc.samples * RO.keystone$RelativeAbundance

RO.keystone.bet<-RO.keystone[(RO.keystone$Betweenness > 0),] # only samples with >0 betweeness


# Restored Forest
AK.between <- as.data.frame(V(ig.AK)$betweenness,normalized=TRUE)
AK.degree <- as.data.frame(V(ig.AK)$degree,mode="all",normalized=TRUE)
AK.rel.abund <- as.data.frame(V(ig.AK)$rel.abund)
AK.fam <- as.data.frame(V(ig.AK)$Family)
AK.species <- as.data.frame(net_tax$Species[match(V(ig.AK)$name,net_tax$OTU)])
AK.no.samples <- as.data.frame(V(ig.AK)$no.samples)
sample_data(AK.physeq)
AK.total.samples <- as.data.frame(rep(237),times=162)
AK.keystone <- cbind(AK.between,AK.degree,AK.rel.abund,AK.fam,AK.species,AK.no.samples,AK.total.samples)
rownames(AK.keystone) <- V(ig.AK)$name
colnames(AK.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family","Species","No.Samples","TotalSamples")
AK.keystone$No.Samples <- as.numeric(AK.keystone$No.Samples)
AK.keystone$perc.samples <- AK.keystone$No.Samples / AK.keystone$TotalSamples
AK.keystone$Prevalence <- AK.keystone$perc.samples * AK.keystone$RelativeAbundance

AK.keystone.bet<-AK.keystone[(AK.keystone$Betweenness > 0),] # only samples with >0 betweeness


# plot RO
RO.keystone.plot <- ggplot(RO.keystone,aes(x=Degree,y=Betweenness)) +
    geom_point(aes(size=Prevalence,colour=Family),position="jitter") +
    ylab("log(Betweenness Centrality-normalized)") + 
    xlab("Node degree") +
    scale_x_continuous(limits=c(0,12), breaks=seq(0,12, by=3)) +
    scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x)10^x),
                       labels=trans_format("log10",math_format(10^.x)),limits=c(10^0,10^4)) +
    theme(text=element_text(colour="black",size=12)) + 
    theme(axis.text.x=element_text(hjust=1,colour="black",size=12)) +
    theme(axis.text.y=element_text(colour="black",size=12)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background=element_blank()) +
    scale_colour_manual(values=pal, limits=GroupRO) +
    guides(colour=guide_legend(override.aes=list(size=3))) +
    theme(legend.key=element_blank())+ theme(legend.text=element_text(size=8)) +
    theme(legend.key.size = unit(.4, "cm"))+
    ggtitle("A) Remnant Forest")
RO.keystone.plot

# plot AK
AK.keystone.plot <- ggplot(AK.keystone,aes(x=Degree,y=Betweenness)) +
    geom_point(aes(size=Prevalence,colour=Family),position="jitter") +
    ylab("") + 
    xlab("Node degree") +
    scale_x_continuous(limits=c(0,12), breaks=seq(0,12, by=3)) +
    scale_y_continuous(trans=log10_trans(),breaks=trans_breaks("log10",function(x)10^x),
                       labels=trans_format("log10",math_format(10^.x)),limits=c(10^0,10^4)) +
    theme(text=element_text(colour="black",size=12)) + 
    theme(axis.text.x=element_text(hjust=1,colour="black",size=12)) +
    theme(axis.text.y=element_text(colour="black",size=12)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background=element_blank()) +
    scale_colour_manual(values=pal, limits=GroupRO) +
    guides(colour=guide_legend(override.aes=list(size=3))) +
    theme(legend.key=element_blank())+ theme(legend.text=element_text(size=8)) +
    theme(legend.key.size = unit(.4, "cm"))+
    ggtitle("B) Restored Forest")
AK.keystone.plot

# extract legend from plot
KS.legend <- get_legend(
    # create some space to the left of the legend
    AK.keystone.plot + theme(legend.box.margin = margin(0, 0, 0, 12)))

KS.plots<- plot_grid(RO.keystone.plot + theme(legend.position = "none"), 
                          AK.keystone.plot + theme(legend.position = "none"), 
                          ncol=2,nrow=1)
plot_grid(KS.plots, KS.legend, rel_widths = c(2, 1)) # legend column 1/2 size as first object
dev.copy()
ggsave("figures/keystone.fig.pdf", width = 7, height = 4)










####### just using ESVs, subsampled
### Subset phyloseq object into different hab types--using ESVs

#First create a vector of otu file rownames and call it "keep", Then subset sample file using the keep vector


# for RO
RO.physeq.E = subset_samples(haka_soil_physeq,HabitatType == "Remnant Forest")
RO.physeq.ESV = filter_taxa(RO.physeq.E, function(x) sum(x > 2) > (0.20*length(x)), TRUE) # remove samples not seen at least 2 times and in 30% of samples

# for AK
AK.physeq.E = subset_samples(haka_soil_physeq,HabitatType == "Restored Forest")
AK.physeq.ESV = filter_taxa(AK.physeq.E, function(x) sum(x > 2) > (0.20*length(x)), TRUE) # remove samples not seen at least 2 times and in 30% of samples

######## Run models

### Construct networks using Spiec Easi (save in case R crashes as an RDS file)
se.RO.ESV<- spiec.easi(RO.physeq.ESV, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
se.AK.ESV<- spiec.easi(AK.physeq.ESV, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))

### Convert to igraph objects for plotting and computing graph statistics.
#Bring in edge weights
ig.RO.ESV<-adj2igraph(se.RO.ESV$refit$stars, vertex.attr=list(name=taxa_names(RO.physeq.ESV)))
ig.AK.ESV<-adj2igraph(se.AK.ESV$refit$stars, vertex.attr=list(name=taxa_names(AK.physeq.ESV)))

V(ig.RO.ESV)$label<-V(ig.RO.ESV)$name
V(ig.RO.ESV)$frame.color = "black"

V(ig.AK.ESV)$label<-V(ig.AK.ESV)$name
V(ig.AK.ESV)$frame.color = "black"


#########  RO network in igrpah
#Calculate number of samples each AM fungal species (node) is detected, and add as vertex attribute
RO_dat =fast_melt(RO.physeq.ESV)
RO_no.samples = RO_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO_no.samples = as.data.frame(RO_no.samples)
V(ig.RO.ESV)$no.samples = as.numeric(RO_no.samples$no.samples[match(V(ig.RO.ESV)$name, RO_no.samples$TaxaID)])
ig.RO.ESV = delete_vertices(ig.RO.ESV, V(ig.RO.ESV)$no.samples == 0)

#Calculate relative abundance of each AM fungal species (node), and add as vertex attribute
RO.rel_abund <- summarize_taxa(RO.physeq.ESV,"Species")
V(ig.RO)$rel.abund = as.numeric(RO.rel_abund$meanRA[match(V(ig.RO)$name, RO_no.samples$TaxaID)])

#Important nodes
#First calcuate betweenness (nodes w/ high betweenness centrality represent key connector species in network)
V(ig.RO)$betweenness = igraph::betweenness(ig.RO,weights=E(ig.RO),normalized=FALSE)
#Then calculate node degree (nodes with high degree represent network hubs)
V(ig.RO)$degree = igraph::degree(ig.RO,mode="all",normalized=FALSE)
#Size nodes according to the two metrics
V(ig.RO)$size = (V(ig.RO)$betweenness*V(ig.RO)$degree)/100

#Colour by Family. First need to assign colours to match previously made graphs
net_tax=read.csv("data/haka_soil_taxonomy.csv", header = TRUE,row.names = 1)
net_tax$OTU<-row.names(net_tax)
V(ig.RO)$Family=as.character(net_tax$Family[match(V(ig.RO)$name, net_tax$OTU)])










# COMES LATER
connectedness_within_AK <- t.test(subset(obs_fungal_networks,
                                         HabitatType == "Restored Forest")$Connectedness,
                                  subset(obs_fungal_networks,
                                         HabitatType == "Restored Forest")$Connectedness,
                                  paired=FALSE, var.equal=FALSE)
connectedness_within_AK

#Density
Density_between_hab_roots <- t.test(subset(fungal_networks,
                                           HabitatType == "Restored Forest" & SampleType == "roots")$Density,
                                    subset(fungal_networks,
                                           HabitatType == "Remnant Forest" & SampleType == "roots")$Density,
                                    paired=FALSE, var.equal=FALSE)
Density_between_hab_roots

Density_between_hab_soil <- t.test(subset(fungal_networks,
                                          HabitatType == "Restored Forest" & SampleType == "soil")$Density,
                                   subset(fungal_networks,
                                          HabitatType == "Remnant Forest" & SampleType == "soil")$Density,
                                   paired=FALSE, var.equal=FALSE)
Density_between_hab_soil

Density_within_RO <- t.test(subset(fungal_networks,
                                   HabitatType == "Remnant Forest" & SampleType == "roots")$Density,
                            subset(fungal_networks,
                                   HabitatType == "Remnant Forest" & SampleType == "soil")$Density,
                            paired=FALSE, var.equal=FALSE)
Density_within_RO

Density_within_AK <- t.test(subset(fungal_networks,
                                   HabitatType == "Restored Forest" & SampleType == "roots")$Density,
                            subset(fungal_networks,
                                   HabitatType == "Restored Forest" & SampleType == "soil")$Density,
                            paired=FALSE, var.equal=FALSE)
Density_within_AK



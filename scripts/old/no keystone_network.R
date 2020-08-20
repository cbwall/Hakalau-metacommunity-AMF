###################
###########
############
###########
#######
## test network topology without "keystone hubs"


### Subset phyloseq object into different hab types and then by soil--using virtual taxa
# remnant RO.physeq without keystone
RO.physeq.nokey<-subset_taxa(RO.physeq, Species != "VTX00272")

# restored AK.physeq without keystone
AK.physeq.nokey<-subset_taxa(AK.physeq, Species != "VTX00272")

### Construct networks using Spiec Easi (save in case R crashes as an RDS file)
se.RO.nokey<- spiec.easi(RO.physeq.nokey, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
saveRDS(se.RO.nokey,"output/spiec_easi_files/se.RO.nokey.RData")
se.RO.nokey<- readRDS("output/spiec_easi_files/se.RO.nokey.RData")

se.AK.nokey<- spiec.easi(AK.physeq.nokey, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=100))
saveRDS(se.AK.nokey,"output/spiec_easi_files/se.AK.nokey.RData")
se.AK.nokey <- readRDS("output/spiec_easi_files/se.AK.nokey.RData")

### Convert to igraph objects for plotting and computing graph statistics.
ig.RO.nokey<-adj2igraph(se.RO.nokey$refit$stars, vertex.attr=list(name=taxa_names(RO.physeq.nokey)))
ig.AK.nokey<-adj2igraph(se.AK.nokey$refit$stars, vertex.attr=list(name=taxa_names(AK.physeq.nokey)))

V(ig.RO.nokey)$label<-V(ig.RO.nokey)$name
V(ig.RO.nokey)$frame.color = "black"

V(ig.AK.nokey)$label<-V(ig.RO.nokey)$name
V(ig.AK.nokey)$frame.color = "black"

#########  RO network in igrpah
#Calculate number of samples each AM fungal species (node) is detected, and add as vertex attribute
RO_dat.nokey =fast_melt(RO.physeq.nokey)
RO_no.samples.nokey = RO_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
RO_no.samples.nokey = as.data.frame(RO_no.samples.nokey)
V(ig.RO.nokey)$no.samples = as.numeric(RO_no.samples.nokey$no.samples[match(V(ig.RO.nokey)$name, RO_no.samples.nokey$TaxaID)])
ig.RO.nokey = delete_vertices(ig.RO.nokey, V(ig.RO.nokey)$no.samples == 0)

#Calculate relative abundance of each AM fungal species (node), and add as vertex attribute
RO.rel_abund.nokey <- summarize_taxa(RO.physeq.nokey,"Species")
V(ig.RO.nokey)$rel.abund = as.numeric(RO.rel_abund.nokey$meanRA[match(V(ig.RO.nokey)$name, RO_no.samples.nokey$TaxaID)])

#Important nodes
#First calcuate betweenness (nodes w/ high betweenness centrality represent key connector species in network)
V(ig.RO.nokey)$betweenness = igraph::betweenness(ig.RO.nokey,weights=E(ig.RO.nokey),normalized=FALSE)
#Then calculate node degree (nodes with high degree represent network hubs)
V(ig.RO.nokey)$degree = igraph::degree(ig.RO.nokey,mode="all",normalized=FALSE)
#Size nodes according to the two metrics
V(ig.RO.nokey)$size = (V(ig.RO.nokey)$betweenness*V(ig.RO.nokey)$degree)/100

#Colour by Family. First need to assign colours to match previously made graphs
net_tax.nokey=read.csv("data/haka_soil_taxonomy.csv", header = TRUE,row.names = 1)
net_tax.nokey$OTU<-row.names(net_tax.nokey)
V(ig.RO.nokey)$Family=as.character(net_tax.nokey$Family[match(V(ig.RO.nokey)$name, net_tax.nokey$OTU)])
V(ig.RO.nokey)$Genus=as.character(net_tax.nokey$Genus[match(V(ig.RO.nokey)$name, net_tax.nokey$OTU)])

# if want to set colors manually with names....
pal<-brewer.pal(10,"Paired")

V(ig.RO.nokey)$color=V(ig.RO.nokey)$Family
V(ig.RO.nokey)$color=gsub("Acaulosporaceae",pal[1] ,V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Ambisporaceae",pal[2],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Archaeosporaceae",pal[3],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Archaeosporales_uncultured",pal[4],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Claroideoglomeraceae",pal[5],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Diversisporaceae",pal[6],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Geosiphonaceae",pal[7],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Gigasporaceae",pal[8],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Glomeraceae",pal[9],V(ig.RO.nokey)$color)
V(ig.RO.nokey)$color=gsub("Paraglomeraceae",pal[10],V(ig.RO.nokey)$color)

######### RO network plot
GroupRO.nokey <- levels(as.factor(V(ig.RO.nokey)$Family))
ColorRO.nokey <- levels(as.factor(V(ig.RO.nokey)$color))

par(mfrow=c(1,1))
par(mar=c(0,2,2,0))
set.seed(606)
plot(ig.RO.nokey, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.RO.nokey)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr, vertex.label=V(ig.RO.nokey)$Family)
legend('topright',legend=GroupRO.nokey, pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.6, inset=c(0.6, 0.025))

######### 


#########  AK network in igrpah
#Calculate number of samples each AM fungal species (node) is detected, and add as vertex attribute
AK_dat.nokey =fast_melt(AK.physeq.nokey)
AK_no.samples.nokey = AK_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
AK_no.samples.nokey = as.data.frame(AK_no.samples.nokey)
V(ig.AK.nokey)$no.samples = as.numeric(AK_no.samples.nokey$no.samples[match(V(ig.AK.nokey)$name, AK_no.samples.nokey$TaxaID)])
ig.AK.nokey = delete_vertices(ig.AK.nokey, V(ig.AK.nokey)$no.samples == 0)

#Calculate relative abundance of each AM fungal species (node), and add as vertex attribute
AK.rel_abund.nokey <- summarize_taxa(AK.physeq.nokey,"Species")
V(ig.AK.nokey)$rel.abund = as.numeric(AK.rel_abund.nokey$meanRA[match(V(ig.AK.nokey)$name, AK_no.samples.nokey$TaxaID)])

#Important nodes
#First calcuate betweenness (nodes w/ high betweenness centrality represent key connector species in network)
V(ig.AK.nokey)$betweenness = igraph::betweenness(ig.AK.nokey,weights=E(ig.AK.nokey),normalized=FALSE)
#Then calculate node degree (nodes with high degree represent network hubs)
V(ig.AK.nokey)$degree = igraph::degree(ig.AK.nokey,mode="all",normalized=FALSE)
#Size nodes according to the two metrics
V(ig.AK.nokey)$size = (V(ig.AK.nokey)$betweenness*V(ig.AK.nokey)$degree)/100

#Colour by Family
V(ig.AK.nokey)$Family=as.character(net_tax.nokey$Family[match(V(ig.AK.nokey)$name, net_tax.nokey$OTU)])
V(ig.AK.nokey)$Genus=as.character(net_tax.nokey$Genus[match(V(ig.AK.nokey)$name, net_tax.nokey$OTU)]) # add Genus to datafram

# if want to set colors manually with names....
V(ig.AK.nokey)$color=V(ig.AK.nokey)$Family
V(ig.AK.nokey)$color=gsub("Acaulosporaceae",pal[1] ,V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Ambisporaceae",pal[2],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Archaeosporaceae",pal[3],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Archaeosporales_uncultured",pal[4],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Claroideoglomeraceae",pal[5],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Diversisporaceae",pal[6],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Geosiphonaceae",pal[7],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Gigasporaceae",pal[8],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Glomeraceae",pal[9],V(ig.AK.nokey)$color)
V(ig.AK.nokey)$color=gsub("Paraglomeraceae",pal[10],V(ig.AK.nokey)$color)

#########  AK network plot
GroupAK.nokey <- levels(as.factor(V(ig.AK.nokey)$Family))
ColorAK.nokey <- levels(as.factor(V(ig.AK.nokey)$color))

par(mfrow=c(1,1))
par(mar=c(0,2,2,0))
set.seed(606)
plot(ig.AK.nokey, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.AK.nokey)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=V(ig.AK.nokey)$Family)
legend('topright',legend=GroupAK.nokey, pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.6, inset=c(0.6, 0.025))

# combined plot
par(mar=c(0.75,0.75,0.75,0.75))
par(mfrow=c(1,2))
set.seed(606)
plot(ig.RO.nokey, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.RO.nokey)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Remnant Forest")
set.seed(606)
plot(ig.AK.nokey, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.AK.nokey)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Restored Forest")
legend('topright', title="AM fungal Family",'right', legend=GroupAK.nokey, pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.5, inset=c(0.0, 0.1))


dev.copy(pdf, "figures/Network.nokeystone.pdf", height=7, width=8)
dev.off() 
###################
###########
############
###########
#######


##############################################################################################################
### Co-occurrence networks by plot (habitat type x plot) for overall network characteristics ###
##############################################################################################################

RO1.physeq.nokey = subset_samples(RO.physeq.nokey, HabitatType == "Remnant Forest" & Plot =="RO1")
RO2.physeq.nokey = subset_samples(RO.physeq.nokey, HabitatType == "Remnant Forest" & Plot =="RO2")
RO3.physeq.nokey = subset_samples(RO.physeq.nokey, HabitatType == "Remnant Forest" & Plot =="RO3")
RO4.physeq.nokey = subset_samples(RO.physeq.nokey, HabitatType == "Remnant Forest" & Plot =="RO4")
RO5.physeq.nokey = subset_samples(RO.physeq.nokey, HabitatType == "Remnant Forest" & Plot =="RO5")
RO6.physeq.nokey = subset_samples(RO.physeq.nokey, HabitatType == "Remnant Forest" & Plot =="RO6")

AK1.physeq.nokey = subset_samples(AK.physeq.nokey, HabitatType == "Restored Forest"  & Plot == "AK1")
AK2.physeq.nokey = subset_samples(AK.physeq.nokey, HabitatType == "Restored Forest"  & Plot == "AK2")
AK3.physeq.nokey = subset_samples(AK.physeq.nokey, HabitatType == "Restored Forest"  & Plot == "AK3")
AK4.physeq.nokey = subset_samples(AK.physeq.nokey, HabitatType == "Restored Forest"  & Plot == "AK4")
AK5.physeq.nokey = subset_samples(AK.physeq.nokey, HabitatType == "Restored Forest"  & Plot == "AK5")
AK6.physeq.nokey = subset_samples(AK.physeq.nokey, HabitatType == "Restored Forest"  & Plot == "AK6")


#### Make Spiec-Easi objects and SAVE ####
se.RO1.nokey<- spiec.easi(RO1.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO1.nokey,"output/spiec_easi_files/se.RO1.nokey.RData")
se.RO2.nokey<- spiec.easi(RO2.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO2.nokey,"output/spiec_easi_files/se.RO2.nokey.RData")
se.RO3.nokey<- spiec.easi(RO3.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO3.nokey,"output/spiec_easi_files/se.RO3.nokey.RData")
se.RO4.nokey<- spiec.easi(RO4.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO4.nokey,"output/spiec_easi_files/se.RO4.nokey.RData")
se.RO5.nokey<- spiec.easi(RO5.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO5.nokey,"output/spiec_easi_files/se.RO5.nokey.RData")
se.RO6.nokey<- spiec.easi(RO6.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.RO6.nokey,"output/spiec_easi_files/se.RO6.nokey.RData")


se.AK1.nokey<- spiec.easi(AK1.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK1.nokey,"output/spiec_easi_files/se.AK1.nokey.RData")
se.AK2.nokey<- spiec.easi(AK2.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK2.nokey,"output/spiec_easi_files/se.AK2.nokey.RData")
se.AK3.nokey<- spiec.easi(AK3.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK3.nokey,"output/spiec_easi_files/se.AK3.nokey.RData")
se.AK4.nokey<- spiec.easi(AK4.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK4.nokey,"output/spiec_easi_files/se.AK4.nokey.RData")
se.AK5.nokey<- spiec.easi(AK5.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK5.nokey,"output/spiec_easi_files/se.AK5.nokey.RData")
se.AK6.nokey<- spiec.easi(AK6.physeq.nokey,method='mb',
                          lambda.min.ratio=1e-2,nlambda=20,pulsar.params=list(rep.num=100))
saveRDS(se.AK6.nokey,"output/spiec_easi_files/se.AK6.nokey.RData")

#### Bring networks back into R (in case R crashes like a mother fucker)    ####
se.RO1.nokey<- readRDS("output/spiec_easi_files/se.RO1.nokey.RData")
se.RO2.nokey<- readRDS("output/spiec_easi_files/se.RO2.nokey.RData")
se.RO3.nokey<- readRDS("output/spiec_easi_files/se.RO3.nokey.RData")
se.RO4.nokey<- readRDS("output/spiec_easi_files/se.RO4.nokey.RData")
se.RO5.nokey<- readRDS("output/spiec_easi_files/se.RO5.nokey.RData")
se.RO6.nokey<- readRDS("output/spiec_easi_files/se.RO6.nokey.RData")

se.AK1.nokey<- readRDS("output/spiec_easi_files/se.AK1.nokey.RData")
se.AK2.nokey<- readRDS("output/spiec_easi_files/se.AK2.nokey.RData")
se.AK3.nokey<- readRDS("output/spiec_easi_files/se.AK3.nokey.RData")
se.AK4.nokey<- readRDS("output/spiec_easi_files/se.AK4.nokey.RData")
se.AK5.nokey<- readRDS("output/spiec_easi_files/se.AK5.nokey.RData")
se.AK6.nokey<- readRDS("output/spiec_easi_files/se.AK6.nokey.RData")

#### Convert to igraph networks ####

ig.RO1.nokey<-adj2igraph(se.RO1.nokey$refit$stars,vertex.attr=list(name=taxa_names(RO1.physeq.nokey)))
RO1_dat.nokey =fast_melt(RO1.physeq.nokey)
RO1_no.samples.nokey = RO1_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
RO1_no.samples.nokey = as.data.frame(RO1_no.samples.nokey)
V(ig.RO1.nokey)$no.samples=as.numeric(RO1_no.samples.nokey$no.samples[match(V(ig.RO1.nokey)$name, RO1_no.samples.nokey$TaxaID)])
ig.RO1.nokey= delete_vertices(ig.RO1.nokey,V(ig.RO1.nokey)$no.samples == 0)

ig.RO2.nokey<-adj2igraph(se.RO2.nokey$refit$stars,vertex.attr=list(name=taxa_names(RO2.physeq.nokey)))
RO2_dat.nokey =fast_melt(RO2.physeq.nokey)
RO2_no.samples.nokey = RO2_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
RO2_no.samples.nokey = as.data.frame(RO2_no.samples.nokey)
V(ig.RO2.nokey)$no.samples=as.numeric(RO2_no.samples.nokey$no.samples[match(V(ig.RO2.nokey)$name, RO2_no.samples.nokey$TaxaID)])
ig.RO2.nokey= delete_vertices(ig.RO2.nokey,V(ig.RO2.nokey)$no.samples == 0)

ig.RO3.nokey<-adj2igraph(se.RO3.nokey$refit$stars,vertex.attr=list(name=taxa_names(RO3.physeq.nokey)))
RO3_dat.nokey =fast_melt(RO3.physeq.nokey)
RO3_no.samples.nokey = RO3_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
RO3_no.samples.nokey = as.data.frame(RO3_no.samples.nokey)
V(ig.RO3.nokey)$no.samples=as.numeric(RO3_no.samples.nokey$no.samples[match(V(ig.RO3.nokey)$name, RO3_no.samples.nokey$TaxaID)])
ig.RO3.nokey= delete_vertices(ig.RO3.nokey,V(ig.RO3.nokey)$no.samples == 0)

ig.RO4.nokey<-adj2igraph(se.RO4.nokey$refit$stars,vertex.attr=list(name=taxa_names(RO4.physeq.nokey)))
RO4_dat.nokey =fast_melt(RO4.physeq.nokey)
RO4_no.samples.nokey = RO4_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
RO4_no.samples.nokey = as.data.frame(RO4_no.samples.nokey)
V(ig.RO4.nokey)$no.samples=as.numeric(RO4_no.samples.nokey$no.samples[match(V(ig.RO4.nokey)$name, RO4_no.samples.nokey$TaxaID)])
ig.RO4.nokey= delete_vertices(ig.RO4.nokey,V(ig.RO4.nokey)$no.samples == 0)

ig.RO5.nokey<-adj2igraph(se.RO5.nokey$refit$stars,vertex.attr=list(name=taxa_names(RO5.physeq.nokey)))
RO5_dat.nokey =fast_melt(RO5.physeq.nokey)
RO5_no.samples.nokey = RO5_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
RO5_no.samples.nokey = as.data.frame(RO5_no.samples.nokey)
V(ig.RO5.nokey)$no.samples=as.numeric(RO5_no.samples.nokey$no.samples[match(V(ig.RO5.nokey)$name, RO5_no.samples.nokey$TaxaID)])
ig.RO5.nokey= delete_vertices(ig.RO5.nokey,V(ig.RO5.nokey)$no.samples == 0)

ig.RO6.nokey<-adj2igraph(se.RO6.nokey$refit$stars,vertex.attr=list(name=taxa_names(RO6.physeq.nokey)))
RO6_dat.nokey =fast_melt(RO6.physeq.nokey)
RO6_no.samples.nokey = RO6_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
RO6_no.samples.nokey = as.data.frame(RO6_no.samples.nokey)
V(ig.RO6.nokey)$no.samples=as.numeric(RO6_no.samples.nokey$no.samples[match(V(ig.RO6.nokey)$name, RO6_no.samples.nokey$TaxaID)])
ig.RO6.nokey= delete_vertices(ig.RO6.nokey,V(ig.RO6.nokey)$no.samples == 0)



### AK igraph
ig.AK1.nokey<-adj2igraph(se.AK1.nokey$refit$stars,vertex.attr=list(name=taxa_names(AK1.physeq.nokey)))
AK1_dat.nokey =fast_melt(AK1.physeq.nokey)
AK1_no.samples.nokey = AK1_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
AK1_no.samples.nokey = as.data.frame(AK1_no.samples.nokey)
V(ig.AK1.nokey)$no.samples=as.numeric(AK1_no.samples.nokey$no.samples[match(V(ig.AK1.nokey)$name, AK1_no.samples.nokey$TaxaID)])
ig.AK1.nokey= delete_vertices(ig.AK1.nokey,V(ig.AK1.nokey)$no.samples == 0)

ig.AK2.nokey<-adj2igraph(se.AK2.nokey$refit$stars,vertex.attr=list(name=taxa_names(AK2.physeq.nokey)))
AK2_dat.nokey =fast_melt(AK2.physeq.nokey)
AK2_no.samples.nokey = AK2_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
AK2_no.samples.nokey = as.data.frame(AK2_no.samples.nokey)
V(ig.AK2.nokey)$no.samples=as.numeric(AK2_no.samples.nokey$no.samples[match(V(ig.AK2.nokey)$name, AK2_no.samples.nokey$TaxaID)])
ig.AK2.nokey= delete_vertices(ig.AK2.nokey,V(ig.AK2.nokey)$no.samples == 0)

ig.AK3.nokey<-adj2igraph(se.AK3.nokey$refit$stars,vertex.attr=list(name=taxa_names(AK3.physeq.nokey)))
AK3_dat.nokey =fast_melt(AK3.physeq.nokey)
AK3_no.samples.nokey = AK3_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
AK3_no.samples.nokey = as.data.frame(AK3_no.samples.nokey)
V(ig.AK3.nokey)$no.samples=as.numeric(AK3_no.samples.nokey$no.samples[match(V(ig.AK3.nokey)$name, AK3_no.samples.nokey$TaxaID)])
ig.AK3.nokey= delete_vertices(ig.AK3.nokey,V(ig.AK3.nokey)$no.samples == 0)

ig.AK4.nokey<-adj2igraph(se.AK4.nokey$refit$stars,vertex.attr=list(name=taxa_names(AK4.physeq.nokey)))
AK4_dat.nokey =fast_melt(AK4.physeq.nokey)
AK4_no.samples.nokey = AK4_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
AK4_no.samples.nokey = as.data.frame(AK4_no.samples.nokey)
V(ig.AK4.nokey)$no.samples=as.numeric(AK4_no.samples.nokey$no.samples[match(V(ig.AK4.nokey)$name, AK4_no.samples.nokey$TaxaID)])
ig.AK4.nokey= delete_vertices(ig.AK4.nokey,V(ig.AK4.nokey)$no.samples == 0)

ig.AK5.nokey<-adj2igraph(se.AK5.nokey$refit$stars,vertex.attr=list(name=taxa_names(AK5.physeq.nokey)))
AK5_dat.nokey =fast_melt(AK5.physeq.nokey)
AK5_no.samples.nokey = AK5_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
AK5_no.samples.nokey = as.data.frame(AK5_no.samples.nokey)
V(ig.AK5.nokey)$no.samples=as.numeric(AK5_no.samples.nokey$no.samples[match(V(ig.AK5.nokey)$name, AK5_no.samples.nokey$TaxaID)])
ig.AK5.nokey= delete_vertices(ig.AK5.nokey,V(ig.AK5.nokey)$no.samples == 0)

ig.AK6.nokey<-adj2igraph(se.AK6.nokey$refit$stars,vertex.attr=list(name=taxa_names(AK6.physeq.nokey)))
AK6_dat.nokey =fast_melt(AK6.physeq.nokey)
AK6_no.samples.nokey = AK6_dat.nokey[, list(no.samples=sum(count>0)),by=TaxaID]
AK6_no.samples.nokey = as.data.frame(AK6_no.samples.nokey)
V(ig.AK6.nokey)$no.samples=as.numeric(AK6_no.samples.nokey$no.samples[match(V(ig.AK6.nokey)$name, AK6_no.samples.nokey$TaxaID)])
ig.AK6.nokey= delete_vertices(ig.AK6.nokey,V(ig.AK6.nokey)$no.samples == 0)

obs_fungal_networks <- as.list(ig.RO1.nokey,ig.RO2.nokey,ig.RO3.nokey,ig.RO4.nokey,ig.RO5.nokey,ig.RO6.nokey,
                               ig.AK1.nokey,ig.AK2.nokey,ig.AK3.nokey,ig.AK4.nokey,ig.AK5.nokey,ig.AK6.nokey)


#########################################
#####  CHARACTERISTICS  of Habitats #####
#########################################
### By habitat VTs

# Centrality-Betweeness
ig.RO_between.nokey <- igraph::betweenness(ig.RO.nokey, weights=E(ig.RO.nokey),normalized=TRUE)
ig.AK_between.nokey <- igraph::betweenness(ig.AK.nokey,weights=E(ig.AK.nokey),normalized=TRUE)

#Connectedness (degree)
ig.RO_degree.nokey <- igraph::degree(ig.RO.nokey, mode="all", normalized=TRUE)
ig.AK_degree.nokey <- igraph::degree(ig.AK.nokey, mode="all", normalized=TRUE)

# Average Path Length: how long are paths around a network?
average.path.length(ig.RO.nokey)
average.path.length(ig.AK.nokey)


#######################################################
#####  CHARACTERISTICS  of Plot means in Habitats #####
#######################################################

### For each plot, then generate means
# Centrality-Betweeness
ig.RO1_between.nokey <- igraph::betweenness(ig.RO1.nokey, weights=E(ig.RO1.nokey), normalized=TRUE)
ig.RO2_between.nokey <- igraph::betweenness(ig.RO2.nokey, weights=E(ig.RO2.nokey), normalized=TRUE)
ig.RO3_between.nokey <- igraph::betweenness(ig.RO3.nokey, weights=E(ig.RO3.nokey), normalized=TRUE)
ig.RO4_between.nokey <- igraph::betweenness(ig.RO4.nokey, weights=E(ig.RO4.nokey), normalized=TRUE)
ig.RO5_between.nokey <- igraph::betweenness(ig.RO5.nokey, weights=E(ig.RO5.nokey), normalized=TRUE)
ig.RO6_between.nokey <- igraph::betweenness(ig.RO6.nokey, weights=E(ig.RO6.nokey), normalized=TRUE)

ig.AK1_between.nokey <- igraph::betweenness(ig.AK1.nokey, weights=E(ig.AK1.nokey), normalized=TRUE)
ig.AK2_between.nokey <- igraph::betweenness(ig.AK2.nokey, weights=E(ig.AK2.nokey), normalized=TRUE)
ig.AK3_between.nokey <- igraph::betweenness(ig.AK3.nokey, weights=E(ig.AK3.nokey), normalized=TRUE)
ig.AK4_between.nokey <- igraph::betweenness(ig.AK4.nokey, weights=E(ig.AK4.nokey), normalized=TRUE)
ig.AK5_between.nokey <- igraph::betweenness(ig.AK5.nokey, weights=E(ig.AK5.nokey), normalized=TRUE)
ig.AK6_between.nokey <- igraph::betweenness(ig.AK6.nokey, weights=E(ig.AK6.nokey), normalized=TRUE)

haka_centrality.nokey <- c(mean(ig.RO1_between.nokey),mean(ig.RO2_between.nokey),mean(ig.RO3_between.nokey),
                     mean(ig.RO4_between.nokey),mean(ig.RO5_between.nokey),mean(ig.RO6_between.nokey),
                     mean(ig.AK1_between.nokey),mean(ig.AK2_between.nokey),mean(ig.AK3_between.nokey),
                     mean(ig.AK4_between.nokey),mean(ig.AK5_between.nokey),mean(ig.AK6_between.nokey))

haka_centrality.nokey <- as.data.frame(haka_centrality.nokey)


#Connectedness (degree)
ig.RO1_degree.nokey <- igraph::degree(ig.RO1.nokey, mode="all", normalized=TRUE)
ig.RO2_degree.nokey <- igraph::degree(ig.RO2.nokey, mode="all", normalized=TRUE)
ig.RO3_degree.nokey <- igraph::degree(ig.RO3.nokey, mode="all", normalized=TRUE)
ig.RO4_degree.nokey <- igraph::degree(ig.RO4.nokey, mode="all", normalized=TRUE)
ig.RO5_degree.nokey <- igraph::degree(ig.RO5.nokey, mode="all", normalized=TRUE)
ig.RO6_degree.nokey <- igraph::degree(ig.RO6.nokey, mode="all", normalized=TRUE)

ig.AK1_degree.nokey <- igraph::degree(ig.AK1.nokey, mode="all", normalized=TRUE)
ig.AK2_degree.nokey <- igraph::degree(ig.AK2.nokey, mode="all", normalized=TRUE)
ig.AK3_degree.nokey <- igraph::degree(ig.AK3.nokey, mode="all", normalized=TRUE)
ig.AK4_degree.nokey <- igraph::degree(ig.AK4.nokey, mode="all", normalized=TRUE)
ig.AK5_degree.nokey <- igraph::degree(ig.AK5.nokey, mode="all", normalized=TRUE)
ig.AK6_degree.nokey <- igraph::degree(ig.AK6.nokey, mode="all", normalized=TRUE)

haka_connectedness.nokey <- c(mean(ig.RO1_degree.nokey),mean(ig.RO2_degree.nokey),mean(ig.RO3_degree.nokey),
                        mean(ig.RO4_degree.nokey),mean(ig.RO5_degree.nokey),mean(ig.RO6_degree.nokey),
                        mean(ig.AK1_degree.nokey),mean(ig.AK2_degree.nokey),mean(ig.AK3_degree.nokey),
                        mean(ig.AK4_degree.nokey),mean(ig.AK5_degree.nokey),mean(ig.AK6_degree.nokey))

haka_connectedness.nokey <- as.data.frame(haka_connectedness.nokey)


## Path length average.path.length(ig.RO1)
ig.RO1_path.nokey <- igraph::average.path.length(ig.RO1.nokey)
ig.RO2_path.nokey <- igraph::average.path.length(ig.RO2.nokey)
ig.RO3_path.nokey <- igraph::average.path.length(ig.RO3.nokey)
ig.RO4_path.nokey <- igraph::average.path.length(ig.RO4.nokey)
ig.RO5_path.nokey <- igraph::average.path.length(ig.RO5.nokey)
ig.RO6_path.nokey <- igraph::average.path.length(ig.RO6.nokey)

ig.AK1_path.nokey <- igraph::average.path.length(ig.AK1.nokey)
ig.AK2_path.nokey <- igraph::average.path.length(ig.AK2.nokey)
ig.AK3_path.nokey <- igraph::average.path.length(ig.AK3.nokey)
ig.AK4_path.nokey <- igraph::average.path.length(ig.AK4.nokey)
ig.AK5_path.nokey <- igraph::average.path.length(ig.AK5.nokey)
ig.AK6_path.nokey <- igraph::average.path.length(ig.AK6.nokey)

haka_path.nokey <- c(mean(ig.RO1_path.nokey),mean(ig.RO2_path.nokey),mean(ig.RO3_path.nokey),
               mean(ig.RO4_path.nokey),mean(ig.RO5_path.nokey),mean(ig.RO6_path.nokey),
               mean(ig.AK1_path.nokey),mean(ig.AK2_path.nokey),mean(ig.AK3_path.nokey),
               mean(ig.AK4_path.nokey),mean(ig.AK5_path.nokey),mean(ig.AK6_path.nokey))

haka_path.nokey <- as.data.frame(haka_path.nokey)

# Dataframe building
plot<-as.data.frame(c("RO1","RO2","RO3","RO4","RO5","RO6",
                      "AK1","AK2","AK3","AK4","AK5","AK6"))
hab_type <- as.data.frame(rep(c("Remnant Forest","Restored Forest"),each=6))
sample_type<- as.data.frame(rep(c("soil"),each=6,times=2))

# add in MEAN traits for each plot, in each network
fungal_networks.nokey <- cbind(plot,hab_type,sample_type,haka_centrality.nokey, haka_connectedness.nokey, haka_path.nokey)
colnames(fungal_networks.nokey) <- c("Plot","HabitatType","SampleType","Centrality","Connectedness", "Path.length")

fungal_networks.nokey %>%
        group_by(HabitatType) %>%
        summarize(means = mean(Centrality), se= se(Centrality))

fungal_networks.nokey %>%
        group_by(HabitatType) %>%
        summarize(means = mean(Connectedness), se= se(Connectedness))

fungal_networks.nokey %>%
        group_by(HabitatType) %>%
        summarize(means = mean(Path.length), se= se(Path.length))

## Welch t-tests
# Centrality by plots
centrality_habitat.nokey <- t.test(subset(fungal_networks.nokey, HabitatType == "Remnant Forest")$Centrality,
                             subset(fungal_networks.nokey, HabitatType == "Restored Forest")$Centrality,
                             paired=FALSE, var.equal=FALSE)
centrality_habitat.nokey

# Connectedness by plots
connectedness_habitat.nokey <- t.test(subset(fungal_networks.nokey, HabitatType == "Remnant Forest")$Connectedness,
                                subset(fungal_networks.nokey, HabitatType == "Restored Forest")$Connectedness,
                                paired=FALSE, var.equal=FALSE)
connectedness_habitat.nokey


# Path Length by plots
path_habitat.nokey <- t.test(subset(fungal_networks.nokey, HabitatType == "Remnant Forest")$Path.length,
                       subset(fungal_networks.nokey, HabitatType == "Restored Forest")$Path.length,
                       paired=FALSE, var.equal=FALSE)
path_habitat.nokey

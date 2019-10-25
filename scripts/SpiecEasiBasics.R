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

# import data, has OTUs in rows and samples in columns. 
otu <- as.matrix(read.csv("data/haka_soil_ESV_table.csv", header = TRUE,row.names = 1))

# remove all OTUs not in atleast 2 samples
otu<-as.matrix(otu[-which(rowMeans(otu) < 2),])

#add a pseudo count and transpose
otu.pc <- (otu)+1

#total sum scaling (relative abundace)
otu.tss <- t(apply(otu.pc, 1, norm_to_total))

# pair rownames from OTU with the taxonomy from 'taxmat'
tax<-merge(otu.tss, taxmat, by = "row.names", all = TRUE) # merge dataframes
tax.OTU<-na.omit(tax)
tax.names<- tax.OTU %>% dplyr::select(Kingdom, Phylum,  Class, Order,  Family,  Genus, Species, Name)

# -------------------------- Models
otu.est.SpEa.MB <- spiec.easi(otu.tss, method='mb', pulsar.params = list(thresh = 0.1)) #mb method
otu.est.SpEa.GL <- spiec.easi(otu.tss, method='glasso', pulsar.params = list(thresh = 0.1)) #glasso fitting

# --------------------------
## Create igraph objects
hak.mb     <- adj2igraph(getRefit(otu.est.SpEa.MB))
hak.gl     <- adj2igraph(getRefit(otu.est.SpEa.GL))

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(otu.tss, 1))+6
am.coord <- layout.fruchterman.reingold(hak.mb)

par(mfrow=c(1,2))
plot(hak.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(hak.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")

# --------------------------
# We can evaluate the weights on edges networks using the terms from the underlying model. 
# SparCC correlations can be used directly, while SpiecEasi networks need to be massaged a bit. 
# Note that since SPIEC-EASI is based on penalized estimators, the edge weights are not directly comparable to SparCC 

secor  <- cov2cor(getOptCov(otu.est.SpEa.GL))
sebeta <- symBeta(getOptBeta(otu.est.SpEa.MB), mode='maxabs')
elist.gl <- summary(triu(secor*getRefit(otu.est.SpEa.GL), k=1))
elist.mb <- summary(sebeta)

par(mfrow=c(1,2))
hist(elist.mb[,3], col='forestgreen', ylim=c(0,300), main="MB")
hist(elist.gl[,3], col='red', ylim=c(0,300), main="glasso")

# Lets look at the degree statistics from the networks inferred by each method.
dd.gl     <- degree.distribution(hak.gl)
dd.mb     <- degree.distribution(hak.mb)

par(mfrow=c(1,1))
plot(0:(length(dd.gl)-1), dd.gl, col="red" , type='b', ylab="Frequency", xlab="Degree", ylim=c(0,0.2), main="Degree Distributions")
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso"),
       col=c("forestgreen", "red"), pch=1, lty=1)

#------------- plots 
pal <- brewer.pal(10, "BrBG")
Group<- tax.names$Family
vertex.col <- pal[Group]
par(mfrow=c(1,1), mar=c(0,1,1,1))

plot(hak.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso",
     edge.color="gray70", vertex.color=vertex.col, col=1:10)
legend('topright',legend=levels(Group), pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.9, pt.cex=1.3, y.intersp = 0.4, inset=c(0.4, 0))





########### Divide up: 
# first, the otu table is 'otu'
# second, the taxonomy table is an exported phyloseq object 'ESV_dataframe'


### RO (Remnant Forest) sites

## subset otu table
# rows as OTUs, columns as samples
ROdat<-otu[,grep("^RO", colnames(otu))] # data with only RO

# subset RO matrix, ESVs in >2 samples, add pseudo value (1), transpose, relative abundance
RO.otu<-as.matrix(ROdat[-which(rowMeans(ROdat) < 2),])
RO.otu.pc <- (RO.otu)+1

# this transposes data 2x, keeps rows as OTUs
RO.otu.tss <- t(apply(RO.otu.pc, 1, norm_to_total)) 

### BUT the columns need to be OTUs!!!

# pair rownames from OTU with the taxonomy from 'taxmat'
RO.tax<-merge(RO.otu.tss, taxmat, by = "row.names", all = TRUE) # merge dataframes
RO.tax<-droplevels(RO.tax)
RO.tax<-na.omit(RO.tax)
RO.tax.names<- RO.tax %>% dplyr::select(Kingdom, Phylum,  Class, Order,  Family,  Genus, Species, Name)

## now transpose the otus in analysis below so that OTUs are columns, samples rows
RO.otu.SE<-t(RO.otu.tss)

# Models with glassofitting
RO.otu.est.SpEa.GL <- spiec.easi(RO.otu.SE, method='glasso', pulsar.params = list(thresh = 0.1)) 
RO.otu.est.SpEa.MB <- spiec.easi(RO.otu.SE, method='mb', pulsar.params = list(thresh = 0.1)) 

# Create igraph objects
igRO.gl <- adj2igraph(getRefit(RO.otu.est.SpEa.GL))
igRO.mb <- adj2igraph(getRefit(RO.otu.est.SpEa.MB))

## set size of vertex proportional to clr-mean
RO.vsize    <- rowMeans(clr(RO.otu.SE, 1))+6
RO.am.coord <- layout.fruchterman.reingold(igRO.mb)

# Lets look at the degree statistics from the networks inferred by each method.
dd.ROgl     <- degree.distribution(igRO.gl)
dd.ROmb     <- degree.distribution(igRO.mb)



########### Divide up: AK sites 
# AK

## subset otu table
AKdat<-otu[,grep("^AK", colnames(otu))] # data with only RO

# subset AK matrix, ESVs in >2 samples, add pseudo value (1), transpose, relative abundance
AK.otu<-as.matrix(AKdat[-which(rowMeans(AKdat) < 2),])
AK.otu.pc <- (AK.otu)+1
# this transposes data 2x, keeps rows as OTUs
AK.otu.tss <- t(apply(AK.otu.pc, 1, norm_to_total)) 

# pair rownames from OTU with the taxonomy from 'taxmat'
AK.tax<-merge(AK.otu.tss, taxmat, by = "row.names", all = TRUE) # merge dataframes
AK.tax<-droplevels(AK.tax)
AK.tax<-na.omit(AK.tax)
AK.tax.names<- AK.tax %>% dplyr::select(Kingdom, Phylum,  Class, Order,  Family,  Genus, Species, Name)

## now transpose the otus in analysis below so that OTUs are columns, samples rows
AK.otu.SE<-t(AK.otu.tss)

# -------------------------- Models
AK.otu.est.SpEa.GL <- spiec.easi(AK.otu.SE, method='glasso', pulsar.params = list(thresh = 0.1))
AK.otu.est.SpEa.MB <- spiec.easi(AK.otu.SE, method='mb', pulsar.params = list(thresh = 0.1)) 

# Create igraph objects
igAK.gl <- adj2igraph(getRefit(AK.otu.est.SpEa.GL))
igAK.mb <- adj2igraph(getRefit(AK.otu.est.SpEa.MB))

## set size of vertex proportional to clr-mean
AK.vsize    <- rowMeans(clr(AK.otu.SE, 1))+6
AK.am.coord <- layout.fruchterman.reingold(igAK.mb)

# Lets look at the degree statistics from the networks inferred by each method.
dd.AKgl     <- degree.distribution(igAK.gl)
dd.AKmb     <- degree.distribution(igAK.mb)


######## combined plots

## networks

pal <- brewer.pal(10, "BrBG")
GroupRO <- t(RO.tax.names$Family)
vertex.col <- pal[GroupRO]
par(mfrow=c(1,1), mar=c(0,1,1,1))

plot(igRO.gl, layout=RO.am.coord, vertex.size=RO.vsize, vertex.label=NA, main="Remnant Forest", 
     vertex.color=vertex.col, col=1:10) #"#336B87"
legend('topright',legend=levels(GroupRO), pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.8, pt.cex=1.2, y.intersp = 0.6, inset=c(0.6, 0.025))
dev.copy(pdf, "figures/ROnetwork.pdf", height=6, width=7)
dev.off() 


GroupAK <- t(AK.tax.names$Family)
vertex.col <- pal[GroupAK]
plot(igAK.gl, layout=AK.am.coord, vertex.size=AK.vsize, vertex.label=NA, main="Restored Forest", 
     vertex.color=vertex.col, col=1:10) #"#336B87"
legend('topright',legend=levels(GroupAK), pch=16, col=pal, border=NA, box.lty=0, bg="transparent", 
       cex=0.8, pt.cex=1.2, y.intersp = 0.6, inset=c(0.6, 0.05))
dev.copy(pdf, "figures/AKnetwork.pdf", height=6, width=7)
dev.off() 


### degree distance
par(mfrow=c(1,2), mar=c(3,3,2,2))

# RO degree distribution
plot(0:(length(dd.ROgl)-1), dd.ROgl, col="red" , type='b', ylab="Frequency", xlab="Degree", 
     ylim=c(0,0.2), xlim=c(0,20), main="RO Degree Distributions")
points(0:(length(dd.ROmb)-1), dd.ROmb, col="forestgreen", type='b')

# AK degree distribution
plot(0:(length(dd.AKgl)-1), dd.AKgl, col="red" , type='b', ylab="Frequency", xlab="Degree",
     ylim=c(0,0.2), xlim=c(0,20), main="AK Degree Distributions")
points(0:(length(dd.AKmb)-1), dd.AKmb, col="forestgreen", type='b')
legend("topleft", c("MB", "glasso"), box.lty=0, bg="transparent", y.intersp = 0.4,
       col=c("forestgreen", "red"), pch=1, lty=1)

dev.copy(pdf, "figures/network.dd.pdf", height=5, width=7)
dev.off() 

##########

# Lets look at the degree statistics from the networks inferred by each method.
d.AKgl<-as.data.frame(degree(igAK.gl, mode = "total")); colnames(d.AKgl)<-"degree"
d.ROgl <- as.data.frame(degree(igRO.mb, mode = "total")); colnames(d.ROgl)<-"degree"

par(mfrow=c(1,1))
plot(0:(length(d.AKgl$degree)-1), d.AKgl$degree, pch=21, col="red", cex=0.5)
points(0:(length(d.ROgl$degree)-1), d.ROgl$degree, col="forestgreen", pch=21, cex=0.5)
legend("topleft", c("RO", "AK"), box.lty=0, bg="transparent", y.intersp = 0.4,
       col=c("red", "forestgreen"), pch=1, lty=1)







### TO TEST
# NODE DEGREE, BETWEENNESS-CENTRALITY TEST
# directed vs. undirected degrees check




############## Using phyloseq

### Subset phyloseq object into different hab types and then by roots and soil
RO.physeq = subset_samples(haka_soil_physeq,HabitatType == "Remnant Forest")
AK.physeq = subset_samples(haka_soil_physeq,HabitatType == "Restored Forest")

### Construct networks using Spiec Easi (save in case R crashes as an RDS file)
se.RO<- spiec.easi(RO.physeq, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.RO,"output/spiec_easi_files/se.RO.RData")
se.RO<- readRDS("output/spiec_easi_files/se.RO.RData")

se.AK<- spiec.easi(AK.physeq, method='mb', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(rep.num=50))
saveRDS(se.AK,"output/spiec_easi_files/se.AK.RData")
se.AK <- readRDS("output/spiec_easi_files/se.AK.RData")

### Convert to igraph objects for plotting and computing graph statistics.
#Bring in edge weights
ig.RO<-adj2igraph(se.RO$refit$stars, vertex.attr=list(name=taxa_names(ro.physeq)))
ig.AK<-adj2igraph(se.AK$refit$stars, vertex.attr=list(name=taxa_names(ak.physeq)))

V(ig.RO)$label<-V(ig.RO)$name
V(ig.RO)$frame.color = "black"

V(ig.AK)$label<-V(ig.RO)$name
V(ig.AK)$frame.color = "black"

# now both are ready for igraph and plotting


#########  RO network in igrpah
#Calculate number of samples each AM fungal species (node) is detected, and add as vertex attribute
RO_dat =fast_melt(RO.physeq)
RO_no.samples = RO_dat[, list(no.samples=sum(count>0)),by=TaxaID]
RO_no.samples = as.data.frame(RO_no.samples)
V(ig.RO)$no.samples = as.numeric(RO_no.samples$no.samples[match(V(ig.RO)$name, RO_no.samples$TaxaID)])
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

pal<-brewer.pal(10,"Paired")
GroupRO <- V(ig.RO)$Family
vertex.col <- pal[V(ig.RO)$Family]

par(mfrow=c(1,1), mar=c(0,1,1,1))

# if want to set colors manually with names....
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
par(mfrow=c(1,1))
par(mar=c(0,2,2,0))
set.seed(606)
plot(ig.RO, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.RO)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA)
legend('topright',legend=levels(as.factor(V(ig.RO)$Family)), pch=16, col=levels(as.factor(V(ig.RO)$color)), border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.6, inset=c(0.6, 0.025))
dev.copy(pdf, "figures/network.RO.mb.pdf", height=5, width=7)
dev.off() 
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

#Colour by Family. First need to assign colours to match previously made graphs
net_tax=read.csv("data/haka_soil_taxonomy.csv", header = TRUE, row.names = 1)
V(ig.AK)$Family=as.character(net_tax$Family[match(V(ig.AK)$name, net_tax$OTU)])


par(mfrow=c(1,1), mar=c(0,1,1,1))

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
par(mfrow=c(1,1))
par(mar=c(0,2,2,0))
set.seed(606)
plot(ig.AK, vertex.size=7, vertex.shape="circle", vertex.color=V(ig.AK)$color,
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA)
legend('topright',legend=levels(as.factor(V(ig.AK)$Family)), pch=16, col=levels(as.factor(V(ig.AK)$color)), border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.6, inset=c(0.6, 0.025))
dev.copy(pdf, "figures/network.AK.mb.pdf", height=5, width=7)
dev.off() 
######### 


#Manuscript plot
par(mar=c(0.75,0.75,0.75,0.75))
par(mfrow=c(1,2))
set.seed(606)
plot(ig.RO, vertex.size=7, vertex.shape="circle",vertex.color=V(ig.ro.r)$color, 
     edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Remnant Forest")
set.seed(606)
plot(ig.AK, vertex.shape="circle",vertex.color=V(ig.ak.r)$color, 
     vertex.size= 7, edge.width=3, edge.color="dark grey",layout=layout_with_fr,
     vertex.label=NA, main="Restored Forest")
legend('topright',legend=levels(as.factor(V(ig.AK)$Family)), pch=16, col=levels(as.factor(V(ig.AK)$color)), border=NA, box.lty=0, bg="transparent", 
       cex=0.7, pt.cex=1.1, y.intersp = 0.4, inset=c(0.1, 0.05))

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
dev.copy(pdf, "figures/dd.mb.habitat.pdf", height=5, width=7)
dev.off() 
######### 

title = "plot_bar"
plot_bar(AK.physeq, "Family", "Abundance")

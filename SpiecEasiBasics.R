library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)
library(car)
library(huge)
library(MASS)
library(Matrix)

# import data, has OTUs in rows and samples in columns. 
otu <- as.matrix(read.csv("data/haka_soil_ESV_table.csv", header = TRUE,row.names = 1))

# remove all OTUs not in atleast 2 samples
otu<-as.matrix(otu[-which(rowMeans(otu) < 2),])

#add a pseudo count and transpose
otu.pc <- (otu)+1

#total sum scaling (relative abundace)
otu.tss <- t(apply(otu.pc, 1, norm_to_total))
             
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
plot(hak.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso",
     edge.color="gray70", vertex.color="coral") # edge = lines, vertex = symbols






########### Divide up: RO sites
# RO

## subset otu table
ROdat<-otu[,grep("^RO", colnames(otu))] # data with only RO

# subset RO matrix, ESVs in >2 samples, add pseudo value (1), transpose, relative abundance
RO.otu<-as.matrix(ROdat[-which(rowMeans(ROdat) < 2),])
RO.otu.pc <- (RO.otu)+1
RO.otu.tss <- t(apply(RO.otu.pc, 1, norm_to_total))

# Models with glassofitting
RO.otu.est.SpEa.GL <- spiec.easi(RO.otu.tss, method='glasso', pulsar.params = list(thresh = 0.1)) 
RO.otu.est.SpEa.MB <- spiec.easi(RO.otu.tss, method='mb', pulsar.params = list(thresh = 0.1)) 

# Create igraph objects
igRO.gl <- adj2igraph(getRefit(RO.otu.est.SpEa.GL))
igRO.mb <- adj2igraph(getRefit(RO.otu.est.SpEa.MB))

## set size of vertex proportional to clr-mean
RO.vsize    <- rowMeans(clr(RO.otu.tss, 1))+6
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
AK.otu.tss <- t(apply(AK.otu.pc, 1, norm_to_total))

# -------------------------- Models
AK.otu.est.SpEa.GL <- spiec.easi(AK.otu.tss, method='glasso', pulsar.params = list(thresh = 0.1))
AK.otu.est.SpEa.MB <- spiec.easi(RO.otu.tss, method='mb', pulsar.params = list(thresh = 0.1)) 

# Create igraph objects
igAK.gl <- adj2igraph(getRefit(AK.otu.est.SpEa.GL))
igAK.mb <- adj2igraph(getRefit(AK.otu.est.SpEa.MB))

## set size of vertex proportional to clr-mean
AK.vsize    <- rowMeans(clr(AK.otu.tss, 1))+6
AK.am.coord <- layout.fruchterman.reingold(igAK.mb)

# Lets look at the degree statistics from the networks inferred by each method.
dd.AKgl     <- degree.distribution(igAK.gl)
dd.AKmb     <- degree.distribution(igAK.mb)

######## combined plots

## networks
par(mfrow=c(1,2), mar=c(0,1,1,1))
plot(igRO.gl, layout=RO.am.coord, vertex.size=RO.vsize, vertex.label=NA, main="Remnant Forest", 
     vertex.color="#88A550") 
plot(igAK.gl, layout=AK.am.coord, vertex.size=AK.vsize, vertex.label=NA, main="Restored Forest", 
     vertex.color="#336B87")
dev.copy(pdf, "figures/network.pdf", height=5, width=7)
dev.off() 


### degree distance
par(mfrow=c(1,2))
# RO degree distribution
plot(0:(length(dd.ROgl)-1), dd.ROgl, col="red" , type='b', ylab="Frequency", xlab="Degree", 
     ylim=c(0,0.2), xlim=c(0,25), main="RO Degree Distributions")
points(0:(length(dd.ROmb)-1), dd.ROmb, col="forestgreen", type='b')

# AK degree distribution
plot(0:(length(dd.AKgl)-1), dd.AKgl, col="red" , type='b', ylab="Frequency", xlab="Degree",
     ylim=c(0,0.2), xlim=c(0,25), main="AK Degree Distributions")
points(0:(length(dd.AKmb)-1), dd.AKmb, col="forestgreen", type='b')
legend("topleft", c("MB", "glasso"), box.lty=0, bg="transparent", y.intersp = 0.4,
       col=c("forestgreen", "red"), pch=1, lty=1)

dev.copy(pdf, "figures/network.dd.pdf", height=5, width=7)
dev.off() 


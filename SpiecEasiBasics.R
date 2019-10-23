library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)
library(car)
library(huge)
library(MASS)
librarY(Matrix)

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




########### Divide up and run RO and AK sites
# RO

## subset otu table
ROdat<-otu[,grep("^RO", colnames(otu))] # data with only RO

# subset RO matrix, ESVs in >2 samples, add pseudo value (1), transpose, relative abundance
RO.otu<-as.matrix(ROdat[-which(rowMeans(ROdat) < 2),])
RO.otu.pc <- (RO.otu)+1
RO.otu.tss <- t(apply(RO.otu.pc, 1, norm_to_total))

# -------------------------- Models
RO.otu.est.SpEa.GL <- spiec.easi(RO.otu.tss, method='glasso', pulsar.params = list(thresh = 0.1)) #glasso fitting

library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(igraph)
library(car)
library(huge)
library(MASS)

# import data
otu <- as.matrix(read.csv("data/haka_soil_ESV_table.csv", header = TRUE,row.names = 1))


#'otu' has OTUs in rows and samples in columns. it contains raw counts and all OTUs are present in 2+ samples.
#add a pseudo count
otu.pc <- t(otu)+1
#total sum scaling (relative abundace)
otu.tss <- t(apply(otu.pc, 1, norm_to_total))
otu.est.SpEa <- spiec.easi(otu.tss, method='glasso', pulsar.params = list(thresh = 0.1))




#plot with nodes as blue circles and all edges grey
plot.igraph(adj2igraph(otu.sel$refit), vertex.shape="circle", vertex.label=NA, vertex.size=4, edge.arrow.mode=0)

#get the optimal parameters in order to color the edges (there's probably a better way to do this)
otu.cov <- otu.sel$opt.cov
otu.cov[which(otu.sel$refit==0)] <- 0
otu.adj <- otu.sel$refit
otu.adj[which(otu.cov<0)] <- -1

edgecol <- recode(otu.cov[upper.tri(otu.cov)], "-1:-0.00001='red'; 0.00001:1='green'; 0='black'; else='grey'")
table(edgecol)
#81% (383/470 edges) positive

#plot with nodes as circles colored by class (vector clscol is a vector of color names for each OTU), and edges colored by positive (green) or negative (red)
plot.igraph(adj2igraph(otu.sel2$refit), vertex.shape="circle", vertex.color=clscol, vertex.label=NA, vertex.size=4, edge.arrow.mode=0, 
            edge.color=edgecol[-which(edgecol=="black")])

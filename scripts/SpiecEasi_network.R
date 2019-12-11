############################################################################################################
### Composite co-occurrence networks (habitat type x sample type) for visual and identifying keystone sp ###
############################################################################################################

###### ###### ###### ###### 
########## packages
###### ###### ###### ###### 

if (!require("pacman")) install.packages("pacman"); library(pacman) # for rapid install if not in library
if (!require("BiocManager")) install.packages("BiocManager") # for buoinformatic packages
if (!require("devtools")) install.packages("devtools") # for developement tools
install_github("zdk123/SpiecEasi")

devtools::install_github('oswaldosantos/ggsn')

pacman::p_load("ade4", "multtest","car", "phyloseq","rhdf5","ggplot2","colorspace","stringi", "geosphere", 
               "ggplot2", "ggmap", "dplyr", "gridExtra", "geosphere", "sf", "raster", "spData", "cowplot", "huge",
               "tmap", "leaflet", "mapview", "shiny", "fossil","igraph","SpiecEasi", "scales", "RgoogleMaps", "devtools", "ggsn", "vegan", "multcomp",
               "dplyr", "grid", "scales", "gridExtra", "emmeans", "Matrix", "MASS", "multcompView", "ggpubr", "Rmisc", "purrr",
               "RVAideMemoire", "RColorBrewer", "vegan")




### Fast Melt

# Define some functions for quick taxa summary

fast_melt = function(physeq,
                     includeSampleVars = character(),
                     omitZero = FALSE){
    require("phyloseq")
    require("data.table")
    # supports "naked" otu_table as `physeq` input.
    otutab = as(otu_table(physeq), "matrix")
    if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
    otudt = data.table(otutab, keep.rownames = TRUE)
    setnames(otudt, "rn", "TaxaID")
    # Enforce character TaxaID key
    otudt[, TaxaIDchar := as.character(TaxaID)]
    otudt[, TaxaID := NULL]
    setnames(otudt, "TaxaIDchar", "TaxaID")
    # Melt count table
    mdt = melt.data.table(otudt, 
                          id.vars = "TaxaID",
                          variable.name = "SampleID",
                          value.name = "count")
    if(omitZero){
        # Omit zeroes and negative numbers
        mdt <- mdt[count > 0]
    }
    # Omit NAs
    mdt <- mdt[!is.na(count)]
    # Calculate relative abundance
    mdt[, RelativeAbundance := count / sum(count), by = SampleID]
    if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
        # If there is a tax_table, join with it. Otherwise, skip this join.
        taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
        setnames(taxdt, "rn", "TaxaID")
        # Enforce character TaxaID key
        taxdt[, TaxaIDchar := as.character(TaxaID)]
        taxdt[, TaxaID := NULL]
        setnames(taxdt, "TaxaIDchar", "TaxaID")
        # Join with tax table
        setkey(taxdt, "TaxaID")
        setkey(mdt, "TaxaID")
        mdt <- taxdt[mdt]
    }
    # includeSampleVars = c("DaysSinceExperimentStart", "SampleType")
    # includeSampleVars = character()
    # includeSampleVars = c()
    # includeSampleVars = c("aksjdflkas") 
    wh.svars = which(sample_variables(physeq) %in% includeSampleVars)
    if( length(wh.svars) > 0 ){
        # Only attempt to include sample variables if there is at least one present in object
        sdf = as(sample_data(physeq), "data.frame")[, wh.svars, drop = FALSE]
        sdt = data.table(sdf, keep.rownames = TRUE)
        setnames(sdt, "rn", "SampleID")
        # Join with long table
        setkey(sdt, "SampleID")
        setkey(mdt, "SampleID")
        mdt <- sdt[mdt]
    }
    setkey(mdt, "TaxaID")
    return(mdt)
}

summarize_taxa = function(physeq, Rank, GroupBy = NULL){
    require("phyloseq")
    require("data.table")
    Rank <- Rank[1]
    if(!Rank %in% rank_names(physeq)){
        message("The argument to `Rank` was:\n", Rank,
                "\nBut it was not found among taxonomic ranks:\n",
                paste0(rank_names(physeq), collapse = ", "), "\n",
                "Please check the list shown above and try again.")
    }
    if(!is.null(GroupBy)){
        GroupBy <- GroupBy[1]
        if(!GroupBy %in% sample_variables(physeq)){
            message("The argument to `GroupBy` was:\n", GroupBy,
                    "\nBut it was not found among sample variables:\n",
                    paste0(sample_variables(physeq), collapse = ", "), "\n",
                    "Please check the list shown above and try again.")
        }
    }
    # Start with fast melt
    mdt = fast_melt(physeq)
    if(!is.null(GroupBy)){
        # Add the variable indicated in `GroupBy`, if provided.
        sdt = data.table(SampleID = sample_names(physeq),
                         var1 = get_variable(physeq, GroupBy))
        setnames(sdt, "var1", GroupBy)
        # Join
        setkey(sdt, SampleID)
        setkey(mdt, SampleID)
        mdt <- sdt[mdt]
    }
    # Summarize
    if(!is.null(GroupBy)){
        summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                               sdRA = sd(RelativeAbundance),
                               minRA = min(RelativeAbundance),
                               maxRA = max(RelativeAbundance)),
                        by = c(Rank, GroupBy)]
    } else {
        Nsamples = nsamples(physeq)
        # No GroupBy argument, can be more precise with the mean, sd, etc.
        summarydt = mdt[, list(meanRA = sum(RelativeAbundance) / Nsamples,
                               sdRA = sd(c(RelativeAbundance, numeric(Nsamples - .N))),
                               minRA = ifelse(test = .N < Nsamples,
                                              yes = 0L, 
                                              no = min(RelativeAbundance)),
                               maxRA = max(RelativeAbundance)),
                        by = c(Rank)]
    }
    return(summarydt)
}

plot_taxa_summary = function(physeq, Rank, GroupBy = NULL){
    require("phyloseq")
    require("data.table")
    require("ggplot2")
    # Get taxa summary table 
    dt1 = summarize_taxa(physeq, Rank = Rank, GroupBy = GroupBy)
    # Set factor appropriately for plotting
    RankCol = which(colnames(dt1) == Rank)
    setorder(dt1, -meanRA)
    dt1[, RankFac := factor(dt1[[Rank]], 
                            levels = rev(dt1[[Rank]]))]
    dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
    dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
    # Set zeroes to one-tenth the smallest value
    ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
    ebarMinFloor <- ebarMinFloor / 10
    dt1[(ebarMin == 0), ebarMin := ebarMinFloor]
    
    pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
        scale_x_log10() +
        xlab("Mean Relative Abundance") +
        ylab(Rank) +
        theme_bw()
    if(!is.null(GroupBy)){
        # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
        pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                    size = 5)
    } else {
        # Don't include error bars for faceted version
        pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                            xmin = ebarMin))
    }
    return(pRank)
}
################



############## Using phyloseq
#####################
### SEQUENCE DATA ###
#####################

###################################
####Import files and pre-process###
###################################

#Import files generated in QIIME
otu <- as.matrix(read.csv("data/haka_soil_ESV_table.csv", header = TRUE,row.names = 1))
haka_otu <- t(otu)
taxmat <- as.matrix(read.csv("data/haka_soil_taxonomy.csv", header = TRUE,row.names = 1))
sample <- haka_meta
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
                                    Na.cation= sample$Na.cation, Longitude= sample$Longitude, 
                                    Latitude= sample$Latitude, stringsAsFactors = FALSE))

sample$Hab.by.hot<-interaction(sample$HabitatType, sample$Host)

row.names(sampledata) <- row.names(sample)


#Change each file to the phyloseq format
OTU = otu_table(haka_otu, taxa_are_rows = FALSE)
physeq = phyloseq(OTU)
TAX = tax_table(taxmat)
haka_soil_physeq = merge_phyloseq(OTU, sampledata, TAX)

#Create new physeq object working at only the species level taxonomy. Collapses all ESVs identified as the same species
haka_VT_soil_physeq <- tax_glom(haka_soil_physeq,"Species")



###################
################### Subset data for network analysis
###################

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
V(ig.RO)$Genus=as.character(net_tax$Genus[match(V(ig.RO)$name, net_tax$OTU)])

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
V(ig.AK)$Genus=as.character(net_tax$Genus[match(V(ig.AK)$name, net_tax$OTU)]) # add Genus to datafram

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


##############################
#####   KEYSTONE SPECIES   ###
##############################

# Remnant Forest
RO.between <- as.data.frame(V(ig.RO)$betweenness,normalized=TRUE)
RO.degree <- as.data.frame(V(ig.RO)$degree,mode="all",normalized=TRUE)
RO.rel.abund <- as.data.frame(V(ig.RO)$rel.abund)
RO.fam <- as.data.frame(V(ig.RO)$Family)
RO.gen <- as.data.frame(V(ig.RO)$Genus)
RO.species <- as.data.frame(net_tax$Species[match(V(ig.RO)$name,net_tax$OTU)])
RO.no.samples <- as.data.frame(V(ig.RO)$no.samples)
sample_data(RO.physeq)
RO.total.samples <- as.data.frame(rep(212),times=130) # 212 sample, 130 taxa
RO.keystone <- cbind(RO.between,RO.degree,RO.rel.abund,RO.fam, RO.gen,RO.species,RO.no.samples,RO.total.samples)
rownames(RO.keystone) <- V(ig.RO)$name
colnames(RO.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family", "Genus", "Species","No.Samples","TotalSamples")
RO.keystone$No.Samples <- as.numeric(RO.keystone$No.Samples)
RO.keystone$perc.samples <- RO.keystone$No.Samples / RO.keystone$TotalSamples
RO.keystone$Prevalence <- RO.keystone$perc.samples * RO.keystone$RelativeAbundance


# Restored Forest
AK.between <- as.data.frame(V(ig.AK)$betweenness,normalized=TRUE)
AK.degree <- as.data.frame(V(ig.AK)$degree,mode="all",normalized=TRUE)
AK.rel.abund <- as.data.frame(V(ig.AK)$rel.abund)
AK.fam <- as.data.frame(V(ig.AK)$Family)
AK.gen <- as.data.frame(V(ig.AK)$Genus)
AK.species <- as.data.frame(net_tax$Species[match(V(ig.AK)$name,net_tax$OTU)])
AK.no.samples <- as.data.frame(V(ig.AK)$no.samples)
sample_data(AK.physeq)
AK.total.samples <- as.data.frame(rep(264),times=130) # 264 samples, 130 taxa
AK.keystone <- cbind(AK.between,AK.degree,AK.rel.abund,AK.fam,AK.gen,AK.species,AK.no.samples,AK.total.samples)
rownames(AK.keystone) <- V(ig.AK)$name
colnames(AK.keystone) <- c("Betweenness","Degree","RelativeAbundance","Family","Genus","Species","No.Samples","TotalSamples")
AK.keystone$No.Samples <- as.numeric(AK.keystone$No.Samples)
AK.keystone$perc.samples <- AK.keystone$No.Samples / AK.keystone$TotalSamples
AK.keystone$Prevalence <- AK.keystone$perc.samples * AK.keystone$RelativeAbundance


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
    scale_colour_manual(values=pal, limits=GroupAK) +
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
plot_grid(KS.plots, KS.legend, rel_widths = c(3, 1)) # legend column 1/3 size as first object
dev.copy()
ggsave("figures/keystone.fig.pdf", width = 7, height = 4)

### keystone species of interest: 
### RO: Claroideoglomeraceae spp. VTX00225 (8 node degrees, 193 Between-Central)
### AK: Acaulosporaceae spp. VTX00272, (11 node degrees, 547 Between-Central)
haka_VT_soil_physeq # original full phyloseq object with VT species

# relative abundance and set low counts to zero

ESV_rel_abund <- transform_sample_counts(haka_VT_soil_physeq,function(x)x/sum(x))
ESV_dataframe<-psmelt(ESV_rel_abund)
ESV_dataframe$Abundance<- ifelse(ESV_dataframe$Abundance < 1e-5, 0, ESV_dataframe$Abundance) 

# make mean relative abundance for each VT at levels of interest
Abund<-aggregate(Abundance~HabitatType+Plot+TreeID+Longitude+Latitude+Family+Genus+Species, data=ESV_dataframe, FUN=mean)

# RO keystone mean relative abundance
Abund.Claro<-Abund[(Abund$Species=="VTX00225"),] # RO
sum(Abund.Claro$Abundance > 0) # 6 trees

# AK keystone mean relative abundance
Abund.Acaul<-Abund[(Abund$Species=="VTX00272"),] # AK 
sum(Abund.Acaul$Abundance > 0) # 72 trees

AK.key.df.plot<-aggregate(Abundance~HabitatType+Plot+Species, data=Abund.Acaul, FUN=mean)
AK.key.dfmean<-aggregate(Abundance~HabitatType+Species, data=Abund.Acaul, FUN=mean)
AK.key.dfsd<-aggregate(Abundance~HabitatType+Species, data=Abund.Acaul, FUN=sd)

AK.key.test <- t.test(subset(AK.key.df.plot, HabitatType == "Remnant Forest")$Abundance,
                             subset(AK.key.df, HabitatType == "Restored Forest")$Abundance,
                             paired=FALSE, var.equal=FALSE)


#### Ubiquity and Abundance plot
 # RO phyloseq object = 212 samples
 # AK phyloseq object = 264 samples

# mean abundance across all samples, separated by habitat type
# in each sample (tree ID) is the VT present
Abund$Pres<-ifelse(Abund$Abundance > 0, 1, 0)
Abund$tot.samp<-ifelse(Abund$HabitatType=="Remnant Forest", 212, 264)
Abund$prop.sampl<-Abund$Pres/Abund$tot.samp

Ubiq.tot<-aggregate(prop.sampl~HabitatType+Family+Genus+Species, data=Abund, FUN=sum) # sum of sample propor.
colnames(Ubiq.tot)[5] <- "Ubiq"

Abund.reorder<-aggregate(Abundance~HabitatType+Family+Genus+Species, data=ESV_dataframe, FUN=mean)
colnames(Abund.reorder)[5] <- "Abund"

Abund.SD<-aggregate(Abundance~HabitatType+Family+Genus+Species, data=ESV_dataframe, FUN=sd)
colnames(Abund.SD)[5] <- "Abund.SD"

Abund.Ubiq<-cbind(Abund.reorder, Abund.SD[5], Ubiq.tot[5])

AbUb.RO<-Abund.Ubiq[(Abund.Ubiq$HabitatType=="Remnant Forest"),]
AbUb.RO<-AbUb.RO[(AbUb.RO$Abund > 0),]

AbUb.AK<-Abund.Ubiq[(Abund.Ubiq$HabitatType=="Restored Forest"),]
AbUb.AK<-AbUb.AK[(AbUb.AK$Abund > 0),]

RO.AbUb.plot <- ggplot(AbUb.RO,aes(x=Abund,y=Ubiq)) +
    geom_point(aes(size=Abund.SD,colour=Family),position="jitter") +
    ylab("Ubiquity (proporiton)") + 
    xlab("Relative abundance") +
    scale_x_continuous(limits=c(0,0.4), breaks=seq(0,1, by=0.2)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1, by=0.2))+
    theme(text=element_text(colour="black",size=10)) + 
    theme(axis.text.x=element_text(hjust=1,colour="black",size=10)) +
    theme(axis.text.y=element_text(colour="black",size=10)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background=element_blank()) +
    scale_colour_manual(values=pal, limits=GroupRO) +
    guides(colour=guide_legend(override.aes=list(size=3))) +
    theme(legend.key=element_blank())+ theme(legend.text=element_text(size=8)) +
    theme(legend.key.size = unit(.4, "cm"))+
    ggtitle("A) Remnant Forest")
RO.AbUb.plot

AK.AbUb.plot <- ggplot(AbUb.AK,aes(x=Abund,y=Ubiq)) +
    geom_point(aes(size=Abund.SD,colour=Family),position="jitter") +
    ylab("") + 
    xlab("Relative abundance") +
    scale_x_continuous(limits=c(0,0.4), breaks=seq(0,1, by=0.2)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1, by=0.2))+
    theme(text=element_text(colour="black",size=10)) + 
    theme(axis.text.x=element_text(hjust=1,colour="black",size=10)) +
    theme(axis.text.y=element_text(colour="black",size=10)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background=element_blank()) +
    scale_colour_manual(values=pal, limits=GroupAK) +
    guides(colour=guide_legend(override.aes=list(size=3))) +
    theme(legend.key=element_blank())+ theme(legend.text=element_text(size=8)) +
    theme(legend.key.size = unit(.4, "cm"))+
    ggtitle("B) Restored Forest")
AK.AbUb.plot


# extract legend from plot
RO.AbUb.leg <- get_legend(
    # create some space to the left of the legend
    RO.AbUb.plot + theme(legend.box.margin = margin(0, 0, 0, 12)))

AbUb.plots<- plot_grid(RO.AbUb.plot + theme(legend.position = "none"), 
                         AK.AbUb.plot + theme(legend.position = "none"), 
                     ncol=2,nrow=1)
plot_grid(AbUb.plots, RO.AbUb.leg, rel_widths = c(3, 1)) # legend column 1/2 size as first object
dev.copy()
ggsave("figures/AbUb.fig.pdf", width = 8, height = 4)

##############################################################################################################
###### What the hell, let's Krig baby
##############################################################################################################
# load packages
library(sp)
library(gstat)

suppressPackageStartupMessages({
    library(dplyr) # for "glimpse"
    library(ggplot2)
    library(scales) # for "comma"
    library(magrittr)
})

# subset AK abundance data
AK.abund<-Abund.Acaul

# bubble plot of abundance
Keystone.bubble<-ggplot(AK.abund, aes(x = Longitude, y = Latitude)) + 
    geom_point(aes(color = HabitatType, size = Abundance), alpha = 0.5)+
    coord_equal() +
    scale_color_manual(values = c("#88A550","#336B87")) + theme_bw()
Keystone.bubble
ggsave("figures/Keystone.bubble.png",width= 8,height=8, plot=Keystone.bubble)



# spatial coordinate object
coordinates(AK.abund)<- ~Longitude +Latitude # coordinates for samples of interest
AK.abund@bbox # extend of binding box

Longitude<-c(-155.298, -155.33)
Latitude<-c(19.815, 19.835)
xy<-cbind(Longitude,Latitude)
S<-SpatialPoints(xy)
bbox(S)

AK.abund@bbox<-bbox(S) # expand binding box 


col.scheme.N <- colorRampPalette(c("white", "chartreuse3",'red'))(20)
Grid.AK.KEY <- spsample(AK.abund, type='regular', n=1e4)
gridded(Grid.AK.KEY) <- TRUE

krig.Key <- krige(Abundance ~ 1, AK.abund, Grid.AK.KEY) # ordinary kriging
plot(variogram(Abundance ~ 1, AK.abund)) # variogram


# SP plot
rv = list("sp.polygons", krig.Key, fill = "chartreuse3", alpha = 0.1)
text1 = list("sp.text", c(-155.31,19.818), "0", cex = .5, which = 1)
text2 = list("sp.text", c(-155.305,19.818), "100 m", cex = .5, which = 1)
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-155.31, 19.820), scale = 50, fill=c("transparent","black"), which = 1)

# levels for overlap
hab.cols<-c("gray20", "gray70")
levs<-as.factor(AK.abund$HabitatType)

spl <- list('sp.points', AK.abund, cex=0.5, pch=21, col=c("gray20", "gray70")[levs])

plot.krig.Key<- spplot(krig.Key["var1.pred"], col.regions=colorRampPalette(col.scheme.N), 
                       sp.layout=spl, col=NA, #  main="Restored Forest Keystone AM-fungi", cex.main=0.8,
                       scales=list(draw = TRUE), xlab= "Longitude", ylab="Latitude", cex.axis=0.8)

plot.krig.Key 
update(plot.krig.Key, key=simpleKey(c("Restored Forest", "Remnant Forest"), points=FALSE, columns=1, cex.main=0.6, cex=0.8, col=c("gray70", "gray20"), space='top'))

dev.copy(jpeg, "figures/plot.krig.test.jpeg", height=600, width=700)
dev.off()



##############################################################################################################
### Co-occurrence networks by plot (habitat type x plot) for overall network characteristics ###
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


#########################################
#####  CHARACTERISTICS  of Habitats #####
#########################################
### By habitat VTs

# Centrality-Betweeness
ig.RO_between <- igraph::betweenness(ig.RO, weights=E(ig.RO),normalized=TRUE)
ig.AK_between <- igraph::betweenness(ig.AK,weights=E(ig.AK),normalized=TRUE)

#Connectedness (degree)
ig.RO_degree <- igraph::degree(ig.RO, mode="all", normalized=TRUE)
ig.AK_degree <- igraph::degree(ig.AK, mode="all", normalized=TRUE)

# Average Path Length: how long are paths around a network?
average.path.length(ig.RO)
average.path.length(ig.AK)


#######################################################
#####  CHARACTERISTICS  of Plot means in Habitats #####
#######################################################

### For each plot, then generate means
# Centrality-Betweeness
ig.RO1_between <- igraph::betweenness(ig.RO1, weights=E(ig.RO1), normalized=TRUE)
ig.RO2_between <- igraph::betweenness(ig.RO2, weights=E(ig.RO2), normalized=TRUE)
ig.RO3_between <- igraph::betweenness(ig.RO3, weights=E(ig.RO3), normalized=TRUE)
ig.RO4_between <- igraph::betweenness(ig.RO4, weights=E(ig.RO4), normalized=TRUE)
ig.RO5_between <- igraph::betweenness(ig.RO5, weights=E(ig.RO5), normalized=TRUE)
ig.RO6_between <- igraph::betweenness(ig.RO6, weights=E(ig.RO6), normalized=TRUE)

ig.AK1_between <- igraph::betweenness(ig.AK1, weights=E(ig.AK1), normalized=TRUE)
ig.AK2_between <- igraph::betweenness(ig.AK2, weights=E(ig.AK2), normalized=TRUE)
ig.AK3_between <- igraph::betweenness(ig.AK3, weights=E(ig.AK3), normalized=TRUE)
ig.AK4_between <- igraph::betweenness(ig.AK4, weights=E(ig.AK4), normalized=TRUE)
ig.AK5_between <- igraph::betweenness(ig.AK5, weights=E(ig.AK5), normalized=TRUE)
ig.AK6_between <- igraph::betweenness(ig.AK6, weights=E(ig.AK6), normalized=TRUE)

haka_centrality <- c(mean(ig.RO1_between),mean(ig.RO2_between),mean(ig.RO3_between),
                        mean(ig.RO4_between),mean(ig.RO5_between),mean(ig.RO6_between),
                        mean(ig.AK1_between),mean(ig.AK2_between),mean(ig.AK3_between),
                        mean(ig.AK4_between),mean(ig.AK5_between),mean(ig.AK6_between))

haka_centrality <- as.data.frame(haka_centrality)


#Connectedness (degree)
ig.RO1_degree <- igraph::degree(ig.RO1, mode="all", normalized=TRUE)
ig.RO2_degree <- igraph::degree(ig.RO2, mode="all", normalized=TRUE)
ig.RO3_degree <- igraph::degree(ig.RO3, mode="all", normalized=TRUE)
ig.RO4_degree <- igraph::degree(ig.RO4, mode="all", normalized=TRUE)
ig.RO5_degree <- igraph::degree(ig.RO5, mode="all", normalized=TRUE)
ig.RO6_degree <- igraph::degree(ig.RO6, mode="all", normalized=TRUE)

ig.AK1_degree <- igraph::degree(ig.AK1, mode="all", normalized=TRUE)
ig.AK2_degree <- igraph::degree(ig.AK2, mode="all", normalized=TRUE)
ig.AK3_degree <- igraph::degree(ig.AK3, mode="all", normalized=TRUE)
ig.AK4_degree <- igraph::degree(ig.AK4, mode="all", normalized=TRUE)
ig.AK5_degree <- igraph::degree(ig.AK5, mode="all", normalized=TRUE)
ig.AK6_degree <- igraph::degree(ig.AK6, mode="all", normalized=TRUE)

haka_connectedness <- c(mean(ig.RO1_degree),mean(ig.RO2_degree),mean(ig.RO3_degree),
                        mean(ig.RO4_degree),mean(ig.RO5_degree),mean(ig.RO6_degree),
                        mean(ig.AK1_degree),mean(ig.AK2_degree),mean(ig.AK3_degree),
                        mean(ig.AK4_degree),mean(ig.AK5_degree),mean(ig.AK6_degree))
        
haka_connectedness <- as.data.frame(haka_connectedness)


## Path length average.path.length(ig.RO1)
ig.RO1_path <- igraph::average.path.length(ig.RO1)
ig.RO2_path <- igraph::average.path.length(ig.RO2)
ig.RO3_path <- igraph::average.path.length(ig.RO3)
ig.RO4_path <- igraph::average.path.length(ig.RO4)
ig.RO5_path <- igraph::average.path.length(ig.RO5)
ig.RO6_path <- igraph::average.path.length(ig.RO6)

ig.AK1_path <- igraph::average.path.length(ig.AK1)
ig.AK2_path <- igraph::average.path.length(ig.AK2)
ig.AK3_path <- igraph::average.path.length(ig.AK3)
ig.AK4_path <- igraph::average.path.length(ig.AK4)
ig.AK5_path <- igraph::average.path.length(ig.AK5)
ig.AK6_path <- igraph::average.path.length(ig.AK6)

haka_path <- c(mean(ig.RO1_path),mean(ig.RO2_path),mean(ig.RO3_path),
                        mean(ig.RO4_path),mean(ig.RO5_path),mean(ig.RO6_path),
                        mean(ig.AK1_path),mean(ig.AK2_path),mean(ig.AK3_path),
                        mean(ig.AK4_path),mean(ig.AK5_path),mean(ig.AK6_path))

haka_path <- as.data.frame(haka_path)

# Dataframe building
plot<-as.data.frame(c("RO1","RO2","RO3","RO4","RO5","RO6",
                      "AK1","AK2","AK3","AK4","AK5","AK6"))
hab_type <- as.data.frame(rep(c("Remnant Forest","Restored Forest"),each=6))
sample_type<- as.data.frame(rep(c("soil"),each=6,times=2))

# add in MEAN traits for each plot, in each network
fungal_networks <- cbind(plot,hab_type,sample_type,haka_centrality, haka_connectedness, haka_path)
colnames(fungal_networks) <- c("Plot","HabitatType","SampleType","Centrality","Connectedness", "Path.length")

fungal_networks %>%
    group_by(HabitatType) %>%
    summarize(means = mean(Centrality), se= se(Centrality))

fungal_networks %>%
    group_by(HabitatType) %>%
    summarize(means = mean(Connectedness), se= se(Connectedness))

fungal_networks %>%
    group_by(HabitatType) %>%
    summarize(means = mean(Path.length), se= se(Path.length))

## Welch t-tests
# Centrality by plots
centrality_habitat <- t.test(subset(fungal_networks, HabitatType == "Remnant Forest")$Centrality,
                                subset(fungal_networks, HabitatType == "Restored Forest")$Centrality,
                                paired=FALSE, var.equal=FALSE)
centrality_habitat

# Connectedness by plots
connectedness_habitat <- t.test(subset(fungal_networks, HabitatType == "Remnant Forest")$Connectedness,
                                subset(fungal_networks, HabitatType == "Restored Forest")$Connectedness,
                                paired=FALSE, var.equal=FALSE)
connectedness_habitat


# Path Length by plots
path_habitat <- t.test(subset(fungal_networks, HabitatType == "Remnant Forest")$Path.length,
                                subset(fungal_networks, HabitatType == "Restored Forest")$Path.length,
                                paired=FALSE, var.equal=FALSE)
path_habitat

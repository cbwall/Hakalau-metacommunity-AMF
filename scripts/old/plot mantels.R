# entire dataset
haka_meta # metadata
hak_otu<-t(otu) # transpose


#separate dataframes
haka_met_mant<-Study.dat[,c(1:32)] # the meta data
hak_otu<-Study.dat[,c(33:1468)] # OTU data

### distance from plots (coded)
plot.dist # all dist
plot.m.dist # mean dist
plot.m.dist$Transect<-c(8,8,8,9,9,9,8,8,8,9,9,9) #add Transect

# set lat-long to mean value for full dataframe


haka_meta$Lat.mean<-ifelse(haka_meta$Plot=="AK1", plot.m.dist[1,3],
                    ifelse(haka_meta$Plot=="AK2", plot.m.dist[2,3],
                    ifelse(haka_meta$Plot=="AK3", plot.m.dist[3,3],
                    ifelse(haka_meta$Plot=="AK4", plot.m.dist[4,3],
                    ifelse(haka_meta$Plot=="AK5", plot.m.dist[5,3],
                    ifelse(haka_meta$Plot=="AK6", plot.m.dist[6,3],
                    ifelse(haka_meta$Plot=="RO1", plot.m.dist[7,3],
                    ifelse(haka_meta$Plot=="RO2", plot.m.dist[8,3],
                    ifelse(haka_meta$Plot=="RO3", plot.m.dist[9,3],
                    ifelse(haka_meta$Plot=="RO4", plot.m.dist[10,3],
                    ifelse(haka_meta$Plot=="RO5", plot.m.dist[11,3],
                     plot.m.dist[12,3])  # last for RO6
                        ))))))))))

haka_meta$Long.mean<-ifelse(haka_meta$Plot=="AK1", plot.m.dist[1,4],
                     ifelse(haka_meta$Plot=="AK2", plot.m.dist[2,4],
                     ifelse(haka_meta$Plot=="AK3", plot.m.dist[3,4],
                     ifelse(haka_meta$Plot=="AK4", plot.m.dist[4,4],
                     ifelse(haka_meta$Plot=="AK5", plot.m.dist[5,4],
                     ifelse(haka_meta$Plot=="AK6", plot.m.dist[6,4],
                     ifelse(haka_meta$Plot=="RO1", plot.m.dist[7,4],
                     ifelse(haka_meta$Plot=="RO2", plot.m.dist[8,4],
                     ifelse(haka_meta$Plot=="RO3", plot.m.dist[9,4],
                     ifelse(haka_meta$Plot=="RO4", plot.m.dist[10,4],
                     ifelse(haka_meta$Plot=="RO5", plot.m.dist[11,4],
                         plot.m.dist[12,4]) # last for RO6
                                  ))))))))))


########## separate distance matrix transect
##### Transect 8

haka_met_mant.T8<-haka_meta[(haka_meta$Transect==8),]

########### separate OTUs accordingly
## subset the OTU file to give only T8 data = AK 1-3 and RO 1-3
# transposed data with only Transect 8
T8dat<-t(otu[, grepl("^AK1", colnames(otu)) | 
               grepl("^AK2", colnames(otu)) |
               grepl("^AK3", colnames(otu)) |
               grepl("^RO1", colnames(otu)) | 
               grepl("^RO2", colnames(otu)) |
               grepl("^RO3", colnames(otu))])

# only rows from T8
samples.T8<-sample[grepl("^AK1", rownames(sample)) |
                     grepl("^AK2", rownames(sample)) |
                     grepl("^AK3", rownames(sample)) |
                     grepl("^RO1", rownames(sample)) |
                     grepl("^RO2", rownames(sample)) |
                     grepl("^RO3", rownames(sample)),]

all.equal(rownames(T8dat), rownames(samples.T8)) # rows match
T8.df<- merge(haka_met_mant.T8, T8dat, by = "row.names", all = TRUE) # all data merged so that row names match
T8.df<-na.omit(T8.df) # remove NAs

T8_meta<-T8.df[,c(1:34)] # the meta data
T8_otu<-T8.df[,c(35:1470)] # OTU data

# Distance to neighboring plots
g_dists_km.T8<- earth.dist(T8_meta[,33:34],dist=TRUE)
g_dists_m.T8<- (g_dists_km.T8 * 1000)
g_dists_m.T8 <- as.matrix(g_dists_m.T8)
rownames(g_dists_m.T8) <- rownames(T8_meta)
colnames(g_dists_m.T8)<- t(rownames(T8_meta))

# rename and make distance matrix
geomat.T8<-g_dists_m.T8
g_dists_m.T8<-as.dist(g_dists_m.T8) 


######## Transect 9
haka_met_mant.T9<-haka_meta[(haka_meta$Transect==9),]

##### subset the OTU file to give only T9 data = AK 4-6 and RO 4-6
# transposed data with only Transect 9
T9dat<-t(otu[, grepl("^AK4", colnames(otu)) | 
               grepl("^AK5", colnames(otu)) |
               grepl("^AK6", colnames(otu)) |
               grepl("^RO4", colnames(otu)) | 
               grepl("^RO5", colnames(otu)) |
               grepl("^RO6", colnames(otu))])

# only rows from T9
samples.T9<-sample[grepl("^AK4", rownames(sample)) |
                     grepl("^AK5", rownames(sample)) |
                     grepl("^AK6", rownames(sample)) |
                     grepl("^RO4", rownames(sample)) |
                     grepl("^RO5", rownames(sample)) |
                     grepl("^RO6", rownames(sample)),]

all.equal(rownames(T9dat), rownames(samples.T9)) # rows match
T9.df<- merge(haka_met_mant.T9, T9dat, by = "row.names", all = TRUE, na.rm=T) # all data merged so that row names match
T9.df<-na.omit(T9.df) # remove NAs
T9_meta<-T9.df[,c(1:34)] # the meta data
T9_otu<-T9.df[,c(35:1470)] # OTU data


# Distance to neighboring plots
g_dists_km.T9<- earth.dist(T9_meta[,33:34],dist=TRUE)
g_dists_m.T9<- (g_dists_km.T9 * 1000)
g_dists_m.T9 <- as.matrix(g_dists_m.T9)
rownames(g_dists_m.T9) <- rownames(T9_meta)
colnames(g_dists_m.T9)<- t(rownames(T9_meta))

# rename and make distance matrix
geomat.T9<-g_dists_m.T9
g_dists_m.T9<-as.dist(g_dists_m.T9) 



############################
############################
#Calculate Bray-Curtis distances on T8 and T9 data
T8_bc_dist = as.dist((vegdist(T8_otu, "bray"))); T8_bc_mat<-as.matrix(vegdist(T8_otu, "bray"))
T9_bc_dist = as.dist((vegdist(T9_otu, "bray"))); T9_bc_mat<-as.matrix(vegdist(T9_otu, "bray"))


# Mantel Test
T8_mantel_test<-mantel.rtest(g_dists_m.T8, T8_bc_dist, nrepet=9999)
T8_mantel_test

T9_mantel_test<-mantel.rtest(g_dists_m.T9, T9_bc_dist, nrepet=9999)
T9_mantel_test

# make df
T8_bcdf<-data.frame(Distance=geomat.T8[lower.tri(geomat.T8)],
                  BrayCurtis= T8_bc_mat[lower.tri(T8_bc_mat)])
summary(T8_bcdf)

# make df
T9_bcdf<-data.frame(Distance=geomat.T9[lower.tri(geomat.T9)],
                    BrayCurtis= T9_bc_mat[lower.tri(T9_bc_mat)])
summary(T9_bcdf)

#### mantel modesl
mantel.mod.T8<-lm(BrayCurtis~log(Distance+1), data=T8_bcdf) # linear model fit
mantel.mod.T9<-lm(BrayCurtis~log(Distance+1), data=T9_bcdf) # linear model fit


### plots
T8_dist_plot<- ggplot(T8_bcdf,aes(x=log(Distance+1), y=BrayCurtis)) +
  geom_point(size=1,alpha=0.5) +
  theme(text=element_text(colour="black",size=12)) + 
  ggtitle("transect 8") +
  stat_smooth(method = "lm", size = 1, se=F) + 
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=12)) +
  theme(axis.text.y=element_text(colour="black",size=12)) +
  theme(legend.title=element_text(colour="black",size=12,face="bold")) +
  theme(legend.text=element_text(colour="black",size=12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

plot(T8_dist_plot)
ggsave("figures/T8_dist_plot.pdf", 
       plot = T8_dist_plot, width = 5, height = 5)


#### T9 plot
T9_dist_plot<- ggplot(T9_bcdf,aes(x=log(Distance+1), y=BrayCurtis)) +
  geom_point(size=1,alpha=0.5) +
  theme(text=element_text(colour="black",size=12)) + 
  stat_smooth(method = "lm", size = 1, se=F) + 
  ggtitle("transect 9") +
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=12)) +
  theme(axis.text.y=element_text(colour="black",size=12)) +
  theme(legend.title=element_text(colour="black",size=12,face="bold")) +
  theme(legend.text=element_text(colour="black",size=12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

plot(T9_dist_plot)
ggsave("figures/T9_dist_plot.pdf", 
       plot = T9_dist_plot, width = 5, height = 5)


##################
####################################
######################################################
# what about plot-plot comparisons within a transect. 
# how does each RO compare to each AK within a transect
# Use AK close, AK further, AK furthest 
# compare to close/further/furthest RO
# generate Bray-Curtis dissimiliarity
######################################################
####################################
##################

# Transect 8 : South 
# plot distance: W --> E (AK1, AK2, AK3, RO2, RO3, RO1)

# Transect 9 : North 
# plot distance: W --> E (AK4, AK5, AK6, RO4, RO5, RO6)

# T8 dat = transposed OTUs
# sample.T8 = samples
# T8.df = dataframe 
# T8_meta = the meta data
# T8_otu = otus

T8dat<-t(otu[, grepl("^AK1", colnames(otu)) | 
               grepl("^AK2", colnames(otu)) |
               grepl("^AK3", colnames(otu)) |
               grepl("^RO1", colnames(otu)) | 
               grepl("^RO2", colnames(otu)) |
               grepl("^RO3", colnames(otu))])

# only rows from T8
samples.T8<-sample[grepl("^AK1", rownames(sample)) |
                     grepl("^AK2", rownames(sample)) |
                     grepl("^AK3", rownames(sample)) |
                     grepl("^RO1", rownames(sample)) |
                     grepl("^RO2", rownames(sample)) |
                     grepl("^RO3", rownames(sample)),]

all.equal(rownames(T8dat), rownames(samples.T8)) # rows match
T8.df<- merge(haka_met_mant.T8, T8dat, by = "row.names", all = TRUE) # all data merged so that row names match
T8.df<-na.omit(T8.df) # remove NAs

T8_meta<-T8.df[,c(1:34)] # the meta data
T8_otu<-T8.df[,c(35:1470)] # OTU data

## T8 single plot dataframes
T8.AK1.df<-T8.df[(T8.df$Plot=="AK1"),]; T8.AK1.meta<- T8.AK1.df[,c(1:34)]; T8.AK1.otu<-T8.AK1.df[,c(35:1470)]
T8.AK2.df<-T8.df[(T8.df$Plot=="AK2"),]; T8.AK2.meta<- T8.AK2.df[,c(1:34)]; T8.AK2.otu<-T8.AK2.df[,c(35:1470)]
T8.AK3.df<-T8.df[(T8.df$Plot=="AK3"),]; T8.AK3.meta<- T8.AK3.df[,c(1:34)]; T8.AK3.otu<-T8.AK3.df[,c(35:1470)]
T8.RO1.df<-T8.df[(T8.df$Plot=="RO1"),]; T8.RO1.meta<- T8.RO1.df[,c(1:34)]; T8.RO1.otu<-T8.RO1.df[,c(35:1470)]
T8.RO2.df<-T8.df[(T8.df$Plot=="RO2"),]; T8.RO2.meta<- T8.RO2.df[,c(1:34)]; T8.RO2.otu<-T8.RO2.df[,c(35:1470)]
T8.RO3.df<-T8.df[(T8.df$Plot=="RO3"),]; T8.RO3.meta<- T8.RO3.df[,c(1:34)]; T8.RO3.otu<-T8.RO3.df[,c(35:1470)]



# Distance to neighboring plots
g_dists_km.T8<- earth.dist(T8_meta[,33:34],dist=TRUE)
g_dists_m.T8<- (g_dists_km.T8 * 1000)
g_dists_m.T8 <- as.matrix(g_dists_m.T8)
rownames(g_dists_m.T8) <- rownames(T8_meta)
colnames(g_dists_m.T8)<- t(rownames(T8_meta))

# rename and make distance matrix
geomat.T8<-g_dists_m.T8
g_dists_m.T8<-as.dist(g_dists_m.T8) 



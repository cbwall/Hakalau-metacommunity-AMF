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
  geom_smooth(method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')), lwd=0.7) +
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=7), 
        axis.text.y=element_text(colour="black",size=7),
        axis.title=element_text(colour="black",size=9),
        legend.title=element_text(colour="black",size=7,face="bold"),
        legend.text=element_text(colour="black",size=7)) +
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
  geom_smooth(method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')), lwd=0.7) +
  ggtitle("transect 9") +
  scale_y_continuous(name="Bray Curtis Dissimilarity",breaks=seq(0,1,0.25),limits=c(0,1)) +
  scale_x_continuous(name="log(Distance (m)+1)",breaks=seq(0,8,1),limits=c(0,8)) +
  theme(axis.text.x=element_text(colour="black",size=7), 
        axis.text.y=element_text(colour="black",size=7),
        axis.title=element_text(colour="black",size=9),
        legend.title=element_text(colour="black",size=7,face="bold"),
        legend.text=element_text(colour="black",size=7)) +
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

# Transect 9 : North 
# plot distance: W --> E (AK4, AK5, AK6, RO4, RO5, RO6)

# Transect 8 : South 
# plot distance: W --> E (AK1, AK2, AK3, RO1, RO2, RO3)


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

# combine for comparisons
T8.RO1.AK1.df<-rbind(T8.RO1.df, T8.AK1.df)
T8.RO1.AK1.meta <-rbind(T8.RO1.meta,T8.AK1.meta)
T8.RO1.AK1.otu<-rbind(T8.RO1.otu,T8.AK1.otu)

T8.RO1.AK2.df<-rbind(T8.RO1.df, T8.AK2.df)
T8.RO1.AK2.meta <-rbind(T8.RO1.meta,T8.AK2.meta)
T8.RO1.AK2.otu<-rbind(T8.RO1.otu,T8.AK2.otu)

T8.RO1.AK3.df<-rbind(T8.RO1.df, T8.AK3.df)
T8.RO1.AK3.meta <-rbind(T8.RO1.meta,T8.AK3.meta)
T8.RO1.AK3.otu<-rbind(T8.RO1.otu,T8.AK3.otu)

T8.RO2.AK1.df<-rbind(T8.RO2.df, T8.AK1.df)
T8.RO2.AK1.meta <-rbind(T8.RO2.meta,T8.AK1.meta)
T8.RO2.AK1.otu<-rbind(T8.RO2.otu,T8.AK1.otu)

T8.RO2.AK2.df<-rbind(T8.RO2.df, T8.AK2.df)
T8.RO2.AK2.meta <-rbind(T8.RO2.meta,T8.AK2.meta)
T8.RO2.AK2.otu<-rbind(T8.RO2.otu,T8.AK2.otu)

T8.RO2.AK3.df<-rbind(T8.RO2.df, T8.AK3.df)
T8.RO2.AK3.meta <-rbind(T8.RO2.meta,T8.AK3.meta)
T8.RO2.AK3.otu<-rbind(T8.RO2.otu,T8.AK3.otu)

T8.RO3.AK1.df<-rbind(T8.RO3.df, T8.AK1.df)
T8.RO3.AK1.meta <-rbind(T8.RO3.meta,T8.AK1.meta)
T8.RO3.AK1.otu<-rbind(T8.RO3.otu,T8.AK1.otu)

T8.RO3.AK2.df<-rbind(T8.RO3.df, T8.AK2.df)
T8.RO3.AK2.meta <-rbind(T8.RO3.meta,T8.AK2.meta)
T8.RO3.AK2.otu<-rbind(T8.RO3.otu,T8.AK2.otu)

T8.RO3.AK3.df<-rbind(T8.RO3.df, T8.AK3.df)
T8.RO3.AK3.meta <-rbind(T8.RO3.meta,T8.AK3.meta)
T8.RO3.AK3.otu<-rbind(T8.RO3.otu,T8.AK3.otu)


# Distance to neighboring plots
######### RO1
# RO1 - AK1
g_dists_km.T8.RO1.AK1<- earth.dist(T8.RO1.AK1.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO1.AK1<- (g_dists_km.T8.RO1.AK1 * 1000)
g_dists_m.T8.RO1.AK1 <- as.matrix(g_dists_m.T8.RO1.AK1)
rownames(g_dists_m.T8.RO1.AK1) <- rownames(T8.RO1.AK1.meta)
colnames(g_dists_m.T8.RO1.AK1)<- t(rownames(T8.RO1.AK1.meta))
# rename and make distance matrix
geomat.T8.RO1.AK1<-g_dists_m.T8.RO1.AK1
g_dists_m.T8.RO1.AK1<-as.dist(g_dists_m.T8.RO1.AK1) 


#### RO1 - AK2
g_dists_km.T8.RO1.AK2<- earth.dist(T8.RO1.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO1.AK2<- (g_dists_km.T8.RO1.AK2 * 1000)
g_dists_m.T8.RO1.AK2 <- as.matrix(g_dists_m.T8.RO1.AK2)
rownames(g_dists_m.T8.RO1.AK2) <- rownames(T8.RO1.AK2.meta)
colnames(g_dists_m.T8.RO1.AK2)<- t(rownames(T8.RO1.AK2.meta))
# rename and make distance matrix
geomat.T8.RO1.AK2<-g_dists_m.T8.RO1.AK2
g_dists_m.T8.RO1.AK2<-as.dist(g_dists_m.T8.RO1.AK2) 


#### RO1 - AK3
g_dists_km.T8.RO1.AK3<- earth.dist(T8.RO1.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO1.AK3<- (g_dists_km.T8.RO1.AK3 * 1000)
g_dists_m.T8.RO1.AK3 <- as.matrix(g_dists_m.T8.RO1.AK3)
rownames(g_dists_m.T8.RO1.AK3) <- rownames(T8.RO1.AK3.meta)
colnames(g_dists_m.T8.RO1.AK3)<- t(rownames(T8.RO1.AK3.meta))
# rename and make distance matrix
geomat.T8.RO1.AK3<-g_dists_m.T8.RO1.AK3
g_dists_m.T8.RO1.AK3<-as.dist(g_dists_m.T8.RO1.AK3) 

################ RO2
# RO2 - AK1
g_dists_km.T8.RO2.AK1<- earth.dist(T8.RO2.AK1.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO2.AK1<- (g_dists_km.T8.RO2.AK1 * 1000)
g_dists_m.T8.RO2.AK1 <- as.matrix(g_dists_m.T8.RO2.AK1)
rownames(g_dists_m.T8.RO2.AK1) <- rownames(T8.RO2.AK1.meta)
colnames(g_dists_m.T8.RO2.AK1)<- t(rownames(T8.RO2.AK1.meta))
# rename and make distance matrix
geomat.T8.RO2.AK1<-g_dists_m.T8.RO2.AK1
g_dists_m.T8.RO2.AK1<-as.dist(g_dists_m.T8.RO2.AK1) 


#### RO2 - AK2
g_dists_km.T8.RO2.AK2<- earth.dist(T8.RO2.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO2.AK2<- (g_dists_km.T8.RO2.AK2 * 1000)
g_dists_m.T8.RO2.AK2 <- as.matrix(g_dists_m.T8.RO2.AK2)
rownames(g_dists_m.T8.RO2.AK2) <- rownames(T8.RO2.AK2.meta)
colnames(g_dists_m.T8.RO2.AK2)<- t(rownames(T8.RO2.AK2.meta))
# rename and make distance matrix
geomat.T8.RO2.AK2<-g_dists_m.T8.RO2.AK2
g_dists_m.T8.RO2.AK2<-as.dist(g_dists_m.T8.RO2.AK2) 


#### RO2 - AK3
g_dists_km.T8.RO2.AK3<- earth.dist(T8.RO2.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO2.AK3<- (g_dists_km.T8.RO2.AK3 * 1000)
g_dists_m.T8.RO2.AK3 <- as.matrix(g_dists_m.T8.RO2.AK3)
rownames(g_dists_m.T8.RO2.AK3) <- rownames(T8.RO2.AK3.meta)
colnames(g_dists_m.T8.RO2.AK3)<- t(rownames(T8.RO2.AK3.meta))
# rename and make distance matrix
geomat.T8.RO2.AK3<-g_dists_m.T8.RO2.AK3
g_dists_m.T8.RO2.AK3<-as.dist(g_dists_m.T8.RO2.AK3) 

################ RO3
# RO3 - AK1
g_dists_km.T8.RO3.AK1<- earth.dist(T8.RO3.AK1.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO3.AK1<- (g_dists_km.T8.RO3.AK1 * 1000)
g_dists_m.T8.RO3.AK1 <- as.matrix(g_dists_m.T8.RO3.AK1)
rownames(g_dists_m.T8.RO3.AK1) <- rownames(T8.RO3.AK1.meta)
colnames(g_dists_m.T8.RO3.AK1)<- t(rownames(T8.RO3.AK1.meta))
# rename and make distance matrix
geomat.T8.RO3.AK1<-g_dists_m.T8.RO3.AK1
g_dists_m.T8.RO3.AK1<-as.dist(g_dists_m.T8.RO3.AK1) 


#### RO3 - AK2
g_dists_km.T8.RO3.AK2<- earth.dist(T8.RO3.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO3.AK2<- (g_dists_km.T8.RO3.AK2 * 1000)
g_dists_m.T8.RO3.AK2 <- as.matrix(g_dists_m.T8.RO3.AK2)
rownames(g_dists_m.T8.RO3.AK2) <- rownames(T8.RO3.AK2.meta)
colnames(g_dists_m.T8.RO3.AK2)<- t(rownames(T8.RO3.AK2.meta))
# rename and make distance matrix
geomat.T8.RO3.AK2<-g_dists_m.T8.RO3.AK2
g_dists_m.T8.RO3.AK2<-as.dist(g_dists_m.T8.RO3.AK2) 


#### RO3 - AK3
g_dists_km.T8.RO3.AK3<- earth.dist(T8.RO3.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.RO3.AK3<- (g_dists_km.T8.RO3.AK3 * 1000)
g_dists_m.T8.RO3.AK3 <- as.matrix(g_dists_m.T8.RO3.AK3)
rownames(g_dists_m.T8.RO3.AK3) <- rownames(T8.RO3.AK3.meta)
colnames(g_dists_m.T8.RO3.AK3)<- t(rownames(T8.RO3.AK3.meta))
# rename and make distance matrix
geomat.T8.RO3.AK3<-g_dists_m.T8.RO3.AK3
g_dists_m.T8.RO3.AK3<-as.dist(g_dists_m.T8.RO3.AK3) 


############################
#Calculate Bray-Curtis distances on T8 between RO to AK

#### RO1
T8_RO1.AK1.bc_dist = as.dist((vegdist(T8.RO1.AK1.otu, "bray")))
T8_RO1.AK1.bc_mat<-as.matrix(vegdist(T8.RO1.AK1.otu, "bray"))

T8_RO1.AK2.bc_dist = as.dist((vegdist(T8.RO1.AK2.otu, "bray")))
T8_RO1.AK2.bc_mat<-as.matrix(vegdist(T8.RO1.AK2.otu, "bray"))

T8_RO1.AK3.bc_dist = as.dist((vegdist(T8.RO1.AK3.otu, "bray")))
T8_RO1.AK3.bc_mat<-as.matrix(vegdist(T8.RO1.AK3.otu, "bray"))

#### RO2
T8_RO2.AK1.bc_dist = as.dist((vegdist(T8.RO2.AK1.otu, "bray")))
T8_RO2.AK1.bc_mat<-as.matrix(vegdist(T8.RO2.AK1.otu, "bray"))

T8_RO2.AK2.bc_dist = as.dist((vegdist(T8.RO2.AK2.otu, "bray")))
T8_RO2.AK2.bc_mat<-as.matrix(vegdist(T8.RO2.AK2.otu, "bray"))

T8_RO2.AK3.bc_dist = as.dist((vegdist(T8.RO2.AK3.otu, "bray")))
T8_RO2.AK3.bc_mat<-as.matrix(vegdist(T8.RO2.AK3.otu, "bray"))

#### RO3
T8_RO3.AK1.bc_dist = as.dist((vegdist(T8.RO3.AK1.otu, "bray")))
T8_RO3.AK1.bc_mat<-as.matrix(vegdist(T8.RO3.AK1.otu, "bray"))

T8_RO3.AK2.bc_dist = as.dist((vegdist(T8.RO3.AK2.otu, "bray")))
T8_RO3.AK2.bc_mat<-as.matrix(vegdist(T8.RO3.AK2.otu, "bray"))

T8_RO3.AK3.bc_dist = as.dist((vegdist(T8.RO3.AK3.otu, "bray")))
T8_RO3.AK3.bc_mat<-as.matrix(vegdist(T8.RO3.AK3.otu, "bray"))

#################
# Mantel Test
T8.RO1.AK1_mantel_test<-mantel.rtest(g_dists_m.T8.RO1.AK1, T8_RO1.AK1.bc_dist, nrepet=9999)
T8.RO1.AK1_mantel_test

T8.RO1.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.RO1.AK2, T8_RO1.AK2.bc_dist, nrepet=9999)
T8.RO1.AK2_mantel_test

T8.RO1.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.RO1.AK3, T8_RO1.AK3.bc_dist, nrepet=9999)
T8.RO1.AK3_mantel_test


### RO2
T8.RO2.AK1_mantel_test<-mantel.rtest(g_dists_m.T8.RO2.AK1, T8_RO2.AK1.bc_dist, nrepet=9999)
T8.RO2.AK1_mantel_test

T8.RO2.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.RO2.AK2, T8_RO2.AK2.bc_dist, nrepet=9999)
T8.RO2.AK2_mantel_test

T8.RO2.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.RO2.AK3, T8_RO2.AK3.bc_dist, nrepet=9999)
T8.RO2.AK3_mantel_test

### RO3
T8.RO3.AK1_mantel_test<-mantel.rtest(g_dists_m.T8.RO3.AK1, T8_RO3.AK1.bc_dist, nrepet=9999)
T8.RO3.AK1_mantel_test

T8.RO3.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.RO3.AK2, T8_RO3.AK2.bc_dist, nrepet=9999)
T8.RO3.AK2_mantel_test

T8.RO3.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.RO3.AK3, T8_RO3.AK3.bc_dist, nrepet=9999)
T8.RO3.AK3_mantel_test

###############
# make df
### RO1
T8_RO1.AK1.bcdf<-data.frame(Distance=geomat.T8.RO1.AK1[lower.tri(geomat.T8.RO1.AK1)],
                    BrayCurtis= T8_RO1.AK1.bc_mat[lower.tri(T8_RO1.AK1.bc_mat)])

T8_RO1.AK2.bcdf<-data.frame(Distance=geomat.T8.RO1.AK2[lower.tri(geomat.T8.RO1.AK2)],
                            BrayCurtis= T8_RO1.AK2.bc_mat[lower.tri(T8_RO1.AK2.bc_mat)])

T8_RO1.AK3.bcdf<-data.frame(Distance=geomat.T8.RO1.AK3[lower.tri(geomat.T8.RO1.AK3)],
                            BrayCurtis= T8_RO1.AK3.bc_mat[lower.tri(T8_RO1.AK3.bc_mat)])

### RO2
T8_RO2.AK1.bcdf<-data.frame(Distance=geomat.T8.RO2.AK1[lower.tri(geomat.T8.RO2.AK1)],
                            BrayCurtis= T8_RO2.AK1.bc_mat[lower.tri(T8_RO2.AK1.bc_mat)])

T8_RO2.AK2.bcdf<-data.frame(Distance=geomat.T8.RO2.AK2[lower.tri(geomat.T8.RO2.AK2)],
                            BrayCurtis= T8_RO2.AK2.bc_mat[lower.tri(T8_RO2.AK2.bc_mat)])

T8_RO2.AK3.bcdf<-data.frame(Distance=geomat.T8.RO2.AK3[lower.tri(geomat.T8.RO2.AK3)],
                            BrayCurtis= T8_RO2.AK3.bc_mat[lower.tri(T8_RO2.AK3.bc_mat)])

### RO3
T8_RO3.AK1.bcdf<-data.frame(Distance=geomat.T8.RO3.AK1[lower.tri(geomat.T8.RO3.AK1)],
                            BrayCurtis= T8_RO3.AK1.bc_mat[lower.tri(T8_RO3.AK1.bc_mat)])

T8_RO3.AK2.bcdf<-data.frame(Distance=geomat.T8.RO3.AK2[lower.tri(geomat.T8.RO3.AK2)],
                            BrayCurtis= T8_RO3.AK2.bc_mat[lower.tri(T8_RO3.AK2.bc_mat)])

T8_RO3.AK3.bcdf<-data.frame(Distance=geomat.T8.RO3.AK3[lower.tri(geomat.T8.RO3.AK3)],
                            BrayCurtis= T8_RO3.AK3.bc_mat[lower.tri(T8_RO3.AK3.bc_mat)])

# mantel test df
prox.mantT8<-as.data.frame(rbind(T8.RO1.AK1_mantel_test[1], T8.RO1.AK2_mantel_test[1],T8.RO1.AK3_mantel_test[1],
             T8.RO2.AK1_mantel_test[1], T8.RO2.AK2_mantel_test[1], T8.RO2.AK3_mantel_test[1],
             T8.RO3.AK1_mantel_test[1], T8.RO3.AK2_mantel_test[1], T8.RO3.AK3_mantel_test[1]))
prox.mantT8$obs<-as.numeric(prox.mantT8$obs)

prox.mantT8$pvalue<-as.numeric(rbind(
  T8.RO1.AK1_mantel_test[5], T8.RO1.AK2_mantel_test[5], T8.RO1.AK3_mantel_test[5],
       T8.RO2.AK1_mantel_test[5], T8.RO2.AK2_mantel_test[5], T8.RO2.AK3_mantel_test[5],
       T8.RO3.AK1_mantel_test[5], T8.RO3.AK2_mantel_test[5], T8.RO3.AK3_mantel_test[5]))
prox.mantT8$comparison<-as.factor(c("T8.RO1.AK1", "T8.RO1.AK2", "T8.RO1.AK3",
                        "T8.RO2.AK1", "T8.RO2.AK2", "T8.RO2.AK3",
                        "T8.RO3.AK1", "T8.RO3.AK2", "T8.RO3.AK3"))
prox.mantT8$AK.group<-as.factor(c("AK far", "AK mid", "AK near", 
                      "AK far", "AK mid", "AK near", 
                      "AK far", "AK mid", "AK near"))
prox.mantT8$Transect<-"8"

################# 
################# 
################# Transect 9
# Transect 9 : North 
# plot distance: W --> E (AK4, AK5, AK6, RO4, RO5, RO6)


# T9 dat = transposed OTUs
# sample.T9 = samples
# T9.df = dataframe 
# T9_meta = the meta data
# T9_otu = otus


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
T9.df<- merge(haka_met_mant.T9, T9dat, by = "row.names", all = TRUE) # all data merged so that row names match
T9.df<-na.omit(T9.df) # remove NAs

T9_meta<-T9.df[,c(1:34)] # the meta data
T9_otu<-T9.df[,c(35:1470)] # OTU data

## T9 single plot dataframes
T9.AK4.df<-T9.df[(T9.df$Plot=="AK4"),]; T9.AK4.meta<- T9.AK4.df[,c(1:34)]; T9.AK4.otu<-T9.AK4.df[,c(35:1470)]
T9.AK5.df<-T9.df[(T9.df$Plot=="AK5"),]; T9.AK5.meta<- T9.AK5.df[,c(1:34)]; T9.AK5.otu<-T9.AK5.df[,c(35:1470)]
T9.AK6.df<-T9.df[(T9.df$Plot=="AK6"),]; T9.AK6.meta<- T9.AK6.df[,c(1:34)]; T9.AK6.otu<-T9.AK6.df[,c(35:1470)]
T9.RO4.df<-T9.df[(T9.df$Plot=="RO4"),]; T9.RO4.meta<- T9.RO4.df[,c(1:34)]; T9.RO4.otu<-T9.RO4.df[,c(35:1470)]
T9.RO5.df<-T9.df[(T9.df$Plot=="RO5"),]; T9.RO5.meta<- T9.RO5.df[,c(1:34)]; T9.RO5.otu<-T9.RO5.df[,c(35:1470)]
T9.RO6.df<-T9.df[(T9.df$Plot=="RO6"),]; T9.RO6.meta<- T9.RO6.df[,c(1:34)]; T9.RO6.otu<-T9.RO6.df[,c(35:1470)]

# combine for comparisons
T9.RO4.AK4.df<-rbind(T9.RO4.df, T9.AK4.df)
T9.RO4.AK4.meta <-rbind(T9.RO4.meta,T9.AK4.meta)
T9.RO4.AK4.otu<-rbind(T9.RO4.otu,T9.AK4.otu)

T9.RO4.AK5.df<-rbind(T9.RO4.df, T9.AK5.df)
T9.RO4.AK5.meta <-rbind(T9.RO4.meta,T9.AK5.meta)
T9.RO4.AK5.otu<-rbind(T9.RO4.otu,T9.AK5.otu)

T9.RO4.AK6.df<-rbind(T9.RO4.df, T9.AK6.df)
T9.RO4.AK6.meta <-rbind(T9.RO4.meta,T9.AK6.meta)
T9.RO4.AK6.otu<-rbind(T9.RO4.otu,T9.AK6.otu)

T9.RO5.AK4.df<-rbind(T9.RO5.df, T9.AK4.df)
T9.RO5.AK4.meta <-rbind(T9.RO5.meta,T9.AK4.meta)
T9.RO5.AK4.otu<-rbind(T9.RO5.otu,T9.AK4.otu)

T9.RO5.AK5.df<-rbind(T9.RO5.df, T9.AK5.df)
T9.RO5.AK5.meta <-rbind(T9.RO5.meta,T9.AK5.meta)
T9.RO5.AK5.otu<-rbind(T9.RO5.otu,T9.AK5.otu)

T9.RO5.AK6.df<-rbind(T9.RO5.df, T9.AK6.df)
T9.RO5.AK6.meta <-rbind(T9.RO5.meta,T9.AK6.meta)
T9.RO5.AK6.otu<-rbind(T9.RO5.otu,T9.AK6.otu)

T9.RO6.AK4.df<-rbind(T9.RO6.df, T9.AK4.df)
T9.RO6.AK4.meta <-rbind(T9.RO6.meta,T9.AK4.meta)
T9.RO6.AK4.otu<-rbind(T9.RO6.otu,T9.AK4.otu)

T9.RO6.AK5.df<-rbind(T9.RO6.df, T9.AK5.df)
T9.RO6.AK5.meta <-rbind(T9.RO6.meta,T9.AK5.meta)
T9.RO6.AK5.otu<-rbind(T9.RO6.otu,T9.AK5.otu)

T9.RO6.AK6.df<-rbind(T9.RO6.df, T9.AK6.df)
T9.RO6.AK6.meta <-rbind(T9.RO6.meta,T9.AK6.meta)
T9.RO6.AK6.otu<-rbind(T9.RO6.otu,T9.AK6.otu)


# Distance to neighboring plots
######### RO4
# RO4 - AK4
g_dists_km.T9.RO4.AK4<- earth.dist(T9.RO4.AK4.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO4.AK4<- (g_dists_km.T9.RO4.AK4 * 1000)
g_dists_m.T9.RO4.AK4 <- as.matrix(g_dists_m.T9.RO4.AK4)
rownames(g_dists_m.T9.RO4.AK4) <- rownames(T9.RO4.AK4.meta)
colnames(g_dists_m.T9.RO4.AK4)<- t(rownames(T9.RO4.AK4.meta))
# rename and make distance matrix
geomat.T9.RO4.AK4<-g_dists_m.T9.RO4.AK4
g_dists_m.T9.RO4.AK4<-as.dist(g_dists_m.T9.RO4.AK4) 


#### RO4 - AK5
g_dists_km.T9.RO4.AK5<- earth.dist(T9.RO4.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO4.AK5<- (g_dists_km.T9.RO4.AK5 * 1000)
g_dists_m.T9.RO4.AK5 <- as.matrix(g_dists_m.T9.RO4.AK5)
rownames(g_dists_m.T9.RO4.AK5) <- rownames(T9.RO4.AK5.meta)
colnames(g_dists_m.T9.RO4.AK5)<- t(rownames(T9.RO4.AK5.meta))
# rename and make distance matrix
geomat.T9.RO4.AK5<-g_dists_m.T9.RO4.AK5
g_dists_m.T9.RO4.AK5<-as.dist(g_dists_m.T9.RO4.AK5) 


#### RO4 - AK6
g_dists_km.T9.RO4.AK6<- earth.dist(T9.RO4.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO4.AK6<- (g_dists_km.T9.RO4.AK6 * 1000)
g_dists_m.T9.RO4.AK6 <- as.matrix(g_dists_m.T9.RO4.AK6)
rownames(g_dists_m.T9.RO4.AK6) <- rownames(T9.RO4.AK6.meta)
colnames(g_dists_m.T9.RO4.AK6)<- t(rownames(T9.RO4.AK6.meta))
# rename and make distance matrix
geomat.T9.RO4.AK6<-g_dists_m.T9.RO4.AK6
g_dists_m.T9.RO4.AK6<-as.dist(g_dists_m.T9.RO4.AK6) 

################ RO5
# RO5 - AK4
g_dists_km.T9.RO5.AK4<- earth.dist(T9.RO5.AK4.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO5.AK4<- (g_dists_km.T9.RO5.AK4 * 1000)
g_dists_m.T9.RO5.AK4 <- as.matrix(g_dists_m.T9.RO5.AK4)
rownames(g_dists_m.T9.RO5.AK4) <- rownames(T9.RO5.AK4.meta)
colnames(g_dists_m.T9.RO5.AK4)<- t(rownames(T9.RO5.AK4.meta))
# rename and make distance matrix
geomat.T9.RO5.AK4<-g_dists_m.T9.RO5.AK4
g_dists_m.T9.RO5.AK4<-as.dist(g_dists_m.T9.RO5.AK4) 


#### RO5 - AK5
g_dists_km.T9.RO5.AK5<- earth.dist(T9.RO5.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO5.AK5<- (g_dists_km.T9.RO5.AK5 * 1000)
g_dists_m.T9.RO5.AK5 <- as.matrix(g_dists_m.T9.RO5.AK5)
rownames(g_dists_m.T9.RO5.AK5) <- rownames(T9.RO5.AK5.meta)
colnames(g_dists_m.T9.RO5.AK5)<- t(rownames(T9.RO5.AK5.meta))
# rename and make distance matrix
geomat.T9.RO5.AK5<-g_dists_m.T9.RO5.AK5
g_dists_m.T9.RO5.AK5<-as.dist(g_dists_m.T9.RO5.AK5) 


#### RO5 - AK6
g_dists_km.T9.RO5.AK6<- earth.dist(T9.RO5.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO5.AK6<- (g_dists_km.T9.RO5.AK6 * 1000)
g_dists_m.T9.RO5.AK6 <- as.matrix(g_dists_m.T9.RO5.AK6)
rownames(g_dists_m.T9.RO5.AK6) <- rownames(T9.RO5.AK6.meta)
colnames(g_dists_m.T9.RO5.AK6)<- t(rownames(T9.RO5.AK6.meta))
# rename and make distance matrix
geomat.T9.RO5.AK6<-g_dists_m.T9.RO5.AK6
g_dists_m.T9.RO5.AK6<-as.dist(g_dists_m.T9.RO5.AK6) 

################ RO6
# RO6 - AK4
g_dists_km.T9.RO6.AK4<- earth.dist(T9.RO6.AK4.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO6.AK4<- (g_dists_km.T9.RO6.AK4 * 1000)
g_dists_m.T9.RO6.AK4 <- as.matrix(g_dists_m.T9.RO6.AK4)
rownames(g_dists_m.T9.RO6.AK4) <- rownames(T9.RO6.AK4.meta)
colnames(g_dists_m.T9.RO6.AK4)<- t(rownames(T9.RO6.AK4.meta))
# rename and make distance matrix
geomat.T9.RO6.AK4<-g_dists_m.T9.RO6.AK4
g_dists_m.T9.RO6.AK4<-as.dist(g_dists_m.T9.RO6.AK4) 


#### RO6 - AK5
g_dists_km.T9.RO6.AK5<- earth.dist(T9.RO6.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO6.AK5<- (g_dists_km.T9.RO6.AK5 * 1000)
g_dists_m.T9.RO6.AK5 <- as.matrix(g_dists_m.T9.RO6.AK5)
rownames(g_dists_m.T9.RO6.AK5) <- rownames(T9.RO6.AK5.meta)
colnames(g_dists_m.T9.RO6.AK5)<- t(rownames(T9.RO6.AK5.meta))
# rename and make distance matrix
geomat.T9.RO6.AK5<-g_dists_m.T9.RO6.AK5
g_dists_m.T9.RO6.AK5<-as.dist(g_dists_m.T9.RO6.AK5) 


#### RO6 - AK6
g_dists_km.T9.RO6.AK6<- earth.dist(T9.RO6.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.RO6.AK6<- (g_dists_km.T9.RO6.AK6 * 1000)
g_dists_m.T9.RO6.AK6 <- as.matrix(g_dists_m.T9.RO6.AK6)
rownames(g_dists_m.T9.RO6.AK6) <- rownames(T9.RO6.AK6.meta)
colnames(g_dists_m.T9.RO6.AK6)<- t(rownames(T9.RO6.AK6.meta))
# rename and make distance matrix
geomat.T9.RO6.AK6<-g_dists_m.T9.RO6.AK6
g_dists_m.T9.RO6.AK6<-as.dist(g_dists_m.T9.RO6.AK6) 


############################
#Calculate Bray-Curtis distances on T9 between RO to AK

#### RO4
T9_RO4.AK4.bc_dist = as.dist((vegdist(T9.RO4.AK4.otu, "bray")))
T9_RO4.AK4.bc_mat<-as.matrix(vegdist(T9.RO4.AK4.otu, "bray"))

T9_RO4.AK5.bc_dist = as.dist((vegdist(T9.RO4.AK5.otu, "bray")))
T9_RO4.AK5.bc_mat<-as.matrix(vegdist(T9.RO4.AK5.otu, "bray"))

T9_RO4.AK6.bc_dist = as.dist((vegdist(T9.RO4.AK6.otu, "bray")))
T9_RO4.AK6.bc_mat<-as.matrix(vegdist(T9.RO4.AK6.otu, "bray"))

#### RO5
T9_RO5.AK4.bc_dist = as.dist((vegdist(T9.RO5.AK4.otu, "bray")))
T9_RO5.AK4.bc_mat<-as.matrix(vegdist(T9.RO5.AK4.otu, "bray"))

T9_RO5.AK5.bc_dist = as.dist((vegdist(T9.RO5.AK5.otu, "bray")))
T9_RO5.AK5.bc_mat<-as.matrix(vegdist(T9.RO5.AK5.otu, "bray"))

T9_RO5.AK6.bc_dist = as.dist((vegdist(T9.RO5.AK6.otu, "bray")))
T9_RO5.AK6.bc_mat<-as.matrix(vegdist(T9.RO5.AK6.otu, "bray"))

#### RO6
T9_RO6.AK4.bc_dist = as.dist((vegdist(T9.RO6.AK4.otu, "bray")))
T9_RO6.AK4.bc_mat<-as.matrix(vegdist(T9.RO6.AK4.otu, "bray"))

T9_RO6.AK5.bc_dist = as.dist((vegdist(T9.RO6.AK5.otu, "bray")))
T9_RO6.AK5.bc_mat<-as.matrix(vegdist(T9.RO6.AK5.otu, "bray"))

T9_RO6.AK6.bc_dist = as.dist((vegdist(T9.RO6.AK6.otu, "bray")))
T9_RO6.AK6.bc_mat<-as.matrix(vegdist(T9.RO6.AK6.otu, "bray"))

#################
# Mantel Test
T9.RO4.AK4_mantel_test<-mantel.rtest(g_dists_m.T9.RO4.AK4, T9_RO4.AK4.bc_dist, nrepet=9999)
T9.RO4.AK4_mantel_test

T9.RO4.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.RO4.AK5, T9_RO4.AK5.bc_dist, nrepet=9999)
T9.RO4.AK5_mantel_test

T9.RO4.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.RO4.AK6, T9_RO4.AK6.bc_dist, nrepet=9999)
T9.RO4.AK6_mantel_test


### RO5
T9.RO5.AK4_mantel_test<-mantel.rtest(g_dists_m.T9.RO5.AK4, T9_RO5.AK4.bc_dist, nrepet=9999)
T9.RO5.AK4_mantel_test

T9.RO5.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.RO5.AK5, T9_RO5.AK5.bc_dist, nrepet=9999)
T9.RO5.AK5_mantel_test

T9.RO5.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.RO5.AK6, T9_RO5.AK6.bc_dist, nrepet=9999)
T9.RO5.AK6_mantel_test

### RO6
T9.RO6.AK4_mantel_test<-mantel.rtest(g_dists_m.T9.RO6.AK4, T9_RO6.AK4.bc_dist, nrepet=9999)
T9.RO6.AK4_mantel_test

T9.RO6.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.RO6.AK5, T9_RO6.AK5.bc_dist, nrepet=9999)
T9.RO6.AK5_mantel_test

T9.RO6.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.RO6.AK6, T9_RO6.AK6.bc_dist, nrepet=9999)
T9.RO6.AK6_mantel_test

###############
# make df
### RO4
T9_RO4.AK4.bcdf<-data.frame(Distance=geomat.T9.RO4.AK4[lower.tri(geomat.T9.RO4.AK4)],
                            BrayCurtis= T9_RO4.AK4.bc_mat[lower.tri(T9_RO4.AK4.bc_mat)])

T9_RO4.AK5.bcdf<-data.frame(Distance=geomat.T9.RO4.AK5[lower.tri(geomat.T9.RO4.AK5)],
                            BrayCurtis= T9_RO4.AK5.bc_mat[lower.tri(T9_RO4.AK5.bc_mat)])

T9_RO4.AK6.bcdf<-data.frame(Distance=geomat.T9.RO4.AK6[lower.tri(geomat.T9.RO4.AK6)],
                            BrayCurtis= T9_RO4.AK6.bc_mat[lower.tri(T9_RO4.AK6.bc_mat)])

### RO5
T9_RO5.AK4.bcdf<-data.frame(Distance=geomat.T9.RO5.AK4[lower.tri(geomat.T9.RO5.AK4)],
                            BrayCurtis= T9_RO5.AK4.bc_mat[lower.tri(T9_RO5.AK4.bc_mat)])

T9_RO5.AK5.bcdf<-data.frame(Distance=geomat.T9.RO5.AK5[lower.tri(geomat.T9.RO5.AK5)],
                            BrayCurtis= T9_RO5.AK5.bc_mat[lower.tri(T9_RO5.AK5.bc_mat)])

T9_RO5.AK6.bcdf<-data.frame(Distance=geomat.T9.RO5.AK6[lower.tri(geomat.T9.RO5.AK6)],
                            BrayCurtis= T9_RO5.AK6.bc_mat[lower.tri(T9_RO5.AK6.bc_mat)])

### RO6
T9_RO6.AK4.bcdf<-data.frame(Distance=geomat.T9.RO6.AK4[lower.tri(geomat.T9.RO6.AK4)],
                            BrayCurtis= T9_RO6.AK4.bc_mat[lower.tri(T9_RO6.AK4.bc_mat)])

T9_RO6.AK5.bcdf<-data.frame(Distance=geomat.T9.RO6.AK5[lower.tri(geomat.T9.RO6.AK5)],
                            BrayCurtis= T9_RO6.AK5.bc_mat[lower.tri(T9_RO6.AK5.bc_mat)])

T9_RO6.AK6.bcdf<-data.frame(Distance=geomat.T9.RO6.AK6[lower.tri(geomat.T9.RO6.AK6)],
                            BrayCurtis= T9_RO6.AK6.bc_mat[lower.tri(T9_RO6.AK6.bc_mat)])

# mantel test df
prox.mantT9<-as.data.frame(rbind(
  T9.RO4.AK4_mantel_test[1], T9.RO4.AK5_mantel_test[1],T9.RO4.AK6_mantel_test[1],
  T9.RO5.AK4_mantel_test[1], T9.RO5.AK5_mantel_test[1], T9.RO5.AK6_mantel_test[1],
  T9.RO6.AK4_mantel_test[1], T9.RO6.AK5_mantel_test[1], T9.RO6.AK6_mantel_test[1]))
prox.mantT9$obs<-as.numeric(prox.mantT9$obs)

prox.mantT9$pvalue<-as.numeric(rbind(
  T9.RO4.AK4_mantel_test[5], T9.RO4.AK5_mantel_test[5], T9.RO4.AK6_mantel_test[5],
  T9.RO5.AK4_mantel_test[5], T9.RO5.AK5_mantel_test[5], T9.RO5.AK6_mantel_test[5],
  T9.RO6.AK4_mantel_test[5], T9.RO6.AK5_mantel_test[5], T9.RO6.AK6_mantel_test[5]))

prox.mantT9$comparison<-as.factor(c("T9.RO4.AK4", "T9.RO4.AK5", "T9.RO4.AK6",
                                    "T9.RO5.AK4", "T9.RO5.AK5", "T9.RO5.AK6",
                                    "T9.RO6.AK4", "T9.RO6.AK5", "T9.RO6.AK6"))

prox.mantT9$AK.group<-as.factor(c("AK far", "AK mid", "AK near", 
                                  "AK far", "AK mid", "AK near", 
                                  "AK far", "AK mid", "AK near"))
prox.mantT9$Transect<-"9"

##############
####################### plots

#combined df
prox.mantT8.T9<-rbind(prox.mantT8, prox.mantT9)
prox.mantT8.T9$AK.group<-factor(prox.mantT8.T9$AK.group, levels = c("AK near", "AK mid", "AK far"))

mantel.proximity<-ggplot(prox.mantT8.T9, aes(x=AK.group, y=obs, color=Transect)) + 
  geom_boxplot(aes(color = Transect), width = 0.5, size = 0.4,
               position = position_dodge(0.7)) +
  geom_dotplot(aes(fill = Transect, color = Transect),
    binaxis='y', stackdir='center', dotsize = 0.5, alpha=0.5,
    position = position_dodge(0.7))+
  scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
  scale_color_manual(values=c("#00AFBB", "#E7B800"))+
  ylab("Observed correlation (r)") +
  xlab("Remnant-to-Restored Forest Proximity") +
  theme(axis.text.x=element_text(colour="black",size=7), 
        axis.text.y=element_text(colour="black",size=7),
        axis.title=element_text(colour="black",size=9),
        legend.key=element_blank())+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background=element_blank())

mantel.proximity
dev.copy(pdf, "figures/mantel.proximity.pdf", width=5, height=6)
dev.off()


######################
######################
################# Transect 8
# Just AKs 
# plot distance: W --> E (AK1, AK2, AK3, RO4, RO5, RO6)


# T8 dat = transposed OTUs
# sample.T8 = samples
# T8.df = dataframe 
# T8_meta = the meta data
# T8_otu = otus

## T8 single plot dataframes
T8.AK1.df<-T8.df[(T8.df$Plot=="AK1"),]; T8.AK1.meta<- T8.AK1.df[,c(1:34)]; T8.AK1.otu<-T8.AK1.df[,c(35:1470)]
T8.AK2.df<-T8.df[(T8.df$Plot=="AK2"),]; T8.AK2.meta<- T8.AK2.df[,c(1:34)]; T8.AK2.otu<-T8.AK2.df[,c(35:1470)]
T8.AK3.df<-T8.df[(T8.df$Plot=="AK3"),]; T8.AK3.meta<- T8.AK3.df[,c(1:34)]; T8.AK3.otu<-T8.AK3.df[,c(35:1470)]

# combine for comparisons
T8.AK1.AK2.df<-rbind(T8.AK1.df,T8.AK2.df)
T8.AK1.AK2.meta <-rbind(T8.AK1.meta,T8.AK2.meta)
T8.AK1.AK2.otu<-rbind(T8.AK1.otu,T8.AK2.otu)

T8.AK1.AK3.df<-rbind(T8.AK1.df,T8.AK3.df)
T8.AK1.AK3.meta <-rbind(T8.AK1.meta,T8.AK3.meta)
T8.AK1.AK3.otu<-rbind(T8.AK1.otu,T8.AK3.otu)

T8.AK2.AK3.df<-rbind(T8.AK2.df,T8.AK3.df)
T8.AK2.AK3.meta <-rbind(T8.AK2.meta,T8.AK3.meta)
T8.AK2.AK3.otu<-rbind(T8.AK2.otu,T8.AK3.otu)


# Distance to neighboring plots
######### 
# AK1 - AK2
g_dists_km.T8.AK1.AK2<- earth.dist(T8.AK1.AK2.meta[,33:34],dist=TRUE)
g_dists_m.T8.AK1.AK2<- (g_dists_km.T8.AK1.AK2 * 1000)
g_dists_m.T8.AK1.AK2 <- as.matrix(g_dists_m.T8.AK1.AK2)
rownames(g_dists_m.T8.AK1.AK2) <- rownames(T8.AK1.AK2.meta)
colnames(g_dists_m.T8.AK1.AK2)<- t(rownames(T8.AK1.AK2.meta))
# rename and make distance matrix
geomat.T8.AK1.AK2<-g_dists_m.T8.AK1.AK2
g_dists_m.T8.AK1.AK2<-as.dist(g_dists_m.T8.AK1.AK2) 


#### AK1-AK3
g_dists_km.T8.AK1.AK3<- earth.dist(T8.AK1.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.AK1.AK3<- (g_dists_km.T8.AK1.AK3 * 1000)
g_dists_m.T8.AK1.AK3 <- as.matrix(g_dists_m.T8.AK1.AK3)
rownames(g_dists_m.T8.AK1.AK3) <- rownames(T8.AK1.AK3.meta)
colnames(g_dists_m.T8.AK1.AK3)<- t(rownames(T8.AK1.AK3.meta))
# rename and make distance matrix
geomat.T8.AK1.AK3<-g_dists_m.T8.AK1.AK3
g_dists_m.T8.AK1.AK3<-as.dist(g_dists_m.T8.AK1.AK3) 


#### AK2-AK3
g_dists_km.T8.AK2.AK3<- earth.dist(T8.AK2.AK3.meta[,33:34],dist=TRUE)
g_dists_m.T8.AK2.AK3<- (g_dists_km.T8.AK2.AK3 * 1000)
g_dists_m.T8.AK2.AK3 <- as.matrix(g_dists_m.T8.AK2.AK3)
rownames(g_dists_m.T8.AK2.AK3) <- rownames(T8.AK2.AK3.meta)
colnames(g_dists_m.T8.AK2.AK3)<- t(rownames(T8.AK2.AK3.meta))
# rename and make distance matrix
geomat.T8.AK2.AK3<-g_dists_m.T8.AK2.AK3
g_dists_m.T8.AK2.AK3<-as.dist(g_dists_m.T8.AK2.AK3) 


############################
#Calculate Bray-Curtis distances on T8 between RO to AK

#### AK1-2-3
T8_AK1.AK2.bc_dist = as.dist((vegdist(T8.AK1.AK2.otu, "bray")))
T8_AK1.AK2.bc_mat<-as.matrix(vegdist(T8.AK1.AK2.otu, "bray"))

T8_AK1.AK3.bc_dist = as.dist((vegdist(T8.AK1.AK3.otu, "bray")))
T8_AK1.AK3.bc_mat<-as.matrix(vegdist(T8.AK1.AK3.otu, "bray"))

T8_AK2.AK3.bc_dist = as.dist((vegdist(T8.AK2.AK3.otu, "bray")))
T8_AK2.AK3.bc_mat<-as.matrix(vegdist(T8.AK2.AK3.otu, "bray"))


#################
# Mantel Test
T8.AK1.AK2_mantel_test<-mantel.rtest(g_dists_m.T8.AK1.AK2, T8_AK1.AK2.bc_dist, nrepet=9999)
T8.AK1.AK2_mantel_test

T8.AK1.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.AK1.AK3, T8_AK1.AK3.bc_dist, nrepet=9999)
T8.AK1.AK3_mantel_test

T8.AK2.AK3_mantel_test<-mantel.rtest(g_dists_m.T8.AK2.AK3, T8_AK2.AK3.bc_dist, nrepet=9999)
T8.AK2.AK3_mantel_test


######################
######################
######################
################# Transect 9
# Just AKs 
# plot distance: W --> E (AK4, AK5, AK6, RO4, RO5, RO6)


# T9 dat = transposed OTUs
# sample.T9 = samples
# T9.df = dataframe 
# T9_meta = the meta data
# T9_otu = otus

## T9 single plot dataframes
T9.AK4.df<-T9.df[(T9.df$Plot=="AK4"),]; T9.AK4.meta<- T9.AK4.df[,c(1:34)]; T9.AK4.otu<-T9.AK4.df[,c(35:1470)]
T9.AK5.df<-T9.df[(T9.df$Plot=="AK5"),]; T9.AK5.meta<- T9.AK5.df[,c(1:34)]; T9.AK5.otu<-T9.AK5.df[,c(35:1470)]
T9.AK6.df<-T9.df[(T9.df$Plot=="AK6"),]; T9.AK6.meta<- T9.AK6.df[,c(1:34)]; T9.AK6.otu<-T9.AK6.df[,c(35:1470)]

# combine for comparisons
T9.AK4.AK5.df<-rbind(T9.AK4.df,T9.AK5.df)
T9.AK4.AK5.meta <-rbind(T9.AK4.meta,T9.AK5.meta)
T9.AK4.AK5.otu<-rbind(T9.AK4.otu,T9.AK5.otu)

T9.AK4.AK6.df<-rbind(T9.AK4.df,T9.AK6.df)
T9.AK4.AK6.meta <-rbind(T9.AK4.meta,T9.AK6.meta)
T9.AK4.AK6.otu<-rbind(T9.AK4.otu,T9.AK6.otu)

T9.AK5.AK6.df<-rbind(T9.AK5.df,T9.AK6.df)
T9.AK5.AK6.meta <-rbind(T9.AK5.meta,T9.AK6.meta)
T9.AK5.AK6.otu<-rbind(T9.AK5.otu,T9.AK6.otu)


# Distance to neighboring plots
######### 
# AK4 - AK5
g_dists_km.T9.AK4.AK5<- earth.dist(T9.AK4.AK5.meta[,33:34],dist=TRUE)
g_dists_m.T9.AK4.AK5<- (g_dists_km.T9.AK4.AK5 * 1000)
g_dists_m.T9.AK4.AK5 <- as.matrix(g_dists_m.T9.AK4.AK5)
rownames(g_dists_m.T9.AK4.AK5) <- rownames(T9.AK4.AK5.meta)
colnames(g_dists_m.T9.AK4.AK5)<- t(rownames(T9.AK4.AK5.meta))
# rename and make distance matrix
geomat.T9.AK4.AK5<-g_dists_m.T9.AK4.AK5
g_dists_m.T9.AK4.AK5<-as.dist(g_dists_m.T9.AK4.AK5) 


#### AK4-AK6
g_dists_km.T9.AK4.AK6<- earth.dist(T9.AK4.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.AK4.AK6<- (g_dists_km.T9.AK4.AK6 * 1000)
g_dists_m.T9.AK4.AK6 <- as.matrix(g_dists_m.T9.AK4.AK6)
rownames(g_dists_m.T9.AK4.AK6) <- rownames(T9.AK4.AK6.meta)
colnames(g_dists_m.T9.AK4.AK6)<- t(rownames(T9.AK4.AK6.meta))
# rename and make distance matrix
geomat.T9.AK4.AK6<-g_dists_m.T9.AK4.AK6
g_dists_m.T9.AK4.AK6<-as.dist(g_dists_m.T9.AK4.AK6) 


#### AK5-AK6
g_dists_km.T9.AK5.AK6<- earth.dist(T9.AK5.AK6.meta[,33:34],dist=TRUE)
g_dists_m.T9.AK5.AK6<- (g_dists_km.T9.AK5.AK6 * 1000)
g_dists_m.T9.AK5.AK6 <- as.matrix(g_dists_m.T9.AK5.AK6)
rownames(g_dists_m.T9.AK5.AK6) <- rownames(T9.AK5.AK6.meta)
colnames(g_dists_m.T9.AK5.AK6)<- t(rownames(T9.AK5.AK6.meta))
# rename and make distance matrix
geomat.T9.AK5.AK6<-g_dists_m.T9.AK5.AK6
g_dists_m.T9.AK5.AK6<-as.dist(g_dists_m.T9.AK5.AK6) 


############################
#Calculate Bray-Curtis distances on T9 between RO to AK

#### AK4-5-6
T9_AK4.AK5.bc_dist = as.dist((vegdist(T9.AK4.AK5.otu, "bray")))
T9_AK4.AK5.bc_mat<-as.matrix(vegdist(T9.AK4.AK5.otu, "bray"))

T9_AK4.AK6.bc_dist = as.dist((vegdist(T9.AK4.AK6.otu, "bray")))
T9_AK4.AK6.bc_mat<-as.matrix(vegdist(T9.AK4.AK6.otu, "bray"))

T9_AK5.AK6.bc_dist = as.dist((vegdist(T9.AK5.AK6.otu, "bray")))
T9_AK5.AK6.bc_mat<-as.matrix(vegdist(T9.AK5.AK6.otu, "bray"))


#################
# Mantel Test
T9.AK4.AK5_mantel_test<-mantel.rtest(g_dists_m.T9.AK4.AK5, T9_AK4.AK5.bc_dist, nrepet=9999)
T9.AK4.AK5_mantel_test

T9.AK4.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.AK4.AK6, T9_AK4.AK6.bc_dist, nrepet=9999)
T9.AK4.AK6_mantel_test

T9.AK5.AK6_mantel_test<-mantel.rtest(g_dists_m.T9.AK5.AK6, T9_AK5.AK6.bc_dist, nrepet=9999)
T9.AK5.AK6_mantel_test





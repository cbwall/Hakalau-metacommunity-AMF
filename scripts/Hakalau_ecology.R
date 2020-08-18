# Purpose: 4 plots to answer reviewers from Molecular Ecology
# 4/10/20
library(plyr)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(plotrix)

# read in data file
## ***** !NOTE! ***** data here is qualitative for rapid assessments and should not be viewed as standardized, rigourous collection of plant community percent cover data. Code below filters data and changes the % cover to presence absence

veg <- read.csv("data/ecology data/Hakalau_Veg.csv", header = T, stringsAsFactors = F, na.strings = "")

# read in corrected plot names
new_names <- read.csv("data/ecology data/NH_New_Name_WP.csv", header = T, stringsAsFactors = F)

# read in cleaned data file
trimmed_data <- read.csv("data/haka_big_transect_data_trimmed.csv",header = T)
haka_data <- trimmed_data %>% filter(NewHabType %in% c("Koa + Understory", "Remnant ohia forest", "Koa + Grass"))

# Clean up ----------------------------------------------------------------

# Fix column names of vegetation data
 veg[1,]  -> colnames(veg)
 veg[2:nrow(veg),]-> veg

#separate site metadata as EnvData
veg[1:12] %>% apply(2,as.factor) %>% as.data.frame() -> EnvData

# set rownames to site ID (WAYPOINT)
row.names(EnvData) <- EnvData$WAYPOINT

# separate out the count data and percent data from vegetation survey
  veg_column_names <- colnames(veg)
  ## count and precent separate out as odds and even columns
  count_cols <- veg_column_names[grepl("_#", veg_column_names)]
  perc_cols   <- veg_column_names[grepl("_%", veg_column_names)]
  
  ## clean count data
  veg[count_cols]-> count_data
  count_data[is.na(count_data)] <-0
  count_data[count_data=="x"]   <- 1
  count_data <- apply(count_data, 2, as.integer) %>% as.data.frame()
  
  row.names(count_data) <- veg$WAYPOINT
  
  ## clean percent data
  veg[perc_cols] -> perc_data
  perc_data[is.na(perc_data)] <-0
  perc_data <- apply(perc_data, 2, as.integer) %>% as.data.frame()
  
  ## combine count and percent data as presence absence
  pa_data <- Reduce("+", list(count_data,perc_data))
  pa_data[pa_data > 0] <- 1
  
  pa_data <- pa_data[colSums(pa_data)>0]
  colnames(pa_data) <- sub("_#","",colnames(pa_data))
  pa_data$WAYPOINT <- veg$WAYPOINT
  pa_data$FOCAL <- veg$FOCALGROUP
  
  ## subset presence absence data to the same plots as the other graphs (haka data)
  pa_data <- pa_data[pa_data$WAYPOINT %in% haka_data$WAYPOINT,]
  
  # remove a couple of understory species
  pa_data <- pa_data %>% select(-PASMOL, -MYOSAN)
  
  # order by Focal group and set levels for Waypoints
  pa_data <- arrange(pa_data, FOCAL)
  pa_data$WAYPOINT <- factor(pa_data$WAYPOINT, levels = pa_data$WAYPOINT)
  
  # run a clustering algorithm to arrange plants nicely
  plant.order <- colnames(pa_data)[c(hclust(dist(t(pa_data[1:14])))$order)]
  
  ## give simple names to each site
  pa_data <- left_join(pa_data, new_names[2:3], by = c("WAYPOINT" = "Name_OLD"))
 
  n_AK <- length(pa_data$FOCAL[pa_data$FOCAL == "AK"])
  n_RO <- length(pa_data$FOCAL[pa_data$FOCAL == "RO"])
  
  pa_data <- pa_data %>% mutate(clean_name = ifelse(is.na(Name_New),paste(FOCAL, "Survey", c(1:n_AK,1:n_RO) ), Name_New))

  pa_data <- pa_data %>% mutate(bold = ifelse(WAYPOINT %in% new_names$Name_OLD,"bold","plain") )
  

# heatmap of fleshy fruit for all sites -----------------------------------
# all surveyed sites with expanded community
# gather into long format for ggplot
pa_data_long <- pa_data %>% gather(key = Host_Abb, value = pres_ab, 1:14)
  
# Add a column for the different fills in the heatmap
pa_data_long <- pa_data_long %>% mutate(PresAbs_levels =
                    ifelse(pa_data_long$FOCAL == "AK" & pa_data_long$pres_ab == 1, "Present (AK)",
                    ifelse(pa_data_long$FOCAL == "RO" & pa_data_long$pres_ab == 1, "Present (RO)", "Absent")))

  # reorder to match clustering algorithm
pa_data_long$Host_Abb <- factor(pa_data_long$Host_Abb, levels = plant.order)

# filter out non-mycorrhizal plants
pa_data_long <- pa_data_long %>% filter(!(Host_Abb %in% c("LEPTAM","VACSPP","FERN","TFERN")))
pa_data_long<-droplevels(pa_data_long)

# generate a vector of bold and plain names for the heatmap
pa_data$clean_name %>% sort() -> y_axis_names
y_axis_names <- gsub(".. .*","plain",y_axis_names)
y_axis_names <- gsub("[RA].*","bold",y_axis_names)

# write out table of plant names and abbreviations
plant_key<- data.frame(plant_species = c("Cheirodendron trigynum",
                             "Coprosma spp.",
                             "Ilex anomala",
                             "Myrsine lessertiana",
                             "Phyllostegia brevidens",
                             "Clermontia lindsayana",
                             "Rubus argutus",
                             "Rubus hawaiensis",
                             "Metrosideros polymorpha",
                             "Acacia koa"
                             ), plant_abbreviation = colnames(pa_data)[c(1:3,5:9,11:12)])

write.csv(plant_key, "output/fleshy_fruited_plant_names.csv")

######################
## Make Plot!

# re-order species alphabetically
pa_data_long$Host_Abb<-factor(pa_data_long$Host_Abb, levels=c("AKOA", "CHETRI", "CLELIN", "COPSPP", "ILEANO", 
                                                              "MP", "MYRLES", "PHYBRE", "RUBARG", "RUBHAW"))


pa_data_long$Host_Abb.names<-pa_data_long$Host_Abb
pa_data_long$Host_Abb.names<- revalue(pa_data_long$Host_Abb.names,
                                      c("AKOA"="Acacia koa", 
                                        "CHETRI"="Cheirodendron trigynum",
                                        "CLELIN"="Clermontia lindsayana",
                                        "COPSPP"="Coprosma spp.",
                                        "ILEANO"="Ilex anomala  ",
                                        "MP"="Metrosideros polymorpha ",
                                        "MYRLES"="Myrsine lessertiana ",
                                        "PHYBRE"="Phyllostegia brevidens ",
                                        "RUBARG"="Rubus argutus",
                                        "RUBHAW"="Rubus hawaiensis"))


pa_data_long$clean_name<-as.factor(pa_data_long$clean_name)

#reorder
pa_data_long$clean_name<-factor(pa_data_long$clean_name, 
      levels=c("AK1", "AK2", "AK3", "AK4","AK5", "AK6", 
               "AK Survey 1", "AK Survey 2", "AK Survey 3", "AK Survey 4", "AK Survey 5",  
               "AK Survey 6", "AK Survey 7", "AK Survey 8", "AK Survey 9", "AK Survey 10",
               "AK Survey 11", "AK Survey 12", "AK Survey 14", "AK Survey 15",
               "AK Survey 16", "AK Survey 17", "AK Survey 18", "AK Survey 19", "AK Survey 20",
               "AK Survey 23", "AK Survey 24", "AK Survey 26", "AK Survey 27", "AK Survey 30", 
               "RO1", "RO2", "RO3",  "RO4", "RO5", "RO6", 
               "RO Survey 1", "RO Survey 5", "RO Survey 6", "RO Survey 8", "RO Survey 9", "RO Survey 10", 
               "RO Survey 12", "RO Survey 13", "RO Survey 14", "RO Survey 15", "RO Survey 16", 
               "RO Survey 17", "RO Survey 18", "RO Survey 20", "RO Survey 21", "RO Survey 22",
               "RO Survey 23", "RO Survey 24", "RO Survey 25"))
        
p1 <- ggplot(pa_data_long, aes(Host_Abb.names, clean_name)) +
        geom_tile(aes(fill = PresAbs_levels)) +
        scale_fill_manual(values = c("white", "gray", "black"))+
        labs(x= "Plant Species", y= "Plot", fill = "Presence In Plot") +
  theme_classic() +
        theme(axis.text.x = element_text(face="italic", angle = 45, hjust = 1),
        axis.text.y = element_text(size=6))
p1
ggsave("figures/Fig S2.surveyed_heatmap.pdf",width = 6, height = 6)                    




# Qualitative DBH % ----------------------------------------------------------

dbh_data <- haka_data %>% group_by(FOCALGROUP, DBH) %>% summarise(count= n())
colnames(dbh_data)[1] <- "Focal.Group"
dbh_data$DBH<-factor(dbh_data$DBH, levels=c("S", "M", "L"))
dbh_data$Focal.Group<-factor(dbh_data$Focal.Group, levels=c("RO", "AK"))

p2 <-ggplot(dbh_data, aes(fill= Focal.Group, x = DBH, y = count)) +
      geom_bar(position = "dodge", stat = "identity", color = "black") +
      scale_fill_manual(values = c("white", "black")) +
      labs(x = "Tree dbh (qualitative score)", y = "Plots Surveyed") + theme_classic()

p2
ggsave(filename = "figures/Fig. Sxx.Qualitative_dbh.pdf", width = 5, height = 5)


  
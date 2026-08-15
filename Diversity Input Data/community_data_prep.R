#-------------------------------------------------------------------------------
#
#  Community Data preparation  
#                 
#-------------------------------------------------------------------------------

# Author

# Created

# Descritption:
# This script prepares the stream fish community data used to calculate multi-
# dimensional beta diversity. This script first compiles community data from 3
# datasets (Appendix 1 for sources), then filters based on criteria in Appendix  
# 2, and rarifies data to ensure even sampling coverage based on number of  
# stream segments sampled. Finally, data are split into two sampling pools, 
# contemporary, which includes all records, and native-only, which serves as a 
# loose proxy for historic data.


# House Keeping ----------------------------------------------------------------
rm(list = ls())

# Load Libraries
library(dplyr)
library(tidyr)
library(stringr)

# Load community data

all_com <- readRDS()  # Data from state agencies listed in Appendix 1

# Load supplemental community data
giam_olden <- readRDS()  # Data from Giam & Olden (2016)
olden <- readRDS()  # Data from Pool et al. (2010) 

# add state column so data sets connect.
giam_olden$State <- NA
olden$State <- NA

giam_olden <- giam_olden[giam_olden$FishScales == F,]
olden <- olden[olden$GiamOlden16 == F,]

giam_olden <- giam_olden %>% select(-FishScales,-LCRBolden)
olden <- olden %>% select(-FishScales,-GiamOlden16)

# Compile all community data
all_com <- rbind(all_com,olden,giam_olden)


# Select only species level records of stream fish  ----------------------------

### Genus level ids ###

# Dataset includes ids at subspecies, species, provisional species, and genus
# level ids. All of these are fine but genus level

# choose only species and subspecies level ids.
com <- all_com[!all_com$Lowest_ID %in% c("genus"),]   

### Subspecies ###

# These data include some subspecies level ids, these are finer than we need, so
# subspecies will be coarsened to species level

# Change common name to that of parents species
parent_common <- com %>% 
  distinct(Scientific_Name, Common_Name)

# create table that has sci name of subspecies and common name of parent species
sub_species <- com %>% 
  filter(Lowest_ID == "subspecies") %>%   
  distinct(Scientific_Name) %>%   
  mutate(Parent_Name = word(Scientific_Name,1,2)) %>%   
  left_join(parent_common, by =c("Parent_Name" = "Scientific_Name"))  

# Change comman names of subspeces
com <- com %>% 
  left_join(sub_species, by = "Scientific_Name") %>%  
  mutate(
    Common_Name = ifelse(
      Lowest_ID == "subspecies",  
      Common_Name.y,  
      Common_Name.x
      )
    ) %>% 
  select(-Parent_Name,-Common_Name.x,-Common_Name.y) 

# Coarsen sci name of subspecies
com$Scientific_Name[com$Lowest_ID == "subspecies"] <- 
  stringr::word(com$Scientific_Name[com$Lowest_ID == "subspecies"], 1,2)

# Coarsen native status of subspecies
sub_huc_list <- readRDS()  # read in list of all hucs of parent and  subspecies
sub_parent <- names(sub_huc_list)  # create vector of parent names

# Check
unique(com$Native8[
  com$Scientific_Name %in% sub_parent[6] & com$HUC_8 %in% sub_huc_list[[6]]
  ])

# For loop to coarsen native status
for ( i in 1:length(sub_parent)) {
  com$Native8[
    com$Scientific_Name %in% sub_parent[i] & com$HUC_8 %in% sub_huc_list[[i]]
    ] <- T  
}

# Check
unique(com$Native8[
  com$Scientific_Name %in% sub_parent[1] & com$HUC_8 %in% sub_huc_list[[1]]
  ])

### Non freshwater fish ###

# samples include records of sharks and rays, these are outside the scope of
# our assemblages and will be removed

# Remove sharks
com <- com[
  !com$Order %in% c("Myliobatiformes","Carcharhiniformes") |
    is.na(com$Order),
  ]

# Remove brackish species, these are genera not found in freshwater trait 
# data base, verified to be mostly marine by literature review. See Marine  
# Genera identification.Rfor more details
marine_genera <- c(
  "Achirus","Anchoa","Agonostomus","Archosargus","Ariopsis",
  "Awaous","Bagre","Bairdiella","Bothus","Brevoortia","Caranx",
  "Centropomus","Citharichthys", "Ctenogobius","Cymatogaster",
  "Cynoscion","Dormitator", "Eleotris","Elops","Entosphenus",
  "Eucinostomus", "Eugerres","Evorthodus","Gobioides", 
  "Gobiomorus","Gobiosoma","Lagodon","Leiostomus","Lutjanus",
  "Megalops","Membras","Microgadus","Microgobius",
  "Micropogonias","Mugil","Myrophis", "Oligolepis", 
  "Oligoplites", "Ophichthus","Opsanus","Paralichthys",
  "Pogonias","Pomatomus","Proterorhinus","Sciaenops",
  "Strongylura","Syngnathus","Trinectes"
  )

com <- com[!com$Genus %in% marine_genera,]

### Export species names for phylo tree ###

# # Export list of species names. will use this to build phylo tree
# all_names <- as.data.frame(unique(com$Scientific_Name))
# write.csv(all_names, "full_species_list.csv",row.names = F)


# Filter community data  -------------------------------------------------------

### Site characteristics ###

# To ensure comparability of samples, we want recent (2000s and up), 
# electrofishing surveys of wadeable streams (watershed area < 300 and water 
# type = stream river)
com_filter <- com %>% 
  filter(
    Electrofishing == T,
    Year > 2000,
    WatershedArea < 300,
    WaterType == "StreamRiver"
    )

### No native status ###

# some species do not have native range status. To be conservative any site with
# one of these species will be remove to ensure we have native status for whole
# community

# find comids with at least one species without native information
no_native_comid <- unique(
  com_filter$COMID[com_filter$NativeStatusAvailable == F]
  )  

# remove sites without native info for at least one fish

com_filter <- com_filter[!com_filter$COMID %in% no_native_comid,]  


### Choose one sampling site per comid ###

# Most species rich survey
set.seed(123)
com_mrich_rd <- com_filter %>% 
  group_by(COMID,Date,X,Y) %>%   
  summarise(richness = n_distinct(Scientific_Name)) %>%   # number of species
  group_by(COMID) %>%   
  filter(richness == max(richness)) %>%   # most sp rich survey
  slice_sample(n = 1) %>%   # if 2 surveys have  equal richness, choose 1 random
  select(COMID, Date, X, Y) %>%  
  left_join(com_filter, by = c("COMID", "Date", "X", "Y"))


### Rarefy data ####
# beta diversity is sensitive to sample density, so ensure all hucs have an 
# equivalent sampling density

## first, remove all sites with only 1 native species. 2 species are needed per
# site to calculate kernel based functional beta diversity. Native species make
# up historic community, so if two species exist there, at least two will exist 
# in contemporary community
com_mrich_rd_n1rm <- com_mrich_rd %>% ungroup() %>% 
  filter(Native8 == T) %>% 
  group_by(COMID,Date, X, Y) %>% 
  summarise(richness = n_distinct(Scientific_Name))  %>%   
  filter(richness > 1) %>% 
  select(COMID, Date, X, Y) %>%  
  left_join(com_filter, by = c("COMID", "Date", "X", "Y"))  


### HUC 2  rarefication ###
## Load in huc comid join
huc_comid <- readRDS() # spatially joined HUC2 and COMID shapefile
huc_comid <- sf::st_drop_geometry(huc_comid)

huc2_comids <- huc_comid %>%
  group_by(HUC_2) %>%
  summarise(tot_comid = n_distinct(COMID))


# Find number of COMID/HUC2 per densities
site_density <- com_mrich_rd_n1rm %>% 
  mutate(HUC_2 = substr(HUC_12,1,2)) %>% 
  group_by(HUC_2,COMID) %>% 
  summarise(n_comid = n_distinct(COMID)) %>% 
  group_by(HUC_2) %>% 
  summarise(n_comid = n()) %>% 
  left_join(huc2_comids) %>% 
  mutate(site_density_nhd = n_comid/tot_comid)
site_density <- site_density[order(site_density$site_density_nhd),] 
huc2_site_n <- site_density %>% select(HUC_2,site_density_nhd)

# Calculate the number of sites and HUC2s of final rarified data using each 
# HUC2's sample density
for (i in 1:nrow(huc2_site_n)) {
  data <- site_density %>% 
    mutate(rare_sites = site_density_nhd[i]*tot_comid)
  data$final_sites <- ifelse(
    round(data$rare_sites) > data$n_comid,
    NA,
    round(data$rare_sites)
    )
  huc2_site_n[i,"total_sites"] <- sum(data$final_sites,na.rm = T)
  huc2_id <- data$HUC_2[!is.na(data$final_sites)]
  huc2_id <- huc2_id[order(huc2_id)]
  huc2_site_n[i,"total_huc2"] <- length(huc2_id)
  huc2_site_n[i,"huc2_id"] <- paste(huc2_id,collapse = ", ")
}
# Sampling denisty of HUC2 17 (rio grande) balances extent (only excludes 1 
# HUC2) with sample size (1,023 sites)

# Calculate rarified sample densities for each huc2 based on Rio Grande 
# sample density (0.0003943076)
site_density <- com_mrich_rd_n1rm %>% 
  mutate(HUC_2 = substr(HUC_12,1,2)) %>% 
  group_by(HUC_2,COMID) %>% 
  summarise(n_comid = n_distinct(COMID)) %>% 
  group_by(HUC_2) %>% 
  summarise(n_comid = n()) %>% 
  left_join(huc2_comids) %>% 
  mutate(site_density_nhd = n_comid/tot_comid, 
         rare_sites = 0.0003943076*tot_comid) 

# Round sites and exclude HUC2s with insufficient samples
site_density$final_sites <- ifelse(
  round(site_density$rare_sites) > site_density$n_comid,
  NA,
  round(site_density$rare_sites)
  )

# Total number of sites
sum(site_density$final_sites,na.rm = T)


## Rarify data by huc2 ##

# Create dataframe with final site counts for each HUC1
huc2_site_counts <- site_density %>% 
  filter(!is.na(final_sites)) %>% 
  select(HUC_2,final_sites)

## For loop to select random comids in each HUC2 based on number sites
sample_list2 <- list() # list to contain number of sites for each huc2

for (i in 1:nrow(huc2_site_counts)) {
  comids <- unique(
    com_mrich_rd_n1rm$COMID[
      substr(com_mrich_rd_n1rm$HUC_12,1,2) == huc2_site_counts$HUC_2[i]
      ]
    )  # pull out all comids from each huc
  
  # random sample sites based on sites calculated based on density, add to list
  set.seed(123)
  sample_list2[[i]] <- sample(comids,huc2_site_counts$final_sites[i])  
}
rare_samples2 <- unlist(sample_list2)  # change list to vector

# only retain rarified sites
com_rare2 <- com_mrich_rd_n1rm[com_mrich_rd_n1rm$COMID %in% rare_samples2,]  

# # Export species names form filtered and rarified data for tree building and
# # trait imputation
# filtered_species <- unique(com_rare2$Scientific_Name)
# write.csv(filtered_species,"filtered_species_list.csv")

# Change to community-species matrix -------------------------------------------

# This creates the contemporary and native only data sets. contemporary data 
# set includes all native and nonnative species. Native data set include all 
# species native at HUC8 level. These data can be used to compare to the 
# contemporary data set. This is not a true historic data set. This only tells 
# what stream segment was like before introductions. It gives no information on 
# species extirpation or translocations of species within HUC8 watersheds. This 
# script will use a for loop to create both data sets int the same manner. 

# 1rst item tell it to include all specie, 2nd to include just native species
TF_list <- list(c(T,F,NA),c(T))  

# create a list to contain the contemporary and native data frames
data <- rep(list(NA),2)  
file_name <- c("mod","his")  

# For loop to create modern and historic data set
for (i in 1:2) {
  com_final <- com_rare2[com_rare2$Native8 %in% TF_list[[i]],]  
  
  # Create community species matrices
  data[[i]] <- com_final %>% 
    group_by(COMID,HUC_12, Scientific_Name) %>%  
    summarise(presence = n_distinct(COMID)) %>%  
    # Change into site species matrix
    pivot_wider(
      names_from = Scientific_Name, 
      values_from = presence, 
      values_fill = 0
      )  
  
  # Format data frame
  data[[i]] <- as.data.frame(data[[i]])  
  rownames(data[[i]]) <- data[[i]]$COMID  
  data[[i]] <- data[[i]][,-1]  
}

# Keep only COMIDS that appear in both historic and modern data
mod_com <- data[[1]][rownames(data[[1]]) %in% rownames(data[[2]]),]
his_com <-data[[2]][rownames(data[[2]]) %in% rownames(data[[1]]),]


# Check for sites that have no species. This should be removed by only 
# including sites present in both historic and modern data set but this code 
# will identify if some coding error resulted a site having species remove
nrow(mod_com[rowSums(mod_com[,-1]) == 0,-1])
nrow(his_com[rowSums(his_com[,-1]) == 0, -1])
# no rows are empty


# Export data  -----------------------------------------------------------------
saveRDS(mod_com,"mod_com_diversity_input.rds")  # save to project folder folder
saveRDS(his_com,"his_com_diversity_input.rds")   # save to project folder folder
saveRDS(com_rare2,"raw_community_diversity_input.rds")



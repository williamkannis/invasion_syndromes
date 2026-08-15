#-------------------------------------------------------------------------------
#
#  Trait data preparation  
#             
#-------------------------------------------------------------------------------

# Author:
  
# Created:

# Description:
# This script has code to combine Giam & Olden (2016)'s dataset with 
# Frimpong & Angermeier (2009)'s dataset and to estimate missing trait values 
# using ancestral character state with the phyloestimate functions. These are 
# estimated separately for continuous and discrete and merged together. Dummy
# variables were created for discrete trait data and final file exported to be 
# used to calculate functional beta diversity in beta_diversity.R


# House Keeping  ---------------------------------------------------------------
rm(list=ls())

# Load Packages
library(dplyr)
library(tidyverse)
library(picante)

# Species names for those included in dataset and analysis
species <- read.csv("full_species_list.csv")
filter_species <- read.csv("filtered_species_list.csv")[,1]

# Reformat data sets -----------------------------------------------------------


# In the section below, data sets are reformatted so they can be merged. In
# addition, species names are standardized so that they match the community and
#  phylogenetic datasets. Finally, any errors in the datasets are corrected.

# Frimpong & Angermeier (2009)'s formatting ------------------------------------

# Load trait data
trait <- readxl::read_excel("")  # available from 
# https://www.sciencebase.gov/catalog/item/5a7c6e8ce4b00f54eb2318c0 

# Create column for species name
trait$Sci_name <- paste(trait$GENUS,trait$SPECIES, sep = " ")  
trait$Sci_name <- ifelse(
  trait$SPECIES == ".", 
  paste(trait$GENUS, trait$COMMONNAME,trait$OTHERNAMES, sep = " "),
  paste(trait$GENUS,trait$SPECIES, sep = " ")
  )

# Some genera have changed name since dataset was created, change names to  
# match database

# Phoxinus -> Chrosomus
trait$Sci_name <- ifelse(
  stringr::word(trait$Sci_name,1) =="Phoxinus", 
  paste("Chrosomus ", stringr::word(trait$Sci_name,2),sep = ""), 
  trait$Sci_name
  )

# Cichlasoma cyanoguttatum -> Herichthys cyanoguttatus
trait$Sci_name[trait$Sci_name == "Cichlasoma cyanoguttatum"] <- 
  "Herichthys cyanoguttatus"

# Cichlasoma octofasciatum -> Rocio octofasciata"
trait$Sci_name[trait$Sci_name == "Cichlasoma octofasciatum" ] <- 
  "Rocio octofasciata"

# Gila bicolor to Siphateles bicolor
trait$Sci_name[trait$Sci_name == "Gila bicolor"] <- 
  "Siphateles bicolor"

### Reformat names ###

# # Format sci name to match community dataset
# fish_revised_name <- taxize::gnr_resolve(
#   trait$Sci_name, 
#   data_source_ids = 3, 
#   canonical =  TRUE
#   )
# fish_revised_name$exact_match <-ifelse(
#   fish_revised_name$user_supplied_name == fish_revised_name$matched_name2, 
#   T, 
#   F
#   )
# # 4 species needed name changed to match that of dataset, many were 
# # provisional

# Fix incorrect names
trait$Sci_name[trait$Sci_name == "Dorosoma petenese"] <- 
  "Dorosoma petenense"
trait$Sci_name[trait$Sci_name == "Cyprinidae . sawfin shiner"] <- 
  "Notropis sp. Sawfin Shiner"
trait$Sci_name[trait$Sci_name == "Percidae . sunburst darter"] <- 
  "Etheostoma mihileze"
trait$Sci_name[trait$Sci_name == "Percidae . Muscadine darter"] <- 
  "Percina smithvanizi"

# # Check to see if dataset names are up to date. If not, change trait names to
# # match older dataset names
# tree_revised_name <- taxize::gnr_resolve(
#   species, 
#   data_source_ids = 3, 
#   canonical =  TRUE
#   )
# tree_revised_name$exact_match <-ifelse(
#   tree_revised_name$user_supplied_name == tree_revised_name$matched_name2, 
#   T, 
#   F
#   )
# # Two species were spelled correctly in trait but incorrectly in dataset, 
# # change trait names to match incorrect names

# This species is spelled this way in dataset
trait$Sci_name[trait$Sci_name == "Gymnocephalus cernuus"] <- 
  "Gymnocephalus cernua"
trait$Sci_name[trait$Sci_name =="Lampetra appendix"] <- 
  "Lethenteron appendix"

# Format data table to match other trait databases
trait_format <- trait %>% 
  select(
    -SID,-ALTSID,-NOTES,-FID,-GENUS,-SPECIES,-GID,-ITISTSN,-COMMONNAME,
    -OTHERNAMES,-FamilyNumber
    )

# Fix negative vlaues to NA
trait_format <- as.data.frame(trait_format)
trait_format[trait_format<0] <- NA  

# Select only species in dataset
ft <- trait_format[trait_format$Sci_name %in% species,]


# Giam & Olden (2016) formatting -----------------------------------------------

# Load trait data
trait <- readxl::read_excel("") # available upon request from authors of 
# Giam & Olden (2016)

# Format species names column
trait$SCINAME <- stringr::str_to_sentence(trait$SCINAME)

# Some genera have changed name since dataset was created, change names to 
# match database

# Phoxinus -> Chrosomus
trait$SCINAME <- ifelse(
  stringr::word(trait$SCINAME,1) =="Phoxinus",
  paste("Chrosomus ", stringr::word(trait$SCINAME,2),sep = ""), 
  trait$SCINAME
  )

# Gila bicolor to Siphateles bicolor
trait$SCINAME[trait$SCINAME == "Gila bicolor"] <- 
  "Siphateles bicolor"

# Cichlasoma octofasciatum -> Rocio octofasciata"
trait$SCINAME[trait$SCINAME == "Cichlasoma octofasciatum" ] <- 
  "Rocio octofasciata"

### Reformat names ###

# # Format sci name to match community dataset
# fish_revised_name <- taxize::gnr_resolve(
#   trait$SCINAME, 
#   data_source_ids = 3, 
#   canonical =  TRUE
#   )
# fish_revised_name$exact_match <-ifelse(
#   fish_revised_name$user_supplied_name == fish_revised_name$matched_name2, 
#   T, 
#   F
#   )
# # Only two names not revised. One is a species not in taxize database, the 
# # other is spelled the same as in dataset,

# Check to see if dataset names are up to date. If not, change trait names to 
# match dataset
# tree_revised_name <- taxize::gnr_resolve(
#   species, 
#   data_source_ids = 3, 
#   canonical =  TRUE
#   )
# tree_revised_name$exact_match <-ifelse(
#   tree_revised_name$user_supplied_name == tree_revised_name$matched_name2, 
#   T, 
#   F
#   )

# This species is spelled this way in dataset
trait$SCINAME[trait$SCINAME == "Gymnocephalus cernuus"] <- 
  "Gymnocephalus cernua"
trait$SCINAME[trait$SCINAME =="Lampetra appendix"] <- 
  "Lethenteron appendix"

# Format data table to match other trait databases
trait_format <- trait %>% 
  rename("Sci_name"= "SCINAME") %>% 
  select(-FID,-FAMILY,-Abrev,-PrimarySource)

# Check for errors
summary(trait_format)
# Egg size is listed as character, check to see that all non numeric 
# values are NAs
unique(trait_format$EggSize[is.na(as.numeric(trait_format$EggSize))])
# one value is 1..715, remove extra period and change to numeric
trait_format$EggSize[
  trait_format$EggSize == "1..715" & !is.na(trait_format$EggSize)
  ] <- 1.715
# Change to numeric
trait_format$EggSize <- as.numeric(trait_format$EggSize)

# select only species in dataset
jt <- trait_format[trait_format$Sci_name %in% species,]


# Trait compilation ------------------------------------------------------------

# Merge datasets
traits <- jt %>% 
  full_join(ft)

# Select traits for beta diveristy calculation
traits_select <- traits %>% select(
  Sci_name,MAXBODYL,AGEMAT,LONG,FECUND,EggSize,
  PARENTC,TROPHICG,SUBPREF,VERTP,SLOWCURR, 
  MODCURR, FASTCURR, MAXTEMP
  )
traits_select <- as.data.frame(traits_select)
row.names(traits_select) <- traits_select$Sci_name
traits_select <- select(traits_select, -Sci_name)


# Missing trait imputation -----------------------------------------------------

# In this section, missing trait values or species are imputed based on 
# phylogenetic relationships. If trait state cannot be imputed  in this manner,
# they will be manually filled in based on literature review.

# Import Phylo Tree
tree <- readRDS("phylo_tree.rds")  # super tree created in Appendix 3
tree_fix <- ape::multi2di(tree)  # allows tree to be used for ace!

# Trim tree to include all species in trait database and to include all missing 
# species from analysis
tree_species <- tree_fix$tip.label[
  tree_fix$tip.label %in% c(row.names(traits_select),filter_species)
  ]  
tree_fix <- ape::keep.tip(tree_fix,tree_species)


# Separate data into cont and discrete traits
cont_columns <- c("MAXBODYL","AGEMAT","LONG","FECUND","EggSize","MAXTEMP")
traits_con <- select(traits_select, all_of(cont_columns))
traits_disc <- select(traits_select, -cont_columns)


# Cont trait imputation --------------------------------------------------------

# list of data frames that include the estimated and se of missing traits
phy_trait_list <- rep(list(NA), ncol(traits_con)) 

# list of data frames with estimated traits
final_trait_list <- phy_trait_list 

# Continuous traits
for (i in 1:ncol(traits_con)){
  # make list of each trait without NA values
  trait <- traits_con[!is.na(traits_con[,i]),][i] 
  
  # Estimate NA values
  phy_trait_list[[i]] <- phyEstimate(phy = tree_fix, trait = trait) 
  
  # pull out estimated values
  est_trait <- phy_trait_list[[i]]["estimate"]  
  
  # Change column names so trait values can be added to main data
  colnames(est_trait) <- colnames(trait) 
  
  # add phly estimated traits to trait data
  final_trait_list[[i]] <- rbind(trait, est_trait)  
}

# Combine traits into one data frame
merge_by_rownames <- function(df1, df2) {
  full_join(
    df1 %>% rownames_to_column('rowname'), 
    df2 %>% rownames_to_column('rowname'), 
    by = 'rowname'
    ) %>% 
    column_to_rownames('rowname')
}  # Function for combining by row name
final_traits_con <- purrr::reduce(final_trait_list,merge_by_rownames)

# remove all species not in analysis
final_traits_con_sub <- final_traits_con[
  row.names(final_traits_con) %in% filter_species,
  ]


# Discrete traits imputation function fix  -------------------------------------

# function has error in its code for choosing best state. This is not 
# a mathematical error but data handling. I went ahead and modified the source
# code to fix this. It is important to note, no changes to ancestral 
# state estimations were made

phyEstimateDisc_BA_fix <- function(
    phy, trait, best.state = TRUE, cutoff = 0.5, ...) {
  # trait should be a data.frame or vector with names matching phylogeny
  
  if (is.vector(trait) || !is.factor(trait)) {
    
    trait <- data.frame(trait)
    
    trait[, 1] <- factor(trait[, 1])
  }  ### FIXED ###
  trait.orig <- trait
  
  # given a tree with a novel taxa on it (taxa with no trait value)
  sppObs <- row.names(trait)
  sppUnobs <- phy$tip.label[!(phy$tip.label %in% sppObs)]
  trtlevels <- levels(trait[, 1])
  res <- as.data.frame(
    matrix(
      nrow = length(sppUnobs), 
      ncol = length(trtlevels), 
      dimnames = list(sppUnobs, trtlevels)
      )
    )
  
  # estimate support for different states for each novel taxon
  for (i in sppUnobs) {
    # for each novel species, prune all but measured + that species
    tree <- ape::drop.tip(phy, subset(sppUnobs, sppUnobs != i))
    
    # root the tree at the novel species (leave root as trichotomy)
    tree <- ape::root(tree, i, resolve.root = FALSE)
    
    ## Mon IDE me dit que ces variables sont jamais utilisés
    
    # record branch length leading to novel species in rerooted tree
    ## edge <- Nnode(tree) - 1 + which(tree$tip.label == i)
    ## bl <- tree$edge.length[edge]
    
    # prune novel species and match new pruned tree <-> trait data
    tree <- ape::drop.tip(tree, i)
    trait <- trait.orig[tree$tip.label,]
    
    # calculate value at root node and impute to novel species
    est <- ace(trait, tree, type = "discrete")
    val <- est$lik.anc[1,]
    
    res[i,] <- val
  }
  
  # estimate the best-supported state for each taxon
  if (best.state) {
    beststate <- as.data.frame(matrix(nrow = dim(res)[1], ncol = 2))
    colnames(beststate) <- c("estimated.state", "estimated.state.support")
    rownames(beststate) <- rownames(res)
    
    for (i in 1:dim(res)[1]) {
      #if >=cutoff % taxa have same label assign a consensus taxon to node
      
      # best <- -sort(unlist-(res[i,]))[1]  # line with bug
      best <- -sort(unlist(-res[i,]))[1]  ### FIXED LINE ###
      # print(res[i,])
      ## la formuler originale
      ##best <- -sort(-(res[i, ]))[1]
      
      if (best >= cutoff) {
        beststate[i, 1] <- names(best)
        beststate[i, 2] <- best
      }
      else
      {
        beststate[i, 1] <- NA
        beststate[i, 2] <- NA
      }
    }
  }
  
  #return the output
  if (best.state) {
    return(cbind(as.matrix(res), beststate))
  } else {
    return(as.matrix(res))
  }
}

# Discrete trait imputation  ---------------------------------------------------

# List of full estimation outputs for each trait
phy_trait_list_disc <- rep(list(NA), ncol(traits_disc))  

# list of data frame including best state of each missing trait
final_trait_list_disc <- phy_trait_list_disc  

# Discrete traits
for (i in 1:ncol(traits_disc)){ 
  # make list of each trait without NA values
  trait <- traits_disc[!is.na(traits_disc[,i]),][i]  
  
  # Estimate NA values
  phy_trait_list_disc[[i]] <- 
    phyEstimateDisc_BA_fix(phy = tree_fix, trait = trait)  
  est_trait <- phy_trait_list_disc[[i]]["estimated.state"]
  
  # Change column names so trait values can be added to main data
  colnames(est_trait) <- colnames(trait) 
  
  # add phly estimated traits to trait data
  final_trait_list_disc[[i]] <- rbind(trait, est_trait)  
}

# Combine traits into one data frame
final_traits_disc <- purrr::reduce(final_trait_list_disc,merge_by_rownames)  

# remove all species not in analysis
final_traits_disc_sub <- final_traits_disc[
  row.names(final_traits_disc) %in% filter_species,
  ]


# Literature imputation  -------------------------------------------------------

# Some values did not have a best state, these need to be filled in 
# literature review.

## Parental care ##
row.names(final_traits_disc_sub[is.na(final_traits_disc_sub$PARENTC),])
# for parental care, there is 1: "Margariscus nachtriebi"
final_traits_disc_sub$PARENTC[
  row.names(final_traits_disc_sub) == "Margariscus nachtriebi"
  ] <- 2  # was once part of Margariscus margarita https://explorer.natureserve.org/Taxon/ELEMENT_GLOBAL.2.869003/Margariscus_nachtriebi

## Substrate ##
row.names(final_traits_disc_sub[is.na(final_traits_disc_sub$SUBPREF),])
# 11 species with NA for substrate
final_traits_disc_sub$SUBPREF[
  row.names(final_traits_disc_sub) == "Dionda nigrotaeniata"
  ] <- "Vegetation"  # filament alagea https://fishbase.mnhn.fr/FieldGuide/FieldGuideSummary.php?GenusName=Dionda&SpeciesName=nigrotaeniata&pda=1&sps=
final_traits_disc_sub$SUBPREF[
  row.names(final_traits_disc_sub) == "Labidesthes sicculus"
  ] <- "Vegetation"  # prefers vegation, weedy lakes https://nas.er.usgs.gov/queries/FactSheet.aspx?SpeciesID=318
final_traits_disc_sub$SUBPREF[
  row.names(final_traits_disc_sub) == "Menidia audens"
  ] <- "Various"  # sand or gravel https://www.fishbase.se/summary/Menidia-audens.html
final_traits_disc_sub$SUBPREF[
  row.names(final_traits_disc_sub) == "Menidia beryllina"
  ] <- "Various"  # sand or gravel https://www.fishbase.se/summary/Menidia-beryllina
final_traits_disc_sub$SUBPREF[
  row.names(final_traits_disc_sub) == "Neogobius melanostomus"
  ] <- "Various" # open sanyd area, rocks, vegation https://nas.er.usgs.gov/queries/factsheet.aspx?SpeciesID=713
final_traits_disc_sub$SUBPREF[
  row.names(final_traits_disc_sub) == "Notropis amplamala"
  ] <- "Sand"  # https://fishbase.mnhn.fr/summary/Ericymba-amplamala.html

### Change factor columns to dummy variables ###
final_traits_disc_binary <- final_traits_disc_sub

# Names of columns that need to become dummy variables
dummy_col <- c("TROPHICG", "SUBPREF", "VERTP")


# For loop to add dummy variables for each level (i) of each factor (k)
for (k in dummy_col) {
  
  # pull out unique levels for each factor
  level <- unique(final_traits_disc_binary[,k])  
  for (i in level){
    # add binary column for each level of each factor
    final_traits_disc_binary[,i] <- ifelse(
      final_traits_disc_binary[,k] == i,
      1,
      0
      )  
  }
}

# Remove factor columns
final_traits_disc_binary <- final_traits_disc_binary %>% select(-dummy_col)

# Change traits to numeric
final_traits_disc_binary <- final_traits_disc_binary %>% 
  mutate_at(1:ncol(final_traits_disc_binary),as.numeric)


# Export  ----------------------------------------------------------------------

# Combine continuous and binary traits
final_traits <- merge_by_rownames(
  final_traits_con_sub,
  final_traits_disc_binary
  )

# Export Trait table
saveRDS(final_traits, "trait_diversity_input.rds")

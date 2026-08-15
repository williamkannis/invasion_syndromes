#-------------------------------------------------------------------------------
#
#  Community invadedness data prep
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 04/28/2026

# Description: Uses raw community data to calculate species richness and
# community invadedness. Exports file to summarize invadedness per syndrome.


# House keeping ----------------------------------------------------------------
rm(list = ls())

# Directories
com_dir <- "Diversity Input Data"
out_dir <- "Analysis_data"

# Load packages
library(dplyr)
library(tidyr)

# Load data
com_df <- readRDS(file.path(com_dir,"raw_community_diversity_input.rds")) %>%
  ungroup()


# Community invadedness --------------------------------------------------------

# Pivot data frame to create columns for the number of each type of species
inv_df <- com_df %>%
  group_by(COMID,Native8) %>%
  summarise(richness = n_distinct(Scientific_Name)) %>%
  pivot_wider(
    names_from = Native8,
    values_from = richness,
    values_fill = 0
  ) %>%
  rename(
    Native_rich = "TRUE",
    Nonnative_rich = "FALSE"
  ) %>%
  mutate(prop_nonnative = Nonnative_rich/(Native_rich+Nonnative_rich))


# Export -----------------------------------------------------------------------
saveRDS(inv_df,file.path(out_dir,"com_invaded.rds"))


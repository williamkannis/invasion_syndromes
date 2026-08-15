#-------------------------------------------------------------------------------
#
#  Syndrome cluster analysis
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 04/28/2026

# Description:


# House keeping ----------------------------------------------------------------
rm(list = ls())

# Directories
beta_dir <- "~/Documents/School/PhD/Analysis/beta diversity/final_outputs"
raw_dir <- "raw_data"
fig_dir <- "figure"

# Load packages
library(dplyr)
library(tibble)
library(mclust)
library(tidyr)
library(ggplot2)
library(ggforce)

# Load data
lcbd_df <-
  readRDS(file.path(beta_dir,"delta_lcbd.rds"))
com_raw_study <-
  readRDS(file.path(raw_dir,"raw_community_diversity_input.rds"))


#  Prepare data ----------------------------------------------------------------
lcbd_clst <- lcbd_df %>%
  column_to_rownames("COMID") %>%
  dplyr::select(contains("_es")) %>%
  filter(if_any(everything(),~ abs(.x) >= 1.96))


# Clustering -------------------------------------------------------------------

# Identify number of clusters via model selection
BIC <- mclustBIC(lcbd_clst,G=1:20)
plot(BIC)
summary(BIC)

## Export version of BIC Plot ##

# png(
#   "figure/bic_plot.png",
#   width = 1800,
#   height = 1200
#   )
# par(cex = 2.5)
# plot(BIC)
# dev.off()

# Perform clustering based on model selection
clst <- Mclust(lcbd_clst,x = BIC)
summary(clst, parameters = TRUE)
plot(clst, what = "classification")
abline(h = 1.96);abline(h = -1.96)
abline(v= 0)

png(
  "Figure/clst_plot.png",
  width = 1800,
  height = 1200
  )
par(cex = 20)
plot(clst, what = "classification",cex = 1.5)
dev.off()
mclust.options("classPlotColors")

# Cluster plots  ---------------------------------------------------------------

# Combinations of dimensions to plot
combo_list <- list(
  c(7,1),  # tax total - fun total
  c(7,4),  # tax total - phy total
  c(4,1)   # phy total - fun total
)

# Create plots for each combination
lapply(combo_list, function(combo){

  # Extract axis names
  mc <- clst
  x_name <- colnames(mc$data)[combo[1]]
  y_name <- colnames(mc$data)[combo[2]]

  # Create plotting data
  plot_df <- data.frame(
    x = mc$data[, x_name],
    y = mc$data[, y_name],
    classification = factor(mc$classification)
  ) %>%
    mutate(
      cluster = case_when(
        classification == 5 ~ 1,
        classification == 6 ~ 2,
        classification == 4 ~ 3,
        classification == 2 ~ 4,
        classification == 3 ~ 5,
        classification == 1 ~ 6,
        T~NA
      ),
      cluster = as.factor(cluster)
    )


  # Create plot
  p <- ggplot(plot_df,
              aes(x = x,
                  y = y,
                  colour = cluster,
                  shape = cluster)) +
    geom_point(size = 5) +
    stat_ellipse(level = 0.95, linewidth = 1.5) +
    geom_vline(
      xintercept = c(-1.96, 1.96),
      linetype = "dashed",
      linewidth = 1.5
    ) +
    geom_hline(
      yintercept = c(-1.96, 1.96),
      linetype = "dashed",
      linewidth = 1.5
    ) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
      x = "",
      y = "",
      colour = ""
    ) +
    theme_bw(base_size = 16) +
    scale_x_continuous(
      breaks = c(-1.96, 0, 1.96)
    ) +
    scale_y_continuous(
      breaks = c(-1.96, 0, 1.96)
    )+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 20),
      legend.position = "none",
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 2
      )
      # legend.title = element_text(size = 16),
      # legend.text = element_text(size = 14)
    )
  print(p)

  # Export file name
  file_name <- paste0("pair_clst_",x_name,"_",y_name,".png")

  # Export
  ggsave(
    filename = file.path(fig_dir,file_name),
    plot = p,
    width = 6,
    height = 6,
    dpi = 600
  )
}
)


# Cluster LCBD summaries ------------------------------------------------------

## Attach cluster data to raw data ##
diff_cluster_data <- data.frame(
  COMID = names(clst$classification), # create dataframe for comid and clusters
  classification = clst$classification
  ) %>%
  mutate(
    cluster = case_when(
      classification == 5 ~ 1,
      classification == 6 ~ 2,
      classification == 4 ~ 3,
      classification == 2 ~ 4,
      classification == 3 ~ 5,
      classification == 1 ~ 6,
      T~NA
    )
  )


raw_data <- lcbd_clst # Create datframe with raw delta lcbd
raw_data$COMID <- row.names(raw_data)  # create column for COMID
raw_data <- diff_cluster_data %>% left_join(raw_data)  # comine cluster and raw data

# Update cluster names

# Significant cat function
es_sig <-function(x) {
  ifelse(x >= 1.96, 3,
         ifelse(x >= 1.64,2,
                ifelse(x >= 0,1,
                       ifelse(x > -1.64,-1,
                              ifelse(x > -1.96,-2,-3)))))

}


# # Mean delta LCBD
# cluster_means_raw <- raw_data %>%
#   group_by(cluster) %>%
#   summarise(
#     across(-COMID,mean,na.rm=T),
#     n_sites = n()
#     ) %>%
#   mutate(across(everything(), es_sig))


# Median delta LCBD (used)
cluster_med_raw <- raw_data %>%
  group_by(cluster) %>%
  summarise(
    across(-COMID,median,na.rm=T),
    n_sites = n()
  ) %>%
  mutate(across(-n_sites, es_sig)) %>%
  ungroup()



# Cluster Species occurrences summary ------------------------------------------

# Format raw species data
com_raw_study$COMID <- as.character(com_raw_study$COMID)  # change COMID to character
com_raw_study <- com_raw_study %>%
  left_join(diff_cluster_data) %>%  # join species data to cluter data
  mutate(
    cluster = case_when(
      is.na(cluster) ~ 7,
      T~cluster
    )
  )
  filter(!is.na(cluster))

# Nonnaitve occurances per cluster
sp_summary_all <- com_raw_study %>%   # summarize by cluster
  filter(Native8 == F) %>%
  group_by(Scientific_Name) %>%
  summarise(n_sites = n())

sp_summary <- com_raw_study %>%   # summarize by cluster
  filter(Native8 == F) %>%
  group_by(cluster,Scientific_Name) %>%
  summarise(n_sites = n())
n_sites <- com_raw_study %>%  # add number of sites
  group_by(cluster) %>%
  summarise(tot_sites = n_distinct(COMID))
sp_summary <- sp_summary %>%   # calculate proprotion of sites occupied
  left_join(n_sites) %>%
  mutate(prop_sites = n_sites/tot_sites) %>%
  dplyr::select(-tot_sites) %>%
  arrange(cluster,desc(prop_sites))

# Native occurances per cluster
sp_summary_native_all <- com_raw_study %>% # summarize by cluster
  filter(Native8 == T) %>%
  group_by(Scientific_Name) %>%
  summarise(n_sites = n())
sp_summary_native <- com_raw_study %>% # summarize by cluster
  filter(Native8 == T) %>%
  group_by(cluster,Scientific_Name) %>%
  summarise(n_sites = n())

sp_summary_native <- sp_summary_native %>%  # calculate proprotion of sites occupied
  left_join(n_sites) %>%
  mutate(prop_sites = n_sites/tot_sites) %>%
  dplyr::select(-tot_sites) %>%
  arrange(cluster,desc(prop_sites))

# ## Export ##
write.csv(sp_summary, "figure/cluster_nonnative_species.csv")
write.csv(sp_summary_native, "figure/cluster_native_species.csv")

# Cluster species richness summary  --------------------------------------------

# MAKE THIS SCRIPT IN DIERISTY INPUT. USE IT TO MAKE FILE OF RICHNESS AND INVADEDNESS

# Calculate taxonomic richness values for each site

comid_richness <- com_raw_study %>%
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
    ) %>%  # rename columes to reflect native/nonnative richness
  mutate(prop_nonnative = Nonnative_rich/(Native_rich+Nonnative_rich)) %>%
  left_join(diff_cluster_data) %>%
  mutate(
    cluster = case_when(
      is.na(cluster) ~ 7,
      T~cluster
    )
  )
  filter(!is.na(cluster))
comid_richness$cluster <- as.factor(comid_richness$cluster)


rich_sum <- comid_richness %>%
  group_by(cluster) %>%
  summarise(mean_native_rich = mean(Native_rich),
    med_native_rich = median(Native_rich),
    mean_nonnative_rich = mean(Nonnative_rich),
    med_nonnative_rich = median(Nonnative_rich),
    mean_nonnative_prop = mean(prop_nonnative),
    med_nonnative_prop = median(prop_nonnative))


# Cluster richness plots  ------------------------------------------------------


library(ggplot2)
for(i in 1:7){
  plot_data <- comid_richness %>%
    filter(cluster == i) %>%
    rename(
      N = Native_rich,
      NN = Nonnative_rich
    ) %>%
    select(-prop_nonnative,-cluster) %>%
    pivot_longer(!COMID)

  p<-ggplot(plot_data,aes(x= name,y=value,fill=name))+
    geom_boxplot() +
    scale_fill_discrete(direction = -1)+
    ylim(0,35) +
    xlab("")+
    ylab("")+
    theme_classic(
      base_size = 50,
      base_rect_size = 2)+
    theme(
      legend.position = "none",
      axis.text.x = element_blank()
      )

  # print(p)

  plot_file <- paste0("figure/rich_plot_s",i,".png")
  ggsave(
    plot_file,
    p,
    width = 10,
    height = 10,
    dpi = 600
    )

}



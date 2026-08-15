#-------------------------------------------------------------------------------
#
#  Syndrome cluster analysis
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 04/28/2026

# Description: Performs model based clustering to identify invasion
# syndromes based on shared changes in multidimensional LCBD effect sizes.
# Filters sites so only those with at least one dimension with a significant
# change are included in cluster analyses. Creates summary tables and plots.


# House keeping ----------------------------------------------------------------
rm(list = ls())

# Directories
beta_dir <- "Analysis_data"
fig_dir <- "Figures"

# Load packages
library(dplyr)
library(tibble)
library(mclust)
library(tidyr)
library(ggplot2)
library(ggforce)

# Load data
lcbd_df <- readRDS(file.path(beta_dir,"delta_lcbd.rds"))
inv_df <- readRDS(file.path(beta_dir,"com_invaded.rds"))


#  Prepare data ----------------------------------------------------------------

# Select only sites with significant changes in at least one dimension
lcbd_clst <- lcbd_df %>%
  column_to_rownames("COMID") %>%
  filter(if_any(everything(),~ abs(.x) >= 1.96))


# Cluster model selection  -----------------------------------------------------

# Identify number of clusters via model selection
BIC <- mclustBIC(lcbd_clst,G=1:20)
plot(BIC)
summary(BIC)

# BIC Plot for Figure s6.2
png(
  file.path(fig_dir,"figure_s6.2.png"),
  width = 1800,
  height = 1200
  )
par(cex = 2.5)
plot(BIC)
dev.off()

# Clustering -------------------------------------------------------------------

# Perform clustering based on model selection
clst <- Mclust(lcbd_clst,x = BIC)
summary(clst, parameters = TRUE)

# Basic visualization
plot(clst, what = "classification")
abline(h = 1.96);abline(h = -1.96)
abline(v= 0)

# Full cluster plot - Figure s6.1
png(
  file.path(fig_dir,"figure_s6.1.png"),
  width = 1800,
  height = 1200
  )
par(cex = 20)
plot(clst, what = "classification",cex = 1.5)
dev.off()
mclust.options("classPlotColors") # display colors in order to id syndromes


# Cluster plots - Figure 5  ----------------------------------------------------

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
    )
  print(p)

  # Export file name
  file_name <- paste0("figure_5_",x_name,"_",y_name,".png")

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

  # Update cluster names to those in figure 6
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

raw_data <- lcbd_clst %>%
  rownames_to_column("COMID") %>%
  right_join(diff_cluster_data)

# Significant cat function. Coarsens LCBD values to sig levels
es_sig <-function(x) {
  ifelse(
    x >= 1.96,
    3,
    ifelse(
      x >= 1.64,
      2,
      ifelse(
        x >= 0,
        1,
        ifelse(
          x > -1.64,
          -1,
          ifelse(
            x > -1.96,
            -2,
            -3)
          )
        )
      )
    )

}


# Coarsened sig levels based on median delta LCBD ES - used for effects sizes
# figure 6.
cluster_med_raw <- raw_data %>%
  group_by(cluster) %>%
  summarise(
    across(-COMID,median,na.rm=T),
    n_sites = n()
  ) %>%
  mutate(across(-c(n_sites,cluster), es_sig)) %>%
  ungroup()


# Cluster species richness summary  --------------------------------------------

# Calculate taxonomic richness values for each cluster
comid_richness <- inv_df %>%
  mutate(COMID = as.character(COMID)) %>%
  left_join(diff_cluster_data) %>%
  mutate(
    cluster = case_when(
      is.na(cluster) ~ 7,
      T~cluster
    ),
    cluster = as.factor(cluster)
  )

# Cluster summary
comid_richness %>%
  group_by(cluster) %>%
  summarise(mean_native_rich = mean(Native_rich),
    med_native_rich = median(Native_rich),
    mean_nonnative_rich = mean(Nonnative_rich),
    med_nonnative_rich = median(Nonnative_rich),
    mean_nonnative_prop = mean(prop_nonnative),
    med_nonnative_prop = median(prop_nonnative))


# Cluster richness plots - Figure 6 barplots------------------------------------
for(i in 1:7){
  plot_data <- comid_richness %>%
    filter(cluster == i) %>%
    rename(
      N = Native_rich,
      NN = Nonnative_rich
    ) %>%
    select(-prop_nonnative,-cluster,-classification) %>%
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

  plot_file <- paste0("figure_6_rich_s",i,".png")
  ggsave(
    file.path(fig_dir,plot_file),
    p,
    width = 10,
    height = 10,
    dpi = 600
    )

}

# Syndrome species occurrences summary -----------------------------------------
# Tables in appendix 7

# !! requires raw community data - not provided !!
com_df <- readRDS("Diversity Input Data/raw_community_diversity_input.rds")

# Format raw species data
com_df$COMID <- as.character(com_df$COMID)
tbl_data <- com_df %>%
  left_join(diff_cluster_data) %>%
  mutate(
    cluster = case_when(
      is.na(cluster) ~ 7,
      T~cluster
    )
  )

# total sites per cluster
n_sites <- tbl_data %>%
  group_by(cluster) %>%
  summarise(tot_sites = n_distinct(COMID))

# Nonnative occurrences per cluster
nn_summary <- tbl_data %>%
  filter(Native8 == F) %>%
  group_by(cluster,Scientific_Name) %>%
  summarise(n_sites = n()) %>%
  left_join(n_sites) %>%
  mutate(prop_sites = n_sites/tot_sites) %>%
  dplyr::select(-tot_sites) %>%
  arrange(cluster,desc(prop_sites))

# Native occurrences per cluster
nat_summary <- tbl_data %>%
  filter(Native8 == T) %>%
  group_by(cluster,Scientific_Name) %>%
  summarise(n_sites = n())%>%
  left_join(n_sites) %>%
  mutate(prop_sites = n_sites/tot_sites) %>%
  dplyr::select(-tot_sites) %>%
  arrange(cluster,desc(prop_sites))

## Export ##
write.csv(nn_summary, file.path(fig_dir,"table_s7_nonnative.csv"))
write.csv(nat_summary, file.path(fig_dir,"table_s7_native.csv"))


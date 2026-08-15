# ------------------------------------------------------------------------------
#
#  SES output compilation
#
# ------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 04/28/2026

# Description : Compiles and formats the resulting SES, ES, and diagnostic
# stats across files. Also visualizes normality diagnostics to allow users to
# decide between SES or ES values for further analyses, and creates plots
# summarizing raw and effect size beta diversity values.


# House Keeping ----------------------------------------------------------------
rm(list=ls())

# Directories
out_dir <- "HPC copy/ses_outputs"
obs_dir <- "HPC copy/obs_out"
export_dir <- "Analysis_data"
fig_dir <- "Figures"

# Packages
library(dplyr)
library(purrr)
library(ggplot2)


# Load in outputs --------------------------------------------------------------

# all outputs
out_files <- list.files(out_dir)

# Metrics of interest
mets <- c(
  # "alpha",  # native alpha
  # "his_lcbd",  # native LCBD
  # "mod_lcbd",  # contemporary LCBD
  "_d_lcbd",  # change in lcbd
  "his_B",  # native beta diversity
  "mod_B",  # contemporary beta diversity
  "_d_B"  # change in beta diversity
  )

# Extract files for each metric of interest
file_list <- lapply(mets, function (x) out_files[grepl(x,out_files)])
names(file_list) <- mets


# Format outputs  --------------------------------------------------------------
out_list <- lapply(file_list, function(y) {

  # Load in files
  ls <- lapply(y,function(x) readRDS(file.path(out_dir,x)))
  names(ls) <- y

  # Transpose list so BLANK is on outside
  ls_t <- purrr::transpose(ls)

  # Format lcbd matrices
  if(all(grepl("lcbd",y))) {
    out <- lapply(ls_t, function(x){
      ls <- lapply(x, function(y) {
        as.data.frame(y) %>% tibble::rownames_to_column("COMID")
        }
        )
      purrr::reduce(ls,left_join,by=join_by(COMID)) %>%
        tibble::column_to_rownames("COMID")
    })
    return(out)
    }

  if(all(grepl("alpha",y)|grepl("_B",y))){
    out <- lapply(ls_t, function(x){
      df <-as.data.frame(do.call(cbind,x))
      colnames(df) <- sub("_ses_out.rds","",colnames(df))
      df
    })
    return(out)
  }

})

# Observed taxonomic data ------------------------------------------------------

# No null model exists for taxonomic beta diversity, so you need to manually
# merge in observed values

# load in obs tax beta
his_tax_obs <- readRDS(file.path(obs_dir,"his_tax_beta_obs.rds"))
mod_tax_obs <- readRDS(file.path(obs_dir,"mod_tax_beta_obs.rds"))

# Extract beta diversity components into matrix
his_tax_b <- do.call(cbind,lapply(his_tax_obs[[1]],as.vector))
mod_tax_b <- do.call(cbind,lapply(mod_tax_obs[[1]],as.vector))
colnames(his_tax_b) <- paste0("tax_his_",colnames(his_tax_b))
colnames(mod_tax_b) <- paste0("tax_mod_",colnames(mod_tax_b))

# Merge to beta data
out_list$his_B$obs <- cbind(his_tax_b,out_list$his_B$obs)  # Used for summary stat barplot
out_list$mod_B$obs <- cbind(mod_tax_b,out_list$mod_B$obs)  # Used for summary stat barplot

# Transpose output list
out_list_t <- purrr::transpose(out_list)


# Check for skew in null distributions  ----------------------------------------

# Skew plots
for (i in 1:length(out_list_t$skew)) {
  skew_df <- out_list_t$skew[[i]]
  boxplot(skew_df)
  abline(h=2)
  abline(h=-2)
  abline(h=0)
  title(names(out_list_t$skew)[i])
}

# Many null distributions are skewed, use empirical ES rather than SES


# Prepare diversity data for analyses  -----------------------------------------

# Subset to only include variables used in analysis
analysis_list <-out_list["_d_lcbd"]

final_list <-lapply(analysis_list, function (x) {

  # select observed and ES or SES
  sub.vec <- "empirical_es"
  ls <- x[sub.vec]

  # Add ses and es suffixes to column names
  colnames(ls$empirical_es) <- paste0(colnames(ls$empirical_es),"_es")

  # Merge obs, es, and ses of each measure
  ls <- lapply(ls, function(y) {
    as.data.frame(y) %>% tibble::rownames_to_column("COMID")
    })
  purrr::reduce(ls,left_join,by=join_by(COMID))
})

# Export
saveRDS(final_list$`_d_lcbd`,file.path(export_dir,"delta_lcbd.rds"))

# Summary stats and plots  -----------------------------------------------------

# The remainder of code is used to create figures 4 and s6.1-2 and summary
# stats used in the multidimensional beta change results section


# Summary Stats ----------------------------------------------------------------

# Percent change in observed mean pairwise beta
mpw_pct_list <-lapply(c("mod_B","his_B"), function(i){
  time_name = ifelse(i=="mod_B","Contemporary","Native")
  out_list_t$obs[[i]] %>%
    summarise(across(everything(),~mean(.x)))
})
do.call(function(x,y) 100*(1-x/y),mpw_pct_list)


# Observed delta LCBD summary
summary(out_list_t$obs$`_d_lcbd`)

# ES delta LCBD summary
summary(out_list_t$empirical_es$`_d_lcbd`)

# Percentage of sites with significant effect sizes
sig_summary <- lapply(names(out_list_t$empirical_es),function(x) {
  df <- out_list_t$empirical_es[[x]]
  n_obs <- nrow(df)
  dif <- df %>%
    summarise(
      # across(everything(), ~ sum(. >= 1.96, na.rm = TRUE)/n_obs)
      across(everything(), ~ sum(. >= 1.96, na.rm = TRUE))
    )

  non <- df %>%
    summarise(
      # across(everything(), ~ sum(. < 1.96 & . > -1.96, na.rm = TRUE)/n_obs)
      across(everything(), ~ sum(. < 1.96 & . > -1.96, na.rm = TRUE))
    )

  hom <- df %>%
    summarise(
      # across(everything(), ~ sum(. <= -1.96, na.rm = TRUE)/n_obs)
      across(everything(), ~ sum(. <= -1.96, na.rm = TRUE))
    )

  out <- rbind(hom,non,dif)
  row.names(out) <- c("hom","none","dif")
  out

})
names(sig_summary) <- names(out_list_t$empirical_es);sig_summary


# Beta barplot (Fig 4a) --------------------------------------------------------

# Extract and summarize pairwise beta diversity
mpw_list <-lapply(c("mod_B","his_B"), function(i){
  time_name = ifelse(i=="mod_B","Contemporary","Native")
  out_list_t$obs[[i]] %>%
    summarise(across(everything(),~mean(.x))) %>%
    mutate(time = time_name) %>%
    tidyr::pivot_longer(-time)
})

# Combine and format for plotting
mpw_df <- bind_rows(mpw_list) %>%
  mutate(
    dim = sub("_mod","",name),
    dim = sub("_his","",dim),
    facet = substr(name,1,3),
    component = sub(".*B", "", name),
    facet =factor(facet,c("tax","fun","phy")),
    time = factor(time, c("Native", "Contemporary"))
  ) %>%
  filter(component != "total")

# Plot
pair <- ggplot(mpw_df, aes(x = time, y = value, fill = dim)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~facet)+
  theme_classic(base_size = 25)+
  ylab("")+
  theme(legend.position = "none",
        axis.title.x=element_blank())+
  scale_fill_manual(values = c("#BF4354","#CCA6A7",
                               "#4A7B97","#A6B0BD",
                               "#CEA23F","#DCD0A7"))

# Export
ggsave(
  filename = file.path(fig_dir,"figure_4a.png"),
  plot = pair,
  width = 10,
  height = 6,
  dpi = 600
)


# Observed plots (Fig 4c) -----------------------------------------------------

# Subset data
obs_plot_df <- out_list_t$obs[["_d_lcbd"]] %>%
  tibble::rownames_to_column("COMID") %>%
  tidyr::pivot_longer(-COMID) %>%
  mutate(name = sub("_d","",name),
         name = sub("rf","",name))

# Set factor levels
plot_factor <- c("tax_Btotal","tax_Brepl","tax_Brich",
                 "fun_Btotal","fun_Brepl","fun_Brich",
                 "phy_Btotal","phy_Brepl","phy_Brich")
obs_plot_df$name <- factor(obs_plot_df$name,plot_factor)

# Plot
obs_box <- ggplot(
  obs_plot_df,
  aes(x = name, y = value, fill = name,color = name)
) +
  geom_boxplot()+
  stat_summary(
    geom = "crossbar",
    width = 0.75,
    fatten = 2,
    color = "black",
    fun.data = function(x){
      return(c(y = median(x), ymin = median(x), ymax = median(x)))
    }
  )+
  geom_hline(yintercept=0, linetype="dashed", color = "red",size = 1)+
  theme_classic(base_size = 25)+
  theme(
    legend.position = "none",
    panel.border =  element_rect(color = "black", fill = NA, size = 2.5),
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank()) +
  scale_y_continuous(breaks = pretty) +
  scale_fill_manual(values = c("#CEA23F","#CEA23F","#CEA23F",
                               "#BF4354","#BF4354","#BF4354",
                               "#4A7B97","#4A7B97","#4A7B97"))+
  scale_color_manual(values = c("#CEA23F","#CEA23F","#CEA23F",
                                "#BF4354","#BF4354","#BF4354",
                                "#4A7B97","#4A7B97","#4A7B97"))

# Export
ggsave(
  filename = file.path(fig_dir,"figure_4c.png"),
  plot = obs_box,
  width = 10,
  height = 6,
  dpi = 600
)


# Effect size plots (Fig 4b,d)  -----------------------------------------------
plot_metrics <-  c("_d_lcbd","_d_B")
lapply(1:2,function(i){
  # Subset and format data
  es_plot_df <- out_list_t$empirical_es[[plot_metrics[i]]] %>%
    tibble::rownames_to_column("COMID") %>%
    tidyr::pivot_longer(-COMID) %>%
    mutate(name = sub("_d","",name),
           name = sub("rf","",name))

  # Set factor levels
  plot_factor <- c("tax_Btotal","tax_Brepl","tax_Brich",
                   "fun_Btotal","fun_Brepl","fun_Brich",
                   "phy_Btotal","phy_Brepl","phy_Brich")
  es_plot_df$name <- factor(es_plot_df$name,plot_factor)


  # Plot
  es_box <- ggplot(
    es_plot_df,
    aes(x = name, y = value, fill = name,color = name)
  )+
    geom_boxplot()+
    stat_summary(
      geom = "crossbar",
      width = 0.75,
      fatten = 2,
      color = "black",
      fun.data = function(x){
        return(c(y = median(x), ymin = median(x), ymax = median(x)))
      }
    )+
    geom_hline(yintercept=0, linetype="dashed", color = "red",size = 1)+
    geom_hline(yintercept=1.96, linetype="dashed", color = "black",size = 1)+
    geom_hline(yintercept=-1.96, linetype="dashed", color = "black",size = 1)+
    theme_classic(base_size = 25)+
    theme(
      legend.position = "none",
      panel.border =  element_rect(color = "black", fill = NA, size = 2.5),
      axis.title.y=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank()) +
    scale_y_continuous(breaks = pretty) +
    scale_fill_manual(values = c(
      "#CEA23F","#CEA23F","#CEA23F",
      "#BF4354","#BF4354","#BF4354",
      "#4A7B97","#4A7B97","#4A7B97"))+
    scale_color_manual(values = c(
      "#CEA23F","#CEA23F","#CEA23F",
      "#BF4354","#BF4354","#BF4354",
      "#4A7B97","#4A7B97","#4A7B97"))

  #Export
  es_name <- c("figure_4d.png","figure_4b.png")[i]
  ggsave(
    filename = file.path(fig_dir,es_name),
    plot = es_box,
    width = 6.6,
    height = 6,
    dpi = 600
  )
})


#-------------------------------------------------------------------------------
#
#  Taxonomic Beta Null Model Preparation
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 03/24/2026

# Description: Prepares null and observed values for taxonomic beta diversity 
# into format that can be easily summarized. Null model files are compiled into 
# one, and change (delta) between contemporary and native species pool is 
# calculated for observed and each null iteration. Each metric (total, 
# replacement, lcbd, etc.) is exported as a single file for batch processing. 
# This script is for use on a HPC cluster.

# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(dplyr)
library(purrr)
library(parallel)

# Directories
obs_dir <- "HPC/obs_out"
null_dir <- "HPC/null_out"
out_dir <- "HPC/ses_inputs"


# Load in data  ----------------------------------------------------------------

mod_dim <- paste0("mod_tax_beta")
his_dim <- paste0("his_tax_beta")

# Load in observed data
mod_obs <- readRDS(file.path(obs_dir,paste0(mod_dim,"_obs.rds")))
his_obs <- readRDS(file.path(obs_dir,paste0(his_dim,"_obs.rds")))

# Extract null files
mod_null_files <- list.files(null_dir,mod_dim)

# Load in null files
mod_null_list <- lapply(
  mod_null_files, 
  function(x) readRDS(file.path(null_dir,x))
  )

# Combine null processing chunks into one list
mod_null <- unlist(mod_null_list,recursive = FALSE)


# Create input list for batch estimation of SES  -------------------------------

# Create obs input list for all beta components and LCBD matrices
mod_obs_input <- mod_obs[[1]]
his_obs_input <- his_obs[[1]]
mod_obs_input$lcbd <- mod_obs[[2]]
his_obs_input$lcbd <- his_obs[[2]]

# Create null input list for all beta components and LCBD matrices
mod_null_t <- purrr::transpose(mod_null)
mod_null_input <- purrr::transpose(mod_null_t[[1]])
mod_null_input$lcbd <- mod_null_t[[2]]


# Estimate the change in input values
d_obs_input <-  purrr::map2(mod_obs_input,his_obs_input,function(x,y) x-y)
d_null_input <- purrr::map2(mod_null_input,his_obs_input, function (x,y) {
    lapply(x, function(x) x-y)
  }
)

# Combine native only, contemporary, and delta inputs
names(his_obs_input) <- 
  sapply(names(his_obs_input),function(x) paste0("his_",x))
names(mod_obs_input) <- 
  sapply(names(mod_obs_input),function(x) paste0("mod_",x))
names(mod_null_input) <- 
  sapply(names(mod_null_input),function(x) paste0("mod_",x))
names(d_obs_input) <- 
  sapply(names(d_obs_input),function(x) paste0("d_",x))
names(d_null_input) <- 
  sapply(names(d_null_input),function(x) paste0("d_",x))
obs_input <- 
  c(d_obs_input)
null_input <- 
  c(d_null_input)

# Check that lists are in same order
all(names(obs_input) == names(null_input))


# Export -----------------------------------------------------------------------

# Export each metric separately for efficient ses batch processing
lapply(seq_len(length(obs_input)), function(i){
  out_list <- (list(obs=obs_input[[i]],null_list=null_input[[i]]))
  out_name <- paste0("tax_",names(obs_input)[i],"_ses_input.rds")
  saveRDS(out_list,file.path(out_dir,out_name))
})



#-------------------------------------------------------------------------------
#
#  Null model standardized effect size batch script
#
#-------------------------------------------------------------------------------

# Author: William K. Annis

# Created: 03/24/2026

# Description: Calculated standardized effects sizes, empirical p values, and
# diagnostic stats for a diversity metric defined by a shell script argument

# House keeping  ---------------------------------------------------------------
rm(list = ls())

# Packages
library(dplyr)

# Directories
input_dir <- "HPC/ses_inputs"
out_dir <- "HPC/ses_outputs"

# Load in custom SES function
source("HPC/00_null_model_effect_size_function.R")

# Load in data based on shell script argument
arg <- as.numeric(commandArgs(trailingOnly = TRUE))  # shell script argument
input_file <- list.files(input_dir)[arg]
input <- readRDS(file.path(input_dir,input_file))

# Run SES function and export  -------------------------------------------------
out <- ses_fun(
  obs =input$obs,
  null_list =input$null_list,
  diag = T, 
  emp_padding =1
  )

# Export
out_name <-sub("ses_input","ses_out",input_file)
saveRDS(out,file.path(out_dir,out_name))


## ─────────────────────────────────────────────────────────────────────────────
##
## Project: C:/Users/Corey/Documents/Statistics/PhD/Projects/BRAR
##
## Purpose of script: Configure packages, functions, options & themes
##
## Author: Corey Voller
##
## Date Created: 21-02-2025
##
## QC'd by:
## QC date:
##
## ─────────────────────────────────────────────────────────────────────────────
##
## Notes:
##   
##
## ─────────────────────────────────────────────────────────────────────────────
# Packages ---------------------------------------------------------------------
# List of packages to be used
packages <-
  c(
    "magrittr",
    "dplyr",
    "tidyr",
    "ggplot2",
    "gt",
    "grid",
    "stringr",
    "rpact",
    "parallel",
    "latex2exp",
    #"ggpubr",
    "scales",
    "patchwork",
    "data.table"
  )
# Install packages which aren't installed in "packages"
# if (length(packages[!(packages %in% installed.packages()[, "Package"])])) {
#   install.packages(packages[!(packages %in% installed.packages()[, "Package"])])
# }

# Load packages
lapply(packages, library, character.only = TRUE)

# Functions --------------------------------------------------------------------
# source functions from sub folder functions
file.sources = list.files(
  path=file.path(getwd(),"functions"),
  pattern = "\\.R$",
  full.names = TRUE,
  ignore.case = T
)
sapply(file.sources,source)


my_palette <- c(
  "#FED789FF",
  "#023743FF",
  "#72874EFF",
  "#476F84FF",
  "#A4BED5FF",
  "#453947FF"
)

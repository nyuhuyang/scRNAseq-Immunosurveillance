#library(ggplot2)
#library(Seurat)
#library(dplyr)
#library(magrittr)
path <- "Rshiny/MouseTumor/"
#source(paste0(path,"util.R"))
#(load(file="data/MouseTumor_2_20190219.Rda"))
#============== expression Rda ===============
samples <- c("Control","NAM")
#PrepareShiny(object, samples, path)
load("data/MouseTumor.Rda")



# Load libraries.
library(ggplot2)
library(survival)

# Set working directory.
setwd("C:/Users/Eric/Desktop/Bergeron/")

# Read and summarize the clinical variables.
source("SummarizeClinical.R")

#source("SummarizeExpression.R")

# Read in expression info which was summarized by SummarizeExpression.R
expr.data = read.table("output/Expression summary.txt", sep="\t", check.names = FALSE)

# Build a metadata table.
# Load patient barcode -> file name mapping.
file.manifest = read.table("TCGA/file_manifest.txt", sep="\t", header=TRUE)
file.names <- colnames(expr.data)[-c(1,2)]
metadata = file.manifest[match(file.names, as.character(file.manifest$File.Name)),]
metadata$Type = ifelse(substr(metadata$Sample, 14, 14)=="1", "Healthy", "Tumor")

# Get a list of sample barcodes where tumor tissue as present within normal tissue
to.exclude = gsub("[A-Z]$", "",subset(file.annotations, Annotation=="Normal tissue contains tumor")$Item.Barcode)
metadata$Exclude = metadata$Sample %in% to.exclude
metadata = cbind(metadata, selectedDF[match(gsub("-..$", "", metadata$Sample), as.character(selectedDF$patient.barcode)),])


# Entrez IDs of interest
# First set.
#goi <- c(27181, 3568, 1232, 5553, 6036, 8288, 6037, 1178, 6356, 6369, 10344)

# Second set.
#goi <- c(3567, 3596, 6352, 4283, 3627, 90865, 916, 920, 925, 3458, 7124, 3002)

# All genes of interest (eosinophiles)
goi.eosino.1 <- c(27181, 3568, 1232, 5553, 6036, 8288, 6037, 1178, 6356, 6369, 10344)
goi.eosino.2 <- c(3567, 3596, 6352, 4283, 3627, 90865, 916, 920, 925, 3458, 7124, 3002)
         
# All genes of interest (dendrites)
goi.dendrites <- c(6346, 9308, 30835, 1235, 913, 909, 941, 942, 1236, 6367, 958)

# Get gene id for positive controls from Talantov et al, 2010.
pos.control.str = c("ACTG2", "CALD1", "CBX3", "DCHS1", "DKK3", "DPT", "FLNA", "FLNC", "GAS1", "GSN", "HIST1H3D", "LIMS2", "LMOD1", "MT1X", "MYH11", "MYLK", "PDLIM3", "PDLIM7", "RASL12", "SH3BGRL", "SMTN", "SORBS1", "SSBP1", "TNS1")
goi.pos <- expr.data$Entrez.Gene.ID[expr.data$Symbol %in% pos.control.str]

goi.glucose <- c(6513, 6515, 6517, 29988, 66035, 154091, 6523, 3098, 3099, 3101, 2645)

goi.ac.loss <- c(4907, 7779, 170691, 56971, 3269, 116441, 84525, 203102, 2131, 54840)
goi.ac.gain <- c(149233, 170589, 220323, 4088, 23192, 1739, 387129, 7010, 54980)

goi.all <- c(goi.eosino.1, goi.eosino.2, goi.dendrites, goi.pos, goi.glucose, goi.ac.loss, goi.ac.gain)

goi.list = list(Eosiniphiles.1   = goi.eosino.1,
                Eosiniphiles.2   = goi.eosino.2,
                Dendrites        = goi.dendrites,
                Pos.control      = as.integer(as.character(goi.pos)),
                Glucose          = goi.glucose,
                Acetylation.loss = goi.ac.loss,
                Acetylation.gain = goi.ac.gain)

source("Survival analysis.R")
source("GleasonGLM.R")

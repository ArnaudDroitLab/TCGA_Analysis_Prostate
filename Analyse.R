# Load libraries.
library(ggplot2)

# Set working directory.
setwd("C:/Users/Eric/Desktop/Bergeron/")

# Get the list of all expression files.
expr.path = "TCGA/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3"
genes.normalized.file <- dir(expr.path, pattern="*.rsem.genes.normalized_results")

# Read the normalized counts from all files.
results <- NULL
for(file.name in genes.normalized.file) {
   file.content = read.table(file.path(expr.path, file.name), sep="\t", header=TRUE, stringsAsFactors=FALSE)
   #uuid = gsub("unc.edu.", "", gsub(".rsem.genes.normalized_results", "", file.name))$Sample
   uuid = file.name
   
   if(is.null(results)) {
       # For the first file, get the gene IDs as well as the normalized counts.
       id.vector = unlist(strsplit(file.content$gene_id, "|", fixed=TRUE))
       symbols.vector = id.vector[seq(1, length(id.vector) - 1, by=2)]
       entrez.id.vector = id.vector[seq(2, length(id.vector), by=2)]

       results <- data.frame(Symbol = symbols.vector, Entrez.Gene.ID = entrez.id.vector, file.content$normalized_count)
       colnames(results)[3] <- uuid
   } else {
       results[,uuid] <- file.content$normalized_count
   }
}

# Load patient barcode -> file name mapping.
file.manifest = read.table("TCGA/file_manifest.txt", sep="\t", header=TRUE)
file.names <- colnames(results)[-c(1,2)]
barcodes <- file.manifest$Sample[match(file.names, file.manifest$File.Name)]
colnames(results)[-c(1,2)] <- as.character(barcodes)

# Entrez IDs of interest
# First set.
#goi <- c(27181, 3568, 1232, 5553, 6036, 8288, 6037, 1178, 6356, 6369, 10344)

# Second set.
goi <- c(3567, 3596, 6352, 4283, 3627, 90865, 916, 920, 925, 3458, 7124, 3002)

# Generate a boxplot of the genes of interest.
boxplot.matrix <- t(results[results$Entrez.Gene.ID %in% goi, c(-1, -2)])
boxplot.matrix[boxplot.matrix <= 1] <- 1
colnames(boxplot.matrix) <- results$Symbol[as.integer(colnames(boxplot.matrix))]
pdf("Boxplot-2.pdf", width=14, height=7)
boxplot(log2(boxplot.matrix))
dev.off()


resultMatrix <- as.matrix(results[,-c(1,2)])

# Calculate the ECDF curve for all genes, and get its values at 0.5 (on a log2 scale) intervals.
number.of.points = length(seq(0, 23, by=0.5))
resultVec = vector(mode="numeric", nrow(results) * number.of.points)
for(i in 1:nrow(results)) {
   resultVec[ (((i-1)*number.of.points) + 1):(i*number.of.points) ] <- ecdf(log2(resultMatrix[i,]))(seq(0, 23, by=0.5))
}

# Get those results in a data-frame for displaying them using ggplot.
resultDF <- data.frame(ECDF=resultVec,
                       Gene=rep(results[,2], each=number.of.points),
                       Expression=seq(0, 23, by=0.5))


# Do the same thing, but just for the genes of interest. We'll use this DF
# for a second set of curves with a different styling.
number.of.points = length(seq(0, 23, by=0.5))
resultVec = vector(mode="numeric", length(goi) * number.of.points)
for(i in 1:length(goi)) {
   
   resultVec[ (((i-1)*number.of.points) + 1):(i*number.of.points) ] <- ecdf(log2(resultMatrix[as.character(results$Entrez.Gene.ID)==as.character(goi[i]),]))(seq(0, 23, by=0.5))
}

goiDF <- data.frame(ECDF=resultVec,
                    Gene=factor(rep(results$Symbol[match(goi, results$Entrez.Gene.ID)], each=number.of.points)),
                    Expression=seq(0, 23, by=0.5))

# Plot the ECDF curves and the labels.
ecdfPlot <- ggplot() +
    geom_line(data=resultDF, mapping=aes(x=Expression, y=ECDF, group=Gene), alpha=0.02) + 
    geom_line(data=goiDF, mapping=aes(x=Expression, y=ECDF, group=Gene, color=Gene), size=1) +
    geom_text(data=goiDF[goiDF$Expression==0,], mapping=aes(x=Expression, y=ECDF, label=Gene), hjust=1.1) +
    xlim(c(-1.5, 20)) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
         
ggsave("ECDF-2.pdf", ecdfPlot, width=14, height=7)          



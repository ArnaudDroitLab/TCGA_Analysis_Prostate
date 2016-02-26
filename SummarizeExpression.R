# Set work directory.
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
       entrez.id.vector = as.integer(id.vector[seq(2, length(id.vector), by=2)])

       results <- data.frame(Symbol = symbols.vector, Entrez.Gene.ID = entrez.id.vector, file.content$normalized_count)
       colnames(results)[3] <- uuid
   } else {
       results[,uuid] <- file.content$normalized_count
   }
}
rownames(results) = results$Entrez.Gene.ID

write.table(results, "output/Expression summary.txt", sep="\t")
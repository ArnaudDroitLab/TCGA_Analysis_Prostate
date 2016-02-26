# Load libraries.
library(ggplot2)
library(survival)

# Set working directory.
setwd("C:/Users/Eric/Desktop/Bergeron/")

source("SummarizeClinical.R")

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
metadata = file.manifest[match(file.names, file.manifest$File.Name),]
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
goi.eosino <- c(27181, 3568, 1232, 5553, 6036, 8288, 6037, 1178, 6356, 6369, 10344,
                3567, 3596, 6352, 4283, 3627, 90865, 916, 920, 925, 3458, 7124, 3002)
         
# All genes of interest (dendrites)
goi.dendrites <- c(6346, 9308, 30835, 1235, 913, 909, 941, 942, 1236, 6367, 958)
         
         
pos.control.str = "ACTG2
CALD1
CBX3
DCHS1
DKK3
DPT
FLNA
FLNC
GAS1
GSN
HIST1H3D
LIMS2
LMOD1
MT1X
MYH11
MYLK
PDLIM3
PDLIM7
RASL12
SH3BGRL
SMTN
SORBS1
SSBP1
TNS1"
         
pos.control = unlist(strsplit(pos.control.str, split="\n"))
         
goi.pos <- results$Entrez.Gene.ID[results$Symbol %in% pos.control]

goi <- c(goi.eosino, goi.dendrites, as.integer(as.character(goi.pos)))
         
dir.create("output/profiles/", recursive=TRUE, showWarnings=FALSE)
relapse.stats = data.frame(Entrez=character(0), Symbol=character(0), p.val=numeric(0), conf.low=numeric(0), conf.high=numeric(0))
gleason.stats = data.frame(Entrez=character(0), Symbol=character(0), p.val=numeric(0), coef=numeric(0))
cor.stats = data.frame(Entrez=character(0), Symbol=character(0), Expr.Relapse=numeric(0), Expr.Gleason=numeric(0), Gleason.Relapse=numeric(0))
for(gene in goi) {
    # Get expression levels.
    expr.levels = unlist(results[as.character(results$Entrez.Gene.ID) == gene, -c(1,2)])
    
    expDF = data.frame(Expression=expr.levels,
                       TissueType=metadata$Type,
                       Relapse=metadata$years.to.relapse,
                       Gleason=metadata$gleason.category,
                       Exclude=metadata$Exclude)
                       
    symbol = results$Symbol[results$Entrez.Gene.ID == gene]
    dir.create(file.path("output/profiles/", symbol), recursive=TRUE, showWarnings=FALSE)
    
    # Analyze relationship to chemical relapse.
    exp.subset = subset(expDF, Exclude==FALSE & TissueType=="Tumor" & !is.na(Relapse))
    ggplot(exp.subset, aes(x=Relapse, y=log2(Expression + 1))) + geom_point()
    ggsave(paste0("output/profiles/", symbol, "/", gene, " - ", symbol, " vs biochemical.relapse.pdf"))
    
    #t.results = t.test(log2(subset(exp.subset, Relapse=="YES")$Expression + 1),
    #                   log2(subset(exp.subset, Relapse=="NO")$Expression + 1))
    cor.expr.relapse = cor(as.numeric(exp.subset$Relapse), log2(exp.subset$Expression + 1), use="pairwise.complete.obs")
                       
#    relapse.stats = rbind(relapse.stats,
#                          data.frame(Entrez=gene,
#                                     Symbol=symbol,
#                                     p.val=t.results$p.value,
#                                     conf.low=t.results$conf.int[1],
#                                     conf.high=t.results$conf.int[2]))

    # Analyze relationship to gleason score
    exp.subset = subset(expDF, Exclude==FALSE & TissueType=="Tumor" & !is.na(Gleason))
    ggplot(exp.subset, aes(x=Gleason, y=log2(Expression + 1))) + geom_boxplot()
    ggsave(paste0("output/profiles/", symbol, "/", gene, " - ", symbol, " vs gleason.score.pdf"))
    
    expr.glm = glm(Expression ~ as.numeric(Gleason), data=exp.subset)
    cor.expr.gleason = cor(as.numeric(exp.subset$Gleason), log2(exp.subset$Expression + 1), use="pairwise.complete.obs")
    gleason.stats = rbind(gleason.stats,
                          data.frame(Entrez=gene,
                                     Symbol=symbol,
                                     p.val=coef(summary(expr.glm))["as.numeric(Gleason)",4],
                                     coef=expr.glm$coefficients["as.numeric(Gleason)"]))

    # Caculate relapse -> gleason correlation as a basis for comparison.
    cor.gleason.relapse = cor(as.numeric(expDF$Relapse), as.numeric(expDF$Gleason), use="pairwise.complete.obs")
                                     
    cor.stats = rbind(cor.stats, 
                      data.frame(Entrez=gene,
                                 Symbol=symbol,
                                 Expr.Relapse=cor.expr.relapse,
                                 Expr.Gleason=cor.expr.gleason,
                                 Gleason.Relapse=cor.gleason.relapse))
}
write.table(relapse.stats, "output/profiles/Relapse stats.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
write.table(gleason.stats, "output/profiles/Gleason stats.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
write.table(cor.stats, "output/profiles/Correlation stats.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
         
         
# # Generate a boxplot of the genes of interest.
# boxplot.matrix <- t(results[results$Entrez.Gene.ID %in% goi, c(-1, -2)])
# boxplot.matrix[boxplot.matrix <= 1] <- 1
# colnames(boxplot.matrix) <- results$Symbol[as.integer(colnames(boxplot.matrix))]
# pdf("Boxplot-2.pdf", width=14, height=7)
# boxplot(log2(boxplot.matrix))
# dev.off()
# 
# 
# resultMatrix <- as.matrix(results[,-c(1,2)])
# 
# # Calculate the ECDF curve for all genes, and get its values at 0.5 (on a log2 scale) intervals.
# number.of.points = length(seq(0, 23, by=0.5))
# resultVec = vector(mode="numeric", nrow(results) * number.of.points)
# for(i in 1:nrow(results)) {
#    resultVec[ (((i-1)*number.of.points) + 1):(i*number.of.points) ] <- ecdf(log2(resultMatrix[i,]))(seq(0, 23, by=0.5))
# }
# 
# # Get those results in a data-frame for displaying them using ggplot.
# resultDF <- data.frame(ECDF=resultVec,
#                        Gene=rep(results[,2], each=number.of.points),
#                        Expression=seq(0, 23, by=0.5))
# 
# 
# # Do the same thing, but just for the genes of interest. We'll use this DF
# # for a second set of curves with a different styling.
# number.of.points = length(seq(0, 23, by=0.5))
# resultVec = vector(mode="numeric", length(goi) * number.of.points)
# for(i in 1:length(goi)) {
#    
#    resultVec[ (((i-1)*number.of.points) + 1):(i*number.of.points) ] <- ecdf(log2(resultMatrix[as.character(results$Entrez.Gene.ID)==as.character(goi[i]),]))(seq(0, 23, by=0.5))
# }
# 
# goiDF <- data.frame(ECDF=resultVec,
#                     Gene=factor(rep(results$Symbol[match(goi, results$Entrez.Gene.ID)], each=number.of.points)),
#                     Expression=seq(0, 23, by=0.5))
# 
# # Plot the ECDF curves and the labels.
# ecdfPlot <- ggplot() +
#     geom_line(data=resultDF, mapping=aes(x=Expression, y=ECDF, group=Gene), alpha=0.02) + 
#     geom_line(data=goiDF, mapping=aes(x=Expression, y=ECDF, group=Gene, color=Gene), size=1) +
#     geom_text(data=goiDF[goiDF$Expression==0,], mapping=aes(x=Expression, y=ECDF, label=Gene), hjust=1.1) +
#     xlim(c(-1.5, 20)) +
#     theme(axis.line = element_line(colour = "black"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank()) 
#          
# ggsave("ECDF-2.pdf", ecdfPlot, width=14, height=7)          

# Build a linear model combining SIGLEC8 and IL5RA
expDF = data.frame(SIGLEC8=unlist(subset(results, Symbol == "SIGLEC8")[-c(1,2)]),
                   IL5RA=unlist(subset(results, Symbol == "IL5RA")[-c(1,2)]),
                   TissueType=metadata$Type,
                   Relapse=metadata$biochemical.relapse,
                   Gleason=metadata$gleason.score,
                   Exclude=metadata$Exclude)

glm.fit = glm(Relapse ~ log2(SIGLEC8+1) + log2(IL5RA+1), data=subset(expDF, TissueType=="Tumor" & !Exclude & !is.na(Relapse)), family="binomial")
sink("output/SIGLEC8 + IL5RA vs Relapse GLM.txt")
summary(glm.fit)
sink(NULL)
                   
glm.fit = glm(as.numeric(Gleason) ~ log2(SIGLEC8+1) + log2(IL5RA+1), data=subset(expDF, TissueType=="Tumor" & !Exclude & !is.na(Gleason)))
sink("output/SIGLEC8 + IL5RA vs Gleason GLM.txt")
summary(glm.fit)
sink(NULL)                   

# Test various combinations of SIGLEC8/IL5RA dichotomization and see which ones gie the most significant results.
for(i in seq(0.2, 0.8, by=0.1)) {
    expDF.copy1 = expDF
    
    expDF.copy1$SIGLEC8 <- factor((expDF$SIGLEC8 >= quantile(expDF$SIGLEC8, i)))
    for(j in seq(0.2, 0.8, by=0.1)) {
        expDF.copy1$IL5RA <- factor((expDF$IL5RA >= quantile(expDF$IL5RA, j)))
        
        cat(i, ", ", j, ":\n")
        glm.fit = glm(Relapse ~ SIGLEC8 + IL5RA, data=subset(expDF.copy1, TissueType=="Tumor" & !Exclude & !is.na(Relapse)), family="binomial")
        if(all(coef(summary(glm.fit))[,4] < 0.05)) {
            cat("Relapse:\n")
            print(summary(glm.fit))
        }
        
        
        glm.fit = glm(as.numeric(Gleason) ~ SIGLEC8 + IL5RA, data=subset(expDF.copy1, TissueType=="Tumor" & !Exclude & !is.na(Gleason)))
        if(all(coef(summary(glm.fit))[,4] < 0.05)) {
            cat("Gleason:\n")
            print(summary(glm.fit))
        }
    }
}     

                  

                             
goi.expr = log2(t(results[as.character(results$Entrez.Gene.ID) %in% goi, -c(1,2)]) + 1)
colnames(goi.expr) = as.character(results[as.character(results$Entrez.Gene.ID) %in% goi,"Symbol"])
expDF = cbind(goi.expr, metadata)

resultsDF = data.frame(Entrez.ID=numeric(0), Symbol=character(0),
                       uni.pval=numeric(0), uni.conf.low=numeric(0), uni.conf.high=numeric(0),
                       multi.pval=numeric(0), multi.conf.low=numeric(0), multi.conf.high=numeric(0))

for(gene in goi) {
    symbol = results$Symbol[results$Entrez.Gene.ID == gene]

    uni.model.str = paste0("Surv(time.to.event, event.type) ~ ", symbol)
    uni.model = coxph(as.formula(uni.model.str), data=expDF)
    
    multi.model.str = paste0("Surv(time.to.event, event.type) ~ as.numeric(gleason.category) + as.numeric(stage) + ", symbol)
    multi.model = coxph(as.formula(multi.model.str), data=expDF)
    
    resultsDF <- rbind(resultsDF,
                       data.frame(Entrez.ID=gene,
                                  Symbol=symbol,
                                  uni.pval=coefficients(summary(uni.model))[as.character(symbol), "Pr(>|z|)"],
                                  uni.conf.low=summary(uni.model)$conf.int[as.character(symbol), "lower .95"],
                                  uni.conf.high=summary(uni.model)$conf.int[as.character(symbol), "upper .95"],
                                  multi.pval=coefficients(summary(multi.model))[as.character(symbol), "Pr(>|z|)"],
                                  multi.conf.low=summary(multi.model)$conf.int[as.character(symbol), "lower .95"],
                                  multi.conf.high=summary(multi.model)$conf.int[as.character(symbol), "upper .95"]))
}

    

cox.fit = coxph(Surv(time.to.event, event.type) ~ SIGLEC8, data=expDF)

cox.fit = coxph(Surv(time.to.event, event.type) ~ as.numeric(gleason.category) + as.numeric(stage) + CCL1 + CCL22 + CCR6 + CCR7 + CD1A + CD1E + CD209 + CD40 + CD80 + CD83 + CD86, data=expDF)
cox.fit = coxph(Surv(time.to.event, event.type) ~ as.numeric(gleason.category) + as.numeric(stage) + CD80 + CD83 + CD86, data=expDF)

cox.fit = coxph(Surv(time.to.event, event.type) ~ as.numeric(gleason.category) + as.numeric(stage) + ACTG2+ CALD1+CBX3+DCHS1+DKK3+DPT+FLNA+FLNC+GAS1+GSN+HIST1H3D+LIMS2+LMOD1+MT1X+MYH11+MYLK+PDLIM3+PDLIM7+RASL12+SH3BGRL+SMTN+SORBS1+SSBP1+TNS1, data=expDF)

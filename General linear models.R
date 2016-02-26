         
dir.create("output/profiles/", recursive=TRUE, showWarnings=FALSE)
relapse.stats = data.frame(Entrez=character(0), Symbol=character(0), p.val=numeric(0), conf.low=numeric(0), conf.high=numeric(0))
gleason.stats = data.frame(Entrez=character(0), Symbol=character(0), p.val=numeric(0), coef=numeric(0))
cor.stats = data.frame(Entrez=character(0), Symbol=character(0), Expr.Relapse=numeric(0), Expr.Gleason=numeric(0), Gleason.Relapse=numeric(0))
for(gene in goi) {
    # Get expression levels.
    expr.levels = unlist(expr.data[expr.data$Entrez.Gene.ID == gene, -c(1,2)])
    
    expDF = data.frame(Expression=expr.levels,
                       TissueType=metadata$Type,
                       Relapse=metadata$biochemical.relapse,
                       Gleason=metadata$gleason.score,
                       Exclude=metadata$Exclude)
                       
    symbol = expr.data$Symbol[expr.data$Entrez.Gene.ID == gene]
    dir.create(file.path("output/profiles/", symbol), recursive=TRUE, showWarnings=FALSE)
    
    # Analyze relationship to chemical relapse.
    exp.subset = subset(expDF, Exclude==FALSE & TissueType=="Tumor" & !is.na(Relapse))
    ggplot(exp.subset, aes(x=Relapse, y=log2(Expression + 1))) + geom_point()
    ggsave(paste0("output/profiles/", symbol, "/", gene, " - ", symbol, " vs biochemical.relapse.pdf"))
    
    t.expr.data = t.test(log2(subset(exp.subset, Relapse=="YES")$Expression + 1),
                       log2(subset(exp.subset, Relapse=="NO")$Expression + 1))
    cor.expr.relapse = cor(as.numeric(exp.subset$Relapse), log2(exp.subset$Expression + 1), use="pairwise.complete.obs")
                       
    relapse.stats = rbind(relapse.stats,
                          data.frame(Entrez=gene,
                                     Symbol=symbol,
                                     p.val=t.expr.data$p.value,
                                     conf.low=t.expr.data$conf.int[1],
                                     conf.high=t.expr.data$conf.int[2]))

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

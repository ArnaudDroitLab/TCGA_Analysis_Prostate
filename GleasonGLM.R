library(preprocessCore)

#normed = normalize.quantiles(as.matrix(log2(expr.data[, -c(1,2)] + 1)))
goi.expr = log2(t(expr.data[expr.data$Entrez.Gene.ID %in% goi.all, -c(1,2)]) + 1)
#goi.expr = t(normed[expr.data$Entrez.Gene.ID %in% goi.all,])
colnames(goi.expr) = as.character(expr.data[expr.data$Entrez.Gene.ID %in% goi.all, "Symbol"])
expDF = cbind(goi.expr, metadata)

for(i in 1:length(goi.list)) {
    goi = goi.list[[i]]
    label = names(goi.list)[i]

    resultsDF = data.frame(Entrez.ID=numeric(0), Symbol=character(0),
                           GS.3.4.estimate=numeric(0), GS.3.4.pval=numeric(0),
                           GS.4.3.estimate=numeric(0), GS.4.3.pval=numeric(0),
                           GS.8.estimate=numeric(0),   GS.8.pval=numeric(0))
    
    for(gene in goi) {
        symbol = expr.data$Symbol[expr.data$Entrez.Gene.ID == gene]
        output.dir = file.path("output/gleason.glm/", label, symbol)
        dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)
    
        uni.model.str = paste0(symbol, "~ gleason.category.4")
        uni.model = glm(as.formula(uni.model.str), data=expDF[!is.na(expDF$gleason.category),])
        sink(file.path(output.dir, "GLM.txt"))
        print(summary(uni.model))
        sink(NULL)
        

        ggplot(expDF[!is.na(expDF$gleason.category),], mapping=aes_string(y=as.character(symbol), x="gleason.category.4")) + geom_boxplot()
        ggsave(file.path(output.dir, "boxplot.pdf"), width=7, height=7)
        
        resultsDF <- rbind(resultsDF,
                           data.frame(Entrez.ID=gene,
                                      Symbol=symbol,
                                      GS.3.4.estimate =coef(summary(uni.model))[,"Estimate"][2],
                                      GS.3.4.pval     =coef(summary(uni.model))[,"Pr(>|t|)"][2],
                                      GS.4.3.estimate =coef(summary(uni.model))[,"Estimate"][3],
                                      GS.4.3.pval     =coef(summary(uni.model))[,"Pr(>|t|)"][3],
                                      GS.8.estimate   =coef(summary(uni.model))[,"Estimate"][4],
                                      GS.8.pval       =coef(summary(uni.model))[,"Pr(>|t|)"][4]))
    }

    write.table(resultsDF, file=paste0("output/gleason.glm/", label, ".txt"), col.names=TRUE, row.names=FALSE, sep="\t")
}    


# Library survival

# Extract the genes of interest from the 
goi.expr = log2(t(expr.data[expr.data$Entrez.Gene.ID %in% goi.all, -c(1,2)]) + 1)
colnames(goi.expr) = as.character(expr.data[expr.data$Entrez.Gene.ID %in% goi.all, "Symbol"])
expDF = cbind(goi.expr, metadata)

for(i in 1:length(goi.list)) {
    goi = goi.list[[i]]
    label = names(goi.list)[i]

    resultsDF = data.frame(Entrez.ID=numeric(0), Symbol=character(0),
                           uni.pval=numeric(0), uni.conf.low=numeric(0), uni.conf.high=numeric(0),
                           multi.pval=numeric(0), multi.conf.low=numeric(0), multi.conf.high=numeric(0))
    
    for(gene in goi) {
        symbol = expr.data$Symbol[expr.data$Entrez.Gene.ID == gene]
        output.dir = file.path("output/survival/", label, symbol)
        dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)
    
        uni.model.str = paste0("Surv(time.to.event, event.type) ~ ", symbol)
        uni.model = coxph(as.formula(uni.model.str), data=expDF)
        sink(file.path(output.dir, "Univariate model.txt"))
        summary(uni.model)
        sink(NULL)
        
        multi.model.str = paste0("Surv(time.to.event, event.type) ~ as.numeric(gleason.category) + as.numeric(stage) + ", symbol)
        multi.model = coxph(as.formula(multi.model.str), data=expDF)
        sink(file.path(output.dir, "Multivariate model.txt"))
        summary(multi.model)
        sink(NULL)
        
        resultsDF <- rbind(resultsDF,
                           data.frame(Entrez.ID=gene,
                                      Symbol=symbol,
                                      univariate.pval=coefficients(summary(uni.model))[as.character(symbol), "Pr(>|z|)"],
                                      univariate.conf.low=summary(uni.model)$conf.int[as.character(symbol), "lower .95"],
                                      univariate.conf.high=summary(uni.model)$conf.int[as.character(symbol), "upper .95"],
                                      multivariate.pval=coefficients(summary(multi.model))[as.character(symbol), "Pr(>|z|)"],
                                      multivariate.conf.low=summary(multi.model)$conf.int[as.character(symbol), "lower .95"],
                                      multivariate.conf.high=summary(multi.model)$conf.int[as.character(symbol), "upper .95"]))
    }

    write.table(resultsDF, file=paste0("output/Cox regression for ", label, ".txt"), col.names=TRUE, row.names=FALSE, sep="\t")
    
}    

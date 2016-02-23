library(TCGAbiolinks)

setwd("C:/Users/Eric/Desktop/Bergeron/")
dir.create("output")

# Load pateint-focused clinical data suing TCGAbiolinks
patient <- TCGAquery_clinic("PRAD", "clinical_patient")
follow_up <- TCGAquery_clinic("PRAD", "clinical_follow_up_v1.0")
biospecimen_slide <- TCGAquery_clinic("PRAD", "biospecimen_slide")

# biospecimens slides contain both normal and tumor tissue, and
# some patients ahve more than one slide. We need to summarize all of this.
biospecimen <- biospecimen_slide[substr(biospecimen_slide$bcr_sample_barcode, 14, 14)==0,c(3,8:15)]
for(i in 2:ncol(biospecimen)) {
    biospecimen[,i] <- as.numeric(biospecimen[,i])
}

agg.biospecimen = aggregate(biospecimen[,-1], by=list(biospecimen$bcr_sample_barcode), mean, na.rm=TRUE)
colnames(agg.biospecimen)[1] <- "Barcode"
agg.biospecimen$Barcode = gsub("-...$", "", agg.biospecimen$Barcode)

# Laod tumor-sample information separately, since TCGAbiolinks does not suport it
# because it is keyed off the tumor uid rather than the pateint uid.
tumor_sample <- read.table("input/ssf_tumor_samples_prad.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")

# Read file annotation to identify problem cases.
file.annotations <- read.table("TCGA/file_annotations.txt", sep="\t", header=TRUE)

# Identify patients which had no tumors.
tumor.free.patients = subset(file.annotations, Item.Type=="Patient" & Category=="Pathology outside specification")$Item.Barcode

# Identify normal samples which contained tumors.
tumor.in.normal.tissue = subset(file.annotations, Item.Type=="Sample" & Category=="Pathology outside specification")$Item.Barcode


# Put everything into one big, happy data-frame.
allDF <- cbind(patient,
               follow_up[match(patient$bcr_patient_uuid, follow_up$bcr_patient_uuid),],
               agg.biospecimen[match(patient$bcr_patient_barcode, agg.biospecimen$Barcode),],
               tumor_sample[match(patient$bcr_patient_barcode, tumor_sample$bcr_patient_barcode),])

# Remvoe no-tumor patients.
allDF <- subset(allDF, !(bcr_patient_barcode %in% tumor.free.patients))
               
levels_except_NA <- function(x) {
    return(sort(setdiff(unique(x), c("[Not Available]", "[Not Applicable]", "[Discrepancy]", "[Unknown]"))))
}

factorize.gleason.category <- function(primary.pattern, secondary.pattern) {
    results <- rep(NA, length(primary.pattern))
    results[as.numeric(primary.pattern) + as.numeric(secondary.pattern) <= 6] <- "<=6"
    results[as.numeric(primary.pattern) == 3 & as.numeric(secondary.pattern) == 4] <- "=3, 4"
    results[as.numeric(primary.pattern) == 4 & as.numeric(secondary.pattern) == 3] <- "=4, 3"
    results[as.numeric(primary.pattern) + as.numeric(secondary.pattern) >= 8] <- ">=8"
    
    return(factor(results))    
}

factorize.clinical.metastasis <- function(anatomic.site, site.text) {
   results <- rep(NA, length(anatomic.site))
   
   results[anatomic.site == "Bone"] <- "Bone"
   results[anatomic.site == "Other, specify" & (site.text == "Bladder" | site.text == "Lung and Liver" | site.text == "Neck")] <- "Distant"
   results[anatomic.site == "Lung"] <- "Distant"
   results[anatomic.site == "Other, specify" & (site.text == "Biochemical" | site.text == "biochemical recurrence only")] <- "Biochemical"
   
   return(factor(results))
}

factorize.death <- function(vital.status, death.reason) {
   results <- rep(NA, length(vital.status))
   results[vital.status=="Alive"] <- "Alive"
   results[death.reason=="Prostate Cancer"] <- "Dead (Prostate cancer)"
   results[death.reason=="[Unknown]" | death.reason=="Other, non-malignant disease"] <- "Dead (Unknown)"

   return(factor(results))
}

factorize.castrate.resistant <- function(days.to.second.bc.rec) {
    results <- rep("NO", length(days.to.second.bc.rec))
    results[days.to.second.bc.rec!="[Not Available]"] <- "YES"
    
    return(factor(results))
}

factorize_percent <- function(x) {
    cut(as.numeric(x), c(0, 1, 10, 50, 100), right=FALSE, include.lowest=TRUE)
}

allDF$pathologic_T_no_T <- gsub("^T", "", allDF$pathologic_T)

selectedDF <- with(allDF, data.frame(
    # Patient ID
    patient.id = bcr_patient_uuid,
    
    # Données démographiques et cliniques de base (au moment de la chirurgie)
    age              = -as.numeric(days_to_birth)/365,
    follow.up        = as.numeric(days_to_last_followup)/365,
    race             = factor(race, levels_except_NA(race)),
    tumor.type       = factor(histological_type, levels_except_NA(histological_type)),
    psa.preop        = cut(as.numeric(psa_result_preop), c(0, 4,10,20, Inf), right=FALSE, include.lowest=TRUE),
    gleason.score    = factor(as.numeric(gleason_score)),
    gleason.category = factorize.gleason.category(primary_pattern, secondary_pattern),
    stage            = factor(pathologic_T_no_T, levels_except_NA(pathologic_T_no_T)),
    stage.category   = factor(gsub("[abc]$", "", pathologic_T_no_T), levels_except_NA(gsub("[abc]$", "", pathologic_T_no_T))),
    limph.node       = factor(gsub("[1-9]", "x", pathologic_N), levels_except_NA(gsub("[1-9]", "x", pathologic_N))),
    metastases       = factor(gsub("1[a-c]", "x", clinical_M),  levels_except_NA(gsub("1[a-c]", "x", clinical_M))),
    # Capsular penetration?
    # Margin status?
    
    # Données de suivi clinique (issues cliniques)
    biochemical.relapse = factor(biochemical_recurrence, levels_except_NA(biochemical_recurrence)),
    years.to.relapse    = as.numeric(days_to_first_biochemical_recurrence) / 365,
    all.metastasis      = factorize.clinical.metastasis(new_neoplasm_event_occurrence_anatomic_site, new_neoplasm_occurrence_anatomic_site_text),
    castrate.resistant  = factorize.castrate.resistant(days_to_second_biochemical_recurrence),
    #castrate.resistant  = factor(tumor_progression_post_ht, levels_except_NA(tumor_progression_post_ht)),
    death               = factorize.death(vital_status, patient_death_reason),
    
    # Données de base sur l'infiltration tumorale
    lymphocyte.infiltration = factorize_percent(percent_lymphocyte_infiltration),
	monocyte.infiltration   = factorize_percent(percent_monocyte_infiltration),
	necrosis                = factorize_percent(percent_necrosis),
    neutrophil.infiltration = factorize_percent(percent_neutrophil_infiltration),
    normal.cells            = factorize_percent(percent_normal_cells),
    stromal.cells           = factorize_percent(percent_stromal_cells),
    tumor.cells	            = factorize_percent(percent_tumor_cells),
    tumor.nuclei            = factorize_percent(percent_tumor_nuclei)
    ))
    
write.table(selectedDF, "output/clinical.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

sink("output/clinical_summary.txt")
for(i in 2:ncol(selectedDF)) {
   cat(paste(colnames(selectedDF)[i], ":\n"))
#   print(summary(selectedDF[,i]))
   if(is.numeric(selectedDF[,i])) {
       cat("\tMean:", mean(selectedDF[,i], na.rm=TRUE), "\n")
       cat("\tMedian:", median(selectedDF[,i], na.rm=TRUE), "\n")
   } else {
       counts <- c(Unknown=sum(is.na(selectedDF[,i])), table(selectedDF[,i]))
       percents <- counts / nrow(selectedDF) * 100
       showDF <- data.frame(Count=counts, Percent=percents, row.names=names(counts))
       print(showDF)
   }
   cat("\n")
}
sink(NULL)

# ggplot(data=selectedDF) + geom_bar(mapping=aes(x=stage, fill=biochemical.relapse, stat="count"))
# ggplot(data=selectedDF) + geom_bar(mapping=aes(x=stage.category, fill=biochemical.relapse, stat="count"))
# ggplot(data=selectedDF) + geom_bar(mapping=aes(x=tumor.type, fill=biochemical.relapse, stat="count"))
# ggplot(data=selectedDF) + geom_bar(mapping=aes(x=limph.node, fill=biochemical.relapse, stat="count"))
# ggplot(data=selectedDF) + geom_bar(mapping=aes(x=metastases, fill=biochemical.relapse, stat="count"))

# Correlate biochemical relapse with other clinical factors.
# For binary factors, perform hypertests.

# For ordered parameters, perform logistic regression.
for(clinical.variable in c("gleason.score")) {
    model.formula = paste("biochemical.relapse ~ as.numeric(", clinical.variable, ")", sep="")
    logit.model = glm(eval(parse(text=model.formula)), data = selectedDF, family = "binomial")
}

# Extremely naive approach: correlate variables one by one, as numeric values.
results <- c()
for(i in 2:ncol(selectedDF)) {
    results <- c(results, cor(as.numeric(selectedDF$biochemical.relapse), as.numeric(selectedDF[,i]), use="pairwise.complete.obs"))
}
names(results) = colnames(selectedDF)[-1]

# Try to see which clinical variable overlap each other
numeric.matrix <- matrix(0, nrow=nrow(selectedDF), ncol=ncol(selectedDF)-1)
for(i in 2:ncol(selectedDF)) {
    numeric.matrix[,i-1] <- as.numeric(selectedDF[,i])
}
colnames(numeric.matrix) = colnames(selectedDF)[-1]
cor.matrix = cor(numeric.matrix, use="pairwise.complete.obs")
plot(hclust(dist(scale(cor.matrix, center=TRUE, scale=TRUE))))

all.predictive = c("age", "race", "tumor.type", "psa.preop", "gleason.score", "gleason.category", "stage", "stage.category", "limph.node", "metastases", "lymphocyte.infiltration", "monocyte.infiltration" , "necrosis", "neutrophil.infiltration", "normal.cells", "stromal.cells", "tumor.cells", "tumor.nuclei")
cor.matrix = cor(numeric.matrix[,colnames(numeric.matrix) %in% all.predictive], use="pairwise.complete.obs")
plot(hclust(dist(scale(cor.matrix, center=TRUE, scale=TRUE))))

predictive.non.redundant = c("age", "race", "tumor.type", "psa.preop", "gleason.score", "stage", "limph.node", "metastases", "lymphocyte.infiltration", "monocyte.infiltration" , "necrosis", "neutrophil.infiltration", "stromal.cells", "tumor.cells")
cor.matrix = cor(numeric.matrix[,colnames(numeric.matrix) %in% predictive.non.redundant], use="pairwise.complete.obs")
plot(hclust(dist(scale(cor.matrix, center=TRUE, scale=TRUE))))

cor.matrix = cor(cbind(numeric.matrix[,colnames(numeric.matrix) %in% predictive.non.redundant], biochemical.relapse=as.numeric(selectedDF$biochemical.relapse)), use="pairwise.complete.obs")
plot(hclust(dist(scale(cor.matrix, center=TRUE, scale=TRUE))))

test.model = glm(biochemical.relapse ~ as.numeric(psa.preop) + as.numeric(gleason.score) + as.numeric(stage) + as.numeric(limph.node), family="binomial", data=selectedDF)
summary(test.model)

# Look at hyper enrichment for all parametric factors

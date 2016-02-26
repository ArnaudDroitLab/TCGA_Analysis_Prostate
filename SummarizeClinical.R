library(TCGAbiolinks)
library(survival)

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
    patient.barcode = bcr_patient_barcode,
    
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
    years.to.death      = as.numeric(days_to_death) / 365,
    
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


# Summarize clinical data in pseudo table form
sink("output/clinical_summary.txt")
for(i in 3:ncol(selectedDF)) {
   #Print variable name
   cat(paste(colnames(selectedDF)[i], ":\n"))

   # Print either the counts or the mean/median, depending on the data type.
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

# Extremely naive approach: correlate variables one by one, as numeric values.
cor.to.relapse <- c()
for(i in 3:ncol(selectedDF)) {
    cor.to.relapse <- c(cor.to.relapse, cor(as.numeric(selectedDF$biochemical.relapse),
                                            as.numeric(selectedDF[,i]),
                                            use="pairwise.complete.obs"))
}
names(cor.to.relapse) = colnames(selectedDF)[c(-1, -2)]
write.table(data.frame(Variable=names(cor.to.relapse), Correlation=cor.to.relapse), file="output/Naive correlation.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


# Try to see which clinical variable overlap each other
# Transform the factor matrix into a numeric matrix.
numeric.matrix <- matrix(0, nrow=nrow(selectedDF), ncol=ncol(selectedDF)-2)
for(i in 3:ncol(selectedDF)) {
    numeric.matrix[,i-2] <- scale(as.numeric(selectedDF[,i]), center=TRUE, scale=TRUE)
}
colnames(numeric.matrix) = colnames(selectedDF)[c(-1, -2)]

# Cluster all variables.
pdf("output/Cluster all clinical.pdf")
# Do not cluster follow.up since it is orthogonal to years.to.death, and is also useless.
plot(hclust(dist(t(numeric.matrix[,colnames(numeric.matrix)!="years.to.death"]))))
dev.off()

# Cluster only predictive, non-redundant variables.
predictive.non.redundant = c("age", "race", "tumor.type", "psa.preop", "gleason.score", "stage", "limph.node", "metastases", "lymphocyte.infiltration", "monocyte.infiltration" , "necrosis", "neutrophil.infiltration", "stromal.cells", "tumor.cells")
num2.matrix = cbind(numeric.matrix[,colnames(numeric.matrix) %in% predictive.non.redundant], biochemical.relapse=as.numeric(selectedDF$biochemical.relapse))
pdf("output/Cluster predictive non redundant clinical.pdf")
plot(hclust(dist(t(num2.matrix))))
dev.off()

# Build a logistic regression model
test.model = glm(biochemical.relapse ~ as.numeric(psa.preop) + as.numeric(gleason.score) + as.numeric(stage) + limph.node + tumor.type, family="binomial", data=selectedDF)
sink("output/logistic regression of relapse on clinical variables.txt")
print(summary(test.model))
sink(NULL)

# Do Cox proportional hazard regression.
# Code time to event and censoring.
selectedDF$time.to.event <- with(selectedDF, pmin(years.to.relapse, follow.up, years.to.death, na.rm=TRUE))
selectedDF$event.type <- with(selectedDF, ifelse(!is.na(years.to.relapse), 1, 0))

# Fit model
cox.fit = coxph(Surv(time.to.event, event.type) ~ as.numeric(gleason.category) + as.numeric(stage.category) + limph.node + as.numeric(psa.preop) + tumor.type + age, data=selectedDF)
sink("output/Cox regression of clinical variables.txt")
print(summary(cox.fit))
sink(NULL)

# Plot average survival time
pdf(width=7, height=7, file="output/Overall survival.pdf")
plot(survfit(cox.fit), ylim=c(0.5, 1), xlab="Years", ylab="Proportion without recurrence")
dev.off()

# Plot effects of gleason category on time to recurrence.
selectedDF.gleason <- with(selectedDF, data.frame(gleason.category = c(0, 1, 2, 3),
                                                  stage.category   = rep(mean(as.numeric(stage.category),      na.rm=TRUE), 4),
                                                  limph.node = rep(mean(as.numeric(limph.node), na.rm=TRUE), 4),
                                                  psa.preop  = rep(mean(as.numeric(psa.preop),  na.rm=TRUE), 4),
                                                  tumor.type = rep(mean(as.numeric(tumor.type), na.rm=TRUE), 4),
                                                  age        = rep(mean(age, na.rm=TRUE), 4)))
                                                  
pdf(width=7, height=7, file="output/gleason survival.pdf")
plot(survfit(cox.fit, newdata=selectedDF.gleason), conf.int=FALSE, mark.time=FALSE,
    lty=c(1, 2, 3, 4), ylim=c(0.4, 1), xlab="Years",
    ylab="Proportion without recurrence")
legend("bottomleft", legend=levels(selectedDF$gleason.category), lty=c(1 ,2, 3, 4), inset=0.02)
dev.off()

# Plot effects of gleason score on stage.category to recurrence.
selectedDF.stage <- with(selectedDF, data.frame(stage.category = c(0, 1, 2),
                                                gleason.category = rep(mean(as.numeric(gleason.category), na.rm=TRUE), 3),
                                                limph.node       = rep(mean(as.numeric(limph.node),       na.rm=TRUE), 3),
                                                psa.preop        = rep(mean(as.numeric(psa.preop),        na.rm=TRUE), 3),
                                                tumor.type       = rep(mean(as.numeric(tumor.type),       na.rm=TRUE), 3),
                                                age              = rep(mean(age, na.rm=TRUE), 3)))
                                                  
pdf(width=7, height=7, file="output/stage survival.pdf")
plot(survfit(cox.fit, newdata=selectedDF.stage), conf.int=FALSE, mark.time=FALSE,
    lty=c(1, 2, 3),  ylim=c(0.3, 1), xlab="Years",
    ylab="Proportion without recurrence")
legend("bottomleft", legend=levels(selectedDF$stage.category), lty=c(1 ,2, 3), inset=0.02)
dev.off()

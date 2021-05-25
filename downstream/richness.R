library(phyloseq)
library(ggplot2)

metadata <- read.csv("./table.csv")
data <- import_biom("./abundance.biom", parseFunction=parse_taxonomy_default)

metadata <- metadata[!duplicated(metadata),]
rownames(metadata) <- NULL

r <- estimate_richness(data)
metadata$InvSimpson <- r[match(metadata$sample_name, rownames(r)), ]$InvSimpson


patients <- metadata[c("patient_id", "treatment")]
patients <- patients[!duplicated(patients),]
rownames(patients) <- NULL

patients$richness_before <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T1",]$InvSimpson, patients$patient_id, patients$treatment))
patients$richness_after  <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T2",]$InvSimpson, patients$patient_id, patients$treatment))

patients <- patients[complete.cases(patients),]
rownames(patients) <- NULL

write.csv(patients, "richness.csv", row.names=FALSE)

plot_before <- ggplot(patients, aes(treatment, richness_before)) + geom_boxplot()
plot_after  <- ggplot(patients, aes(treatment, richness_after)) + geom_boxplot()

ggsave("before.png", plot=plot_before, width=8, height=5)
ggsave("after.png", plot=plot_after, width=8, height=5)

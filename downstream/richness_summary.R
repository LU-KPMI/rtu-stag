library(ggplot2)
library(gridExtra)

metadata <- read.csv("./sample_metadata.csv")
richness <- read.csv("./richness.csv")

metadata$chao1 <- richness[match(metadata$sample_name, richness$sample_name), ]$chao1
metadata$shannon <- richness[match(metadata$sample_name, richness$sample_name), ]$shannon
metadata$observed <- richness[match(metadata$sample_name, richness$sample_name), ]$observed


patients <- metadata[c("patient_id", "treatment")]
patients <- patients[!duplicated(patients),]
rownames(patients) <- NULL

patients$chao1_before <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T1",]$chao1, patients$patient_id, patients$treatment))
patients$chao1_after  <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T2",]$chao1, patients$patient_id, patients$treatment))

patients$shannon_before <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T1",]$shannon, patients$patient_id, patients$treatment))
patients$shannon_after  <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T2",]$shannon, patients$patient_id, patients$treatment))

patients$observed_before  <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T1",]$observed, patients$patient_id, patients$treatment))
patients$observed_after  <- as.numeric(mapply(function(patient_id, treatment) metadata[metadata$patient_id == patient_id & metadata$treatment == treatment & metadata$time == "T2",]$observed, patients$patient_id, patients$treatment))

patients <- patients[complete.cases(patients),]
rownames(patients) <- NULL

write.csv(patients, "patient_richness.csv", row.names=FALSE, quote = FALSE)

chao1_plot_before <- ggplot(patients, aes(treatment, chao1_before)) + geom_violin() + labs(x = "", y = "Chao1") + ggtitle("Before") + coord_cartesian(ylim = c(0, 1500))
chao1_plot_after  <- ggplot(patients, aes(treatment, chao1_after)) + geom_violin() + labs(x = "", y = "") + ggtitle("After") + coord_cartesian(ylim = c(0, 1500))

shannon_plot_before <- ggplot(patients, aes(treatment, shannon_before)) + geom_violin() + labs(x = "", y = "Shannon") + ggtitle("Before") + coord_cartesian(ylim = c(0, 5))
shannon_plot_after  <- ggplot(patients, aes(treatment, shannon_after)) + geom_violin() + labs(x = "", y = "") + ggtitle("After") + coord_cartesian(ylim = c(0, 5))

observed_plot_before <- ggplot(patients, aes(treatment, observed_before)) + geom_violin() + labs(x = "", y = "Observed") + ggtitle("Before") + coord_cartesian(ylim = c(0, 1500))
observed_plot_after  <- ggplot(patients, aes(treatment, observed_after)) + geom_violin() + labs(x = "", y = "") + ggtitle("After") + coord_cartesian(ylim = c(0, 1500))

ggsave("richness.png", grid.arrange(chao1_plot_before, chao1_plot_after, shannon_plot_before, shannon_plot_after, observed_plot_before, observed_plot_after, ncol=2, nrow=3), width=8, height=8)

library(phyloseq)
library(ggplot2)

data <- import_biom("./abundance.biom", parseFunction=parse_taxonomy_default)

sample_names(data)
rank_names(data)

write.csv(estimate_richness(data), file="alpha-diversity.csv")
write.csv(as.matrix(phyloseq::distance(data, method="bray")), file = "beta-diversity.csv")

total = median(sample_sums(data))
standf = function(x, t=total) round(t * (x / sum(x)))
data = transform_sample_counts(data, standf)

plot_bar(data, fill = "Rank2") +
  geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack")
ggsave("abundance.png", width = 8, height = 5)

plot_bar(data, x="Rank2", fill="Rank2") +
  geom_bar(aes(color=Rank2, fill=Rank2), stat="identity", position="stack")
ggsave("total_abundance.png", width = 8, height = 5)

plot_heatmap(data, method = "NMDS", distance = "bray")
ggsave("heatmap.png", width = 8, height = 5)

plot_richness(data, measures=c("Shannon", "Simpson"))
ggsave("richness.png", width = 8, height = 5)

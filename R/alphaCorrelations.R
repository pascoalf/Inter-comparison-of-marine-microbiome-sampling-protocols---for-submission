## supplementary Figure for alpha correlations

# Prokaryotes 16S
grid.arrange(
  prok_16S_diversity_rarefied_metadata %>%
    ggplot(aes(x = Species_richness, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Shannon index",
         title = "MetaB16SV4V5 (Prokaryotes)",
         subtitle = paste("Pearson correlation", 
                       round(cor(prok_16S_diversity_rarefied_metadata$Species_richness,
                                 prok_16S_diversity_rarefied_metadata$Shannon_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  prok_16S_diversity_rarefied_metadata %>%
    ggplot(aes(x = Species_richness, y = Simpson_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Simpson index",
         title = "MetaB16SV4V5 (Prokaryotes)",
         subtitle = paste("Pearson correlation", 
                       round(cor(prok_16S_diversity_rarefied_metadata$Species_richness,
                                 prok_16S_diversity_rarefied_metadata$Simpson_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  prok_16S_diversity_rarefied_metadata %>%
    ggplot(aes(x = Simpson_index, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Simpson index",
         y = "Shannon index",
         title = "MetaB16SV4V5 (Prokaryotes)",
         subtitle = paste("Pearson correlation", 
                       round(cor(prok_16S_diversity_rarefied_metadata$Simpson_index,
                                 prok_16S_diversity_rarefied_metadata$Shannon_index),
                             digits = 3))) +
    theme(panel.grid = element_blank())
#)
,
## Prokaryotes Metagenomes
#grid.arrange(
  prokaryotic_metagenome_diversity %>%
    ggplot(aes(x = Species_richness, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Shannon index",
         title = "MetaG (Prokaryotes)",
         subtitle = paste("Pearson correlation", 
                       round(cor(prokaryotic_metagenome_diversity$Species_richness,
                                 prokaryotic_metagenome_diversity$Shannon_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  prokaryotic_metagenome_diversity %>%
    ggplot(aes(x = Species_richness, y = Simpson_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Simpson index",
         title = "MetaG (Prokaryotes)",
         subtitle = paste("Pearson correlation", 
                       round(cor(prokaryotic_metagenome_diversity$Species_richness,
                                 prokaryotic_metagenome_diversity$Simpson_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  prokaryotic_metagenome_diversity %>%
    ggplot(aes(x = Simpson_index, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Simpson index",
         y = "Shannon index",
         title = "MetaG (Prokaryotes)",
         subtitle = paste("Pearson correlation", 
                       round(cor(prokaryotic_metagenome_diversity$Simpson_index,
                                 prokaryotic_metagenome_diversity$Shannon_index),
                             digits = 3))) +
    theme(panel.grid = element_blank())
#)
,
## Protists 18S
#grid.arrange(
  protist_18S_diversity %>%
    ggplot(aes(x = Species_richness, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Shannon index",
         title = "MetaB18SV9 (Protists)",
         subtitle = paste("Pearson correlation", 
                       round(cor(protist_18S_diversity$Species_richness,
                                 protist_18S_diversity$Shannon_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  protist_18S_diversity %>%
    ggplot(aes(x = Species_richness, y = Simpson_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Simpson index",
         title = "MetaB18SV9 (Protists)",
         subtitle = paste("Pearson correlation", 
                       round(cor(protist_18S_diversity$Species_richness,
                                 protist_18S_diversity$Simpson_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  protist_18S_diversity %>%
    ggplot(aes(x = Simpson_index, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Simpson index",
         y = "Shannon index",
         title = "MetaB18SV9 (Protists)",
         subtitle = paste("Pearson correlation", 
                       round(cor(protist_18S_diversity$Simpson_index,
                                 protist_18S_diversity$Shannon_index),
                             digits = 3))) +
    theme(panel.grid = element_blank())
#)
,
## Protists metagenomes
#grid.arrange(
  protist_metagenome_diversity %>%
    ggplot(aes(x = Species_richness, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Shannon index",
         title = "MetaG (Protists)",
         subtitle = paste("Pearson correlation", 
                       round(cor(protist_metagenome_diversity$Species_richness,
                                 protist_metagenome_diversity$Shannon_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  protist_metagenome_diversity %>%
    ggplot(aes(x = Species_richness, y = Simpson_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Species richness",
         y = "Simpson index",
         title = "MetaG (Protists)",
         subtitle = paste("Pearson correlation", 
                       round(cor(protist_metagenome_diversity$Species_richness,
                                 protist_metagenome_diversity$Simpson_index),
                             digits = 3)))+
    theme(panel.grid = element_blank()),
  protist_metagenome_diversity %>%
    ggplot(aes(x = Simpson_index, y = Shannon_index)) + 
    geom_point() + 
    geom_smooth(se = FALSE, method = "lm", lty = "dashed", col = "grey") +
    labs(x = "Simpson index",
         y = "Shannon index",
         title = "MetaG (Protists)",
         subtitle = paste("Pearson correlation", 
                       round(cor(protist_metagenome_diversity$Simpson_index,
                                 protist_metagenome_diversity$Shannon_index),
                             digits = 3))) +
    theme(panel.grid = element_blank()), ncol = 3
)

##
qualitative_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

## Volumes from 1L to 100L
# Prokaryotes
# Number of OTU
prok_diversity_metadata %>%
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex")) %>% 
  filter(effected_volume <= 100) %>%
  ggplot(aes(effected_volume,Species_richness,col=device))+
  geom_jitter(width = 0.1, size = 3, alpha = 0.75)+
  facet_grid(facets = c("Sequencing_strategy","size_fraction"),scales = "free")+
  labs(x="Volume (L)",
       y="Number of OTU",
       col = "Filter type",
       title = "OTUS | Volumes from 1L to 100L | Prokaryotes | 16S and MetaG | Size fractions and Whole water")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  scale_color_manual(values = qualitative_colors[c(3,5)])

# Protists
# Number of OTU
prot_diversity_metadata %>% 
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex")) %>%
  filter(!is.na(size_fraction),!is.na(Sequencing_strategy),
         effected_volume <= 100) %>% 
  ggplot(aes(effected_volume,Species_richness,col=device))+
  geom_jitter(width = 0.1, size = 3, alpha = 0.75)+
  facet_grid(facets = c("Sequencing_strategy","size_fraction"),scales = "free")+
  labs(x="Volume (L)",
       y="Number of OTU",
       col="Filter",
       title = "OTUS | All volumes| Protists | 18S and MetaG | Size fractions and Whole water")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  scale_color_manual(values = qualitative_colors[c(3,5)])

## Volumes from 1L to 1000L
# Number of OTU
prok_diversity_metadata %>% 
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex")) %>%
  filter(!is.na(size_fraction),!is.na(Sequencing_strategy)) %>% 
  ggplot(aes(effected_volume,Species_richness,col=device))+
  geom_jitter(width = 0.1, size = 3, alpha = 0.75)+  
  facet_grid(facets = c("Sequencing_strategy","size_fraction"),scales = "free")+
  labs(x="Volume (L)",
       y="Number of OTU",
       colour = "Filter",
       title = "OTUS | All Volumes | Prokaryotes | 16S and MetaG | Size fractions and Whole water")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  scale_color_manual(values = qualitative_colors[c(3,5)])

# Protists
# Number of OTU
prot_diversity_metadata %>% 
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex")) %>%
  filter(!is.na(size_fraction),!is.na(Sequencing_strategy)) %>% 
  ggplot(aes(effected_volume,Species_richness,col=device))+
  geom_jitter(width = 0.1, size = 3, alpha = 0.75)+
  facet_grid(facets = c("Sequencing_strategy","size_fraction"),scales = "free")+
  labs(x="Volume (L)",
       y="Number of OTU",
       colour="Filter type",
       title = "OTUS | All Volumes | Protists | 18S and MetaG | Size fractions and Whole water")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  scale_color_manual(values = qualitative_colors[c(3,5)])

## Filter type (10L, Whole water), sterivex vs membrane
# Prokaryotes
# OTUs
prokaryotes_membrane_vs_sterivex_10L_whole_water <-
  prok_diversity_metadata %>% 
  filter(method=="filter",effected_volume == 10) %>% 
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex"))

prok_wilc_dev_10L_ww <-
  prokaryotes_membrane_vs_sterivex_10L_whole_water %>% 
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ device) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x = "device") %>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotic")

# Plot
prokaryotes_membrane_vs_sterivex_10L_whole_water %>% 
  ggplot(aes(x = device,
             y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = device), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Number of OTU", col = NULL,
       title = "OTUs | Filter type | Prokaryotes | 16S and MetaG")+
  facet_wrap(facets = c("Sequencing_strategy"),scales= "free")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)])
  
# quick statistics
prokaryotes_membrane_vs_sterivex_10L_whole_water %>% 
  group_by(Sequencing_strategy,device) %>% 
  summarise(min(Species_richness),
            max(Species_richness),
            median(Species_richness),
            IQR(Species_richness),
            n())

## Sterivex vs membrane, protists
# OTUS
protist_membrane_vs_sterivex_10L_whole_water <-
  prot_diversity_metadata %>% 
  filter(method == "filter", effected_volume == 10) %>% 
  mutate(device = ifelse(device == "membrane", "Membrane", "Sterivex"))

prot_wilc_m_vs_st_10L_ww <-
  protist_membrane_vs_sterivex_10L_whole_water %>% 
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ device) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x = "device")%>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Protist") 

# Plot (Sterivex vs membrane, 10L, WW, protist)
protist_membrane_vs_sterivex_10L_whole_water %>% 
  ggplot(aes(x = device,
             y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = device), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | Filter type | Protists | 18S and MetaG")+
facet_wrap(facets = c("Sequencing_strategy"), scale= "free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)])

# quick summary
prot_diversity_metadata %>% 
  filter(method=="filter",effected_volume == 10) %>% 
  group_by(Sequencing_strategy,device) %>% 
  summarise(min(Species_richness),
            max(Species_richness), 
            median(Species_richness),
            IQR(Species_richness),
            n())

## Whole water vs size fractions, 10L
# OTUs
prokaryotes_ww_vs_sf_10L_membrane <-
  prok_diversity_metadata %>% 
  filter(effected_volume == 10,device == "membrane")

# Kruskal test
prokaryotes_kruskal_test_ww_vs_sf_10L_membrane <-
  prokaryotes_ww_vs_sf_10L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance()%>% 
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Kruskal-Wallis",
         Taxonomic_group = "Prokaryotes")

prokaryotes_ww_vs_sf_10L_membrane %>% 
  group_by(Sequencing_strategy) %>% 
  kruskal_effsize(Species_richness ~ size_fraction)

# post-hoc test (Dunn)
prokaryotes_dunn_ww_vs_sf_10L_membrane <- 
  prokaryotes_ww_vs_sf_10L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  dunn_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x="size_fraction") %>% 
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Dunn",
         Taxonomic_group = "Prokaryotes")

prokaryotes_dunn_ww_vs_sf_10L_membrane$xmax <- c(2,3,3,2)

## Wilcoxon test for metagenomes
prokaryotes_wilcox_test_ww_vs_sf_10L_membrane <-
  prokaryotes_ww_vs_sf_10L_membrane %>% 
  filter(Sequencing_strategy != "MetaB16SV4V5") %>% 
  wilcox_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(Sequencing_strategy = "Metagenome",
         Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotes")

# Plot
prokaryotes_ww_vs_sf_10L_membrane %>% 
  ggplot(aes(x = size_fraction,
             y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | WW vs size fraction | 10L | Prokaryotes | 16S and MetaG",
       subtitle = get_test_label(prokaryotes_kruskal_test_ww_vs_sf_10L_membrane, type = "text", detailed = TRUE),
       caption = get_pwc_label(prokaryotes_dunn_ww_vs_sf_10L_membrane, type = "expression"))+
  facet_wrap(facets = c("Sequencing_strategy"),scale="free")+
  stat_pvalue_manual(prokaryotes_dunn_ww_vs_sf_10L_membrane,
                     label = "p.adj.signif",
                     y.position = 550,
                     hide.ns = TRUE, 
                     size = 10) + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)])
  
# Protists, whole water vs size fractions, 10L
protist_ww_sf_10L_membrane <- 
  prot_diversity_metadata %>% 
  filter(device == "membrane",effected_volume == 10)

# Kruskall wallis test, protist, ww vs sf, 10L
protist_ww_sf_10L_membrane_kruskall_test <- 
  protist_ww_sf_10L_membrane %>% 
  group_by(Sequencing_strategy) %>%
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(Test = "Kruskal-Wallis",
         Taxonomic_group = "Protist")

# Dunn test not necessary, because Kruskall wallis was not significant.  

#
protist_ww_sf_10L_membrane %>% 
  ggplot(aes(x = size_fraction, y = Species_richness)) +
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | WW vs size fraction | 10L | Protists | 18S and MetaG") + #,
  facet_wrap(facets = c("Sequencing_strategy"), scale= "free") +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)])

# quick overview, protists, 10L, ww vs sf
protist_ww_sf_10L_membrane %>% 
  group_by(Sequencing_strategy,size_fraction) %>% 
  summarise(min(Species_richness),
            max(Species_richness))

## Size fractions, 100L, membrane, prokaryotes
prokaryotes_size_fractions_100L_membrane <- 
  prok_diversity_metadata %>% 
  filter(device == "membrane",effected_volume == 100)

## 16S (3 groups to compare)
# Kruskal test, size fractions, membrane prokaryotes, (16S only)
prokaryotes_size_fractions_100L_membrane_kruskal_test_16S <-
  prokaryotes_size_fractions_100L_membrane %>%
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Kruskal-Wallis",
         Taxonomic_group = "Prokaryotic")

# post-hoc test (Dunn) size fractions 100L prokaryotes, (16S only)
prokaryotes_size_fractions_100L_membrane_dunn_test_16S <- 
  prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  dunn_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x="size_fraction") %>%
  mutate(Sequencing_strategy = "MetaB16SV4V5",
         Test = "Dunn",
         Taxonomic_group = "Prokaryotic")

prokaryotes_size_fractions_100L_membrane_dunn_test_16S$xmin <- c(1,1,2)
prokaryotes_size_fractions_100L_membrane_dunn_test_16S$xmax <- c(2,3,3)

# Size fractions, 100L, membrane, 16S only
prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  ggplot(aes(x = size_fraction, y = Species_richness))+
  geom_boxplot(outlier.shape = "cross", size = 1)+
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 4, alpha =0.75)+
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | Size fraction | 100L | Prokaryotes | 16S",
       subtitle = get_test_label(prokaryotes_size_fractions_100L_membrane_kruskal_test_16S, type = "text", detailed = TRUE),
       caption = get_pwc_label(prokaryotes_size_fractions_100L_membrane_dunn_test_16S, type = "expression"))+
  stat_pvalue_manual(prokaryotes_size_fractions_100L_membrane_dunn_test_16S,
                     label = "p.adj.signif",
                     y.position = 650,
                     hide.ns = TRUE, 
                     size = 10) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)])

## Some quick metrics
prokaryotes_size_fractions_100L_membrane %>% 
  group_by(Sequencing_strategy,size_fraction) %>% 
  summarize(median(Species_richness))

## metagenomes (2 groups to compare)
# Size fractions, 100L, membrane, Metagenomes
prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes <-
  prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaG") %>%
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(x = "size_fraction") %>% 
  mutate(Sequencing_strategy = "Metagenome",
         Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotic")

prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes$xmin <- 1
prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes$xmax <- 2
# Plot
prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaG") %>% 
  ggplot(aes(x = size_fraction, y = Species_richness))+
  geom_boxplot(outlier.shape = "cross", size = 1)+
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | size fractions | Prokaryotes | Metagenome",
       subtitle = get_test_label(prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes, type = "text", detailed = TRUE))+
  stat_pvalue_manual(prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes,
                     label = "p.adj.signif",
                     y.position = 230,
                     hide.ns = FALSE, 
                     size = 10) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)])

## Size fractions, 100L, membrane, protists
protist_size_fractions_100L_membrane <- 
  prot_diversity_metadata %>% 
  filter(device == "membrane",effected_volume == 100)

## Kruskall wallis test
protist_size_fractions_100L_membrane_kruskall <- 
  protist_size_fractions_100L_membrane %>% 
  group_by(Sequencing_strategy) %>% 
  kruskal_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(test = "Kruskal-Wallis",
         Taxonomic_group = "Protist")
  
# post-hoc test for 18S (Dune test, protist, size fractions, 100L, membrane)
protist_size_fractions_100L_membrane_dunn <-
  protist_size_fractions_100L_membrane %>% 
  group_by(Sequencing_strategy) %>% 
  dunn_test(Species_richness ~ size_fraction) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  add_xy_position(group = "size_fraction") %>% 
  mutate(test = "Dunn",
         Taxonomic_group = "Protist")

# Plot
protist_size_fractions_100L_membrane %>% 
  ggplot(aes(x = size_fraction, y = Species_richness)) +
  geom_boxplot(outlier.shape = "cross", size = 1) +
  geom_jitter(aes(color = size_fraction), width = 0.1, size = 4, alpha = 0.75) +
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | Size fraction | 100L | Protists | 18S and MetaG") + #,
  facet_wrap(facets = c("Sequencing_strategy"), scale="free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)])

## 2.5L vs 10L, WW, sterivex
# Prokaryotes
prok_2.5L_vs_10L_wilcoxon <- 
  prok_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L") %>%
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ effected_volume) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Prokaryotic")
  
# plot
prok_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L") %>% 
  ggplot(aes(x = as.factor(effected_volume),
             y = Species_richness))+
  geom_boxplot(outlier.shape = "cross", size = 1)+
  geom_jitter(aes(col = as.factor(effected_volume)), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | 2.5L pooled vs 10L | Sterivex | Prokaryotes | 16S and MetaG")+
  facet_wrap(facets = c("Sequencing_strategy"),scale="free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)])

#quick overview
prok_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L") %>% 
  group_by(Sequencing_strategy,effected_volume) %>% 
  summarise(min(Species_richness),
            max(Species_richness),
            median(Species_richness),
            IQR(Species_richness),
            n()) %>% View()

# Protist 2.5L vs 4*2.5L, sterivex, WW
protist_2.5_vs_10_st_ww <- 
  prot_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L")

# Wilcoxon test
prot_2.5_vs_10_st_ww_wilcox <- 
  protist_2.5_vs_10_st_ww %>% 
  group_by(Sequencing_strategy) %>% 
  wilcox_test(Species_richness ~ planned_volume) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  mutate(Test = "Mann-Whitney",
         Taxonomic_group = "Protist")

protist_2.5_vs_10_st_ww %>% 
  group_by(Sequencing_strategy) %>% 
  wilcox_effsize(Species_richness ~ planned_volume)

# OTUs
protist_2.5_vs_10_st_ww %>% 
  ggplot(aes(x = as.factor(planned_volume), y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = as.factor(planned_volume)), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Number of OTU",
       title = "OTUs | 2.5L pooled vs 10L | Sterivex | Protists | 18S and MetaG") + #,
 facet_wrap(facets = c("Sequencing_strategy"),scale="free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)])
  


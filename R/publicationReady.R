# Publication ready


## Figure 4 --> (Figure 2 revised)
prok_diversity_metadata %>% 
  mutate(device = ifelse(device == "sterivex", "Cartridge membrane filter", device),
         device = ifelse(device == "membrane", "Flat membrane filter", device)) %>%
  filter(!is.na(size_fraction),!is.na(Sequencing_strategy)) %>% 
  ggplot(aes(effected_volume,Species_richness,col=device))+
  geom_jitter(width = 0.1, size = 3, alpha = 0.75)+  
  facet_grid(facets = c("Sequencing_strategy","size_fraction"),scales = "free_y")+
  labs(x="Volume (L) in Log10 scale",
       y="Species richness in Log10 scale",
       colour = "Filter",
       title = "")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"),
        legend.position = "top")+
  scale_color_manual(values = qualitative_colors[c(3,5)]) +
  scale_x_log10()+
  scale_y_log10()


## Figure 5 --> Figure 3 (revised)
# A
figure_3_a <- 
  prokaryotes_ww_vs_sf_10L_membrane %>% 
  ggplot(aes(x = size_fraction,
             y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 2, alpha = 0.75)+
  labs(x="Filter pore size",
       y="Species richness",
       #shape = "Replicate:",
       title = "a",
       #subtitle = get_test_label(prokaryotes_kruskal_test_ww_vs_sf_10L_membrane, type = "text", detailed = TRUE),
       #caption = get_pwc_label(prokaryotes_dunn_ww_vs_sf_10L_membrane, type = "expression")
       )+
  facet_wrap(facets = c("Sequencing_strategy"),scale="free")+
  stat_pvalue_manual(prokaryotes_dunn_ww_vs_sf_10L_membrane,
                     label = "p.adj.signif",
                     y.position = 540,
                     hide.ns = TRUE, 
                     size = 7) + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 10, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)])

# B left
figure_3_b_left <- 
  prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5") %>% 
  ggplot(aes(x = size_fraction, y = Species_richness))+
  geom_boxplot(outlier.shape = "cross", size = 1) +
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 2, alpha =0.75)+
  labs(x="Filter pore size",
       y="Species richness",
      # shape = "Replicate:",
       title = "b",
       #subtitle = paste( 
        # get_test_label(prokaryotes_size_fractions_100L_membrane_kruskal_test_16S, type = "text", detailed = TRUE),
         #"\n MetaB16SV4V5"),
       #caption = get_pwc_label(prokaryotes_size_fractions_100L_membrane_dunn_test_16S, type = "expression")
      )+
  stat_pvalue_manual(prokaryotes_size_fractions_100L_membrane_dunn_test_16S,
                     label = "p.adj.signif",
                     y.position = 640,
                     hide.ns = TRUE, 
                     size = 7) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 10, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)])+
  facet_wrap(~Sequencing_strategy)

# B right
figure_3_b_right <- 
  prokaryotes_size_fractions_100L_membrane %>% 
  filter(Sequencing_strategy == "MetaG") %>% 
  ggplot(aes(x = size_fraction, y = Species_richness))+
  geom_boxplot(outlier.shape = "cross", size = 1)+
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 2, alpha = 0.75)+
  labs(x="Filter pore size",
       y="Species richness",
       #shape = "Replicate:",
       title = ""
      # subtitle = paste(
       #  get_test_label(prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes, type = "text", detailed = TRUE),
        # "\n MetaG")
      )+
  stat_pvalue_manual(mutate(prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes, Sequencing_strategy = "MetaG"),
                     label = "p.adj.signif",
                     y.position = 220,
                     hide.ns = FALSE, 
                     size = 7) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 10, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)]) +
  facet_wrap(~Sequencing_strategy)

# C
figure_3_c <- 
  prokaryotes_membrane_vs_sterivex_10L_whole_water %>% 
  mutate(device = ifelse(device == "Sterivex", "Cartridge \nmembrane \nfilter", device),
         device = ifelse(device == "Membrane", "Flat \nmembrane \nfilter", device)) %>% 
  ggplot(aes(x = device,
             y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = device), width = 0.1, size = 2, alpha = 0.75)+
  labs(x="Filter type",
       y="Species richness",
       #shape = "Replicate:",
       col = NULL,
       title = "c")+
  facet_wrap(facets = c("Sequencing_strategy"),scales= "free")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 10, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)])

# D
figure_3_d <- 
  prok_diversity_metadata %>% 
  filter(device == "sterivex", planned_volume != "1L") %>% 
  ggplot(aes(x = as.factor(effected_volume),
             y = Species_richness))+
  geom_boxplot(outlier.shape = "cross", size = 1)+
  geom_jitter(aes(col = as.factor(effected_volume)), 
              width = 0.1, size = 2, alpha = 0.75)+
  #geom_text_repel(aes(label = replica), size = 2.5)+
  labs(x="Volume (L)",
       y="Species richness",
       title = "d")+
  facet_wrap(facets = c("Sequencing_strategy"),scale="free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 10, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)])


cowplot::plot_grid(figure_3_a + theme(legend.position = "top", 
                                      axis.text.x = element_text(angle = 90)),
                   figure_3_b_left + theme(legend.position = "top"),
                   figure_3_b_right + theme(legend.position = "top"),
                   figure_3_c + theme(legend.position = "top"),
                   figure_3_d + theme(legend.position = "top"),
                   align = "hv",
                   axis = "tblr",
                   #rel_widths = c(1,0.5,0.5,1),
                   label_size = 12)

## Figure 6 --> Figure 4
### Prokaryotes MetaB16S
bd_prok_16S_clean_rarefied <- betadisper(
  vegdist(t(prokaryotes_16S_clean_rarefied)),
  group = pull(select(env_16S,all_of("size_fraction"))),
  type = "centroid")

### Prokaryotes MetaG

bd_prok_meta_clean_rarefied <- betadisper(
  vegdist(t(prokaryotes_metagenome_clean_rarefied)),
  group = pull(select(env_metagenome_prok,all_of("size_fraction"))),
  type = "centroid")

## vector with colors
vol_cols <- distinct(env_metagenome_prok, effected_volume, col_volume) %>% arrange(effected_volume) %>% pull(col_volume)


## Combine
par(mfrow = c(2,2))
# A
plot(mds_prok_16S,
     display = "sites", type="p", 
     main = "a \n Prokaryotes, MetaB16SV4V5")
points(mds_prok_16S,
       display = "sites",
       bg = env_16S$col_volume,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_16S,
     ordiellipse(mds_prok_16S,
                 size_fraction,kind="se",conf=0.95,col="grey"))
with(env_16S,
     ordihull(mds_prok_16S,
              size_fraction,col="grey"))
with(env_16S,
     ordispider(mds_prok_16S,
                size_fraction,label = T,col="grey"))


# B
plot(mds_prok_metagenome,
     display = "sites", type="p", 
     main = "b \n Prokaryotes, MetaG")
points(mds_prok_metagenome,
       display = "sites",
       bg = env_metagenome_prok$col_volume,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_metagenome_prok,
     ordiellipse(mds_prok_metagenome,
                 size_fraction,kind="se",conf=0.95,col="grey"))
with(env_metagenome_prok,
     ordihull(mds_prok_metagenome,
              size_fraction,col="grey"))
with(env_metagenome_prok,
     ordispider(mds_prok_metagenome,
                size_fraction,label = T,col="grey"))

# C
boxplot(bd_prok_16S_clean_rarefied, main = "c \n Prokaryotes, MetaB16SV4V5")

# D
boxplot(bd_prok_meta_clean_rarefied, main = "d \n Prokaryotes, MetaG")
#
dev.off()
plot.new() ### to add color legend in image editor
with(env_metagenome_prok, 
     legend("center", 
            legend = levels(as.factor(effected_volume)),
            bty = "n",
            pt.bg = vol_cols,
            pch = 21, 
            col = "black", 
            title = "Volume (L)"))


#####

# Figure 7 --> Figure 5 revised version
tidy_prok_16S_clean_rarefied_taxa %>%
  group_by(FAKE_rank) %>% 
  complete(size_fraction,
           effected_volume,
           fill = list(Species_richness = 0)) %>% 
  filter(!is.na(effected_volume),
         !is.na(size_fraction)) %>%
  ggplot(aes(x = as.factor(effected_volume), 
             fill = size_fraction, y = Species_richness))+
  geom_col(position = position_dodge(preserve = "total"))+
  facet_wrap(~FAKE_rank, scales = "free")+
  scale_fill_manual(values = qualitative_colors[c(1,3,5,7)])+
  labs(fill = "Pore size (\U00B5m):",
       x = "Volume (L)", 
       y = "Species richness",
       title = "")+
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        axis.text = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 9),
        text = element_text(size = 12))

# Figure 8 -> figure 6 revised
grid.arrange(
  tidy_prok_16S_clean_rarefied_taxa %>% 
    filter(!is.na(effected_volume),!is.na(size_fraction),
           FAKE_rank == "Candidatus Marinimicrobia",
           effected_volume %in% c(100)) %>% 
    ggplot(aes(size_fraction, Species_richness))+
    geom_col(aes(fill = size_fraction, col = size_fraction))+
    labs(y = "Species richness", 
         x = "",
         title = "a") + 
    scale_fill_manual(values = qualitative_colors[c(3,5,7)]) + 
    scale_color_manual(values = qualitative_colors[c(3,5,7)]) + 
    guides(col = "none", fill = "none") + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank())
  ,
  ## Relative abundance analysis of candidate marinimicrobia
  prok_16S_clean_rarefied_taxa %>% 
    filter(FAKE_rank == "Candidatus Marinimicrobia") %>% 
    gather(value = "Reads",key = "Sample",-c(OTU,FAKE_rank)) %>% 
    filter(Reads > 0) %>% 
    mutate(Relative_abundance = Reads*100/sum(Reads)) %>% 
    left_join(metadata, by = c("Sample" = "run_accession")) %>% 
    filter(!is.na(effected_volume),!is.na(size_fraction),
           effected_volume %in% c(100)) %>% 
    ggplot(aes(size_fraction,
               Relative_abundance))+
    geom_boxplot(col="black", size= 1)+
    geom_jitter(aes(col = size_fraction), size = 3)+
    labs(x = "Pore size",
         y = "Relative abundance (%)",
         title = "b")+
    scale_y_continuous(labels=label_number(accuracy = 1))+
    scale_y_log10() + 
    theme(axis.text = element_text(family = "Helvetica"),
          text = element_text(family = "Helvetica")) + 
    scale_fill_manual(values = qualitative_colors[c(3,5,7)]) + 
    scale_color_manual(values = qualitative_colors[c(3,5,7)]) + 
    guides(col = "none", fill = "none") + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank())
)

## Figure 9 --> figure 7 revised
prot_diversity_metadata %>%
  mutate(device = ifelse(device == "sterivex", "Cartridge membrane filter", device),
         device = ifelse(device == "membrane", "Flat membrane filter", device)) %>%
  filter(!is.na(size_fraction),!is.na(Sequencing_strategy)) %>% 
  ggplot(aes(effected_volume,Species_richness,col=device))+
  geom_jitter(width = 0.1, size = 3, alpha = 0.75)+
  facet_grid(facets = c("Sequencing_strategy","size_fraction"), scales = "free_y")+
  labs(x="Volume (L) in Log10 scale",
       y="Species richness in Log10 scale",
       col="Filter",
       title = "")+
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"),
        legend.position = "top")+
  scale_color_manual(values = qualitative_colors[c(3,5)]) + 
  scale_x_log10() + 
  scale_y_log10()

## Figure 10 --> Figure 8 revised

# A 
(figure_8_a <- 
  protist_ww_sf_10L_membrane %>% 
  ggplot(aes(x = size_fraction, y = Species_richness)) +
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = size_fraction), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Filter pore size",
       y="Species richness",
       title = "a") + #,
  facet_wrap(facets = c("Sequencing_strategy"), scale= "free") +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)]))

# B 
(figure_8_b <- 
  protist_size_fractions_100L_membrane %>% 
  ggplot(aes(x = size_fraction, y = Species_richness)) +
  geom_boxplot(outlier.shape = "cross", size = 1) +
  geom_jitter(aes(color = size_fraction), width = 0.1, size = 4, alpha = 0.75) +
  labs(x="Filter pore size",
       y="Species richness",
       title = "b") + #,
  facet_wrap(facets = c("Sequencing_strategy"), scale="free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5,7)]))

# C
(figure_8_c <- 
  protist_membrane_vs_sterivex_10L_whole_water %>% 
    mutate(device = ifelse(device == "Sterivex", "Cartridge \nmembrane \nfilter", device),
           device = ifelse(device == "Membrane", "Flat \nmembrane \nfilter", device)) %>%
  ggplot(aes(x = device,
             y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = device), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Filter type",
       y="Species richness",
       title = "c")+
  facet_wrap(facets = c("Sequencing_strategy"), scale= "free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)]))

# D
(figure_8_d <- 
  protist_2.5_vs_10_st_ww %>% 
  ggplot(aes(x = as.factor(planned_volume), y = Species_richness))+
  geom_boxplot(size = 1)+
  geom_jitter(aes(col = as.factor(planned_volume)), width = 0.1, size = 4, alpha = 0.75)+
  labs(x="Volume (L)",
       y="Species richness",
       title = "d") + #,
  facet_wrap(facets = c("Sequencing_strategy"),scale="free") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(family = "Helvetica"),
        strip.text = element_text(color = "black", size = 12, family = "Helvetica"),
        text = element_text(size = 12, family = "Helvetica"))+
  guides(col = "none")+
  scale_color_manual(values = qualitative_colors[c(3,5)]))

# arrange
cowplot::plot_grid(figure_8_a,
                   figure_8_b,
                   figure_8_c,
                   figure_8_d,
                   axis = "tlrb")

## Figure 11 --> Figure 9 revised

### Protists ###
### Prokaryotes MetaB16S
bd_prot_18S_clean_rarefied <- betadisper(
  vegdist(t(protists_18S_clean_rarefied)),
  group = pull(select(env_18S,all_of("size_fraction"))),
  type = "centroid")

### Prokaryotes MetaG

bd_prot_meta_clean_rarefied <- betadisper(
  vegdist(t(protists_metagenome_clean_rarefied)),
  group = pull(select(env_metagenome_prot,all_of("size_fraction"))),
  type = "centroid")


par(mfrow = c(2,2))
plot(mds_prot_18S,
     display = "sites", type="p", 
     main = "a \n Protists, MetaB18SV9")
points(mds_prot_18S,
       display = "sites",
       bg = env_18S$col_volume,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_18S,
     ordiellipse(mds_prot_18S,
                 size_fraction,kind="se",conf=0.95,col="grey"))
with(env_18S,
     ordihull(mds_prot_18S,
              size_fraction,col="grey"))
with(env_18S,
     ordispider(mds_prot_18S,
                size_fraction,label = T,col="grey"))

## Protists, metagenome
plot(mds_prot_metagenome,
     display = "sites", type="p", 
     main = "b \n Protists, MetaG")
points(mds_prot_metagenome,
       display = "sites",
       bg = env_metagenome_prot$col_volume,
       pch = 21, 
       col = "grey", 
       cex = 2)
with(env_metagenome_prot,
     ordiellipse(mds_prot_metagenome,
                 size_fraction,kind="se",conf=0.95,col="grey"))
with(env_metagenome_prot,
     ordihull(mds_prot_metagenome,
              size_fraction,col="grey"))
with(env_metagenome_prot,
     ordispider(mds_prot_metagenome,
                size_fraction,label = T,col="grey"))

# C
boxplot(bd_prot_18S_clean_rarefied, main = "c \n Protists, MetaB18SV9")

# D
boxplot(bd_prot_meta_clean_rarefied, main = "d \n Protists, MetaG")



## Figure 12 --> Figure 10 revised
tidy_prot_18S_clean_rarefied_taxa %>%
  mutate(FAKE_rank = ifelse(FAKE_rank == "Archaeplastida", "Others", FAKE_rank)) %>%
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>%
  group_by(FAKE_rank) %>% 
  complete(size_fraction,
           effected_volume,
           fill = list(Species_richness = 0)) %>% 
  filter(!is.na(effected_volume),
         !is.na(size_fraction)) %>%
  ggplot(aes(x = as.factor(effected_volume), 
             fill = size_fraction, y = Species_richness))+
  geom_col(position = position_dodge(preserve = "total"))+
  facet_wrap(~FAKE_rank, scales = "free")+
  scale_fill_manual(values = qualitative_colors[c(1,3,5,7)])+
  labs(fill = "Pore size (\U00B5m):",
       x = "Volume (L)", 
       y = "Species richness",
       title = "")+
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        axis.text = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 9),
        text = element_text(size = 12))

## Figure 3 --> Supplementary Figure S2
#
par(mfrow = c(2,2))
# 16S
with(env_all_prok_16S, 
     rarecurve(t(prokaryotes_16S_clean),
               col = sf_col, 
               step = 1000, 
               label = FALSE,
               xlab = "Number of reads", 
               ylab  = "Species richness", 
               lwd = 2,
               main = "a \n MetaB16SV4V5"))
legend(x = "right",
       inset = c(0.1,0.1),
       legend = c(">0.22 µm","0.22-3 µm","3-20 µm",">20 µm"),
       col = colors_sf_rarecurve, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 1,
       title = "Pore size (\U00B5m)")
# 18S
with(env_all_prot_18S, 
     rarecurve(t(protists_18S_clean),
               col = sf_col, step = 1000, 
               label = FALSE,
               xlab = "Number of reads", 
               ylab  = "Species richness", 
               lwd = 2,
               main = "b \n MetaB18SV9"))
legend(x = "right",
       inset = c(0.1,0.1),
       legend = c(">0.22 µm","0.22-3 µm","3-20 µm",">20 µm"),
       col = colors_sf_rarecurve, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 1,
       title = "Pore size (\U00B5m)")
# MetaG, Prokaryotes
with(env_all_metagenomes, 
     rarecurve(t(prokaryotes_metagenome_clean),
               col = sf_col, 
               step = 100, 
               label = FALSE,
               xlab = "Number of reads",
               ylab  = "Species richness",
               lwd = 2,
               main = "c \n MetaG, Prokaryotes"))
legend(x = "right",
       inset = c(0.1,0.1),
       legend = c(">0.22 µm","0.22-3 µm","3-20 µm",">20 µm"),
       col = colors_sf_rarecurve, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 1,
       title = "Pore size (\U00B5m)")
# MetaG, Protists
with(env_all_metagenomes, 
     rarecurve(t(protists_metagenome_clean),
               col = sf_col, step = 100, label = FALSE,
               xlab = "Number of reads",
               ylab  = "Species richness",
               lwd = 2,
               main = "d \n MetaG, Protists"))
legend(x = "right",
       inset = c(0.1,0.1),
       legend = c(">0.22 µm","0.22-3 µm","3-20 µm",">20 µm"),
       col = colors_sf_rarecurve, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 1,
       title = "Pore size (\U00B5m)")

## Supplementary Figure S1 --> Supplementary Figure S3 revised
volumes_legend <- 
  c("1L", "2.5L", "10L", 
    "30-60L", "100-120L", 
    "496L", "716-776L", "1000L")
#
par(mfrow = c(2,2))
# 16S
with(env_all_prok_16S, 
     rarecurve(t(prokaryotes_16S_clean),
               col = vol_col, 
               step = 1000, 
               label = FALSE,
               xlab = "Sequences",
               ylab  = "Species richness",
               lwd=2,
               main = "a \n MetaB16SV4V5"))
legend(x = "right",
       inset = c(0.1,0.1),
       legend = volumes_legend,
       col = colors_volume_rare, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 0.8,
       title = "Volume (L)")
# 18S
with(env_all_prot_18S, 
     rarecurve(t(protists_18S_clean),
               col = vol_col, step = 1000, label = FALSE,
               xlab = "Sequences",ylab  = "Species richness",lwd=2,
               main = "b \n MetaB18SV9"))
legend(x = "right",
       inset = c(0.1,0.1),
       legend = volumes_legend,
       col = colors_volume_rare, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 1,
       title = "Volume (L)")

# MetaG, Prokaryotes
with(env_all_metagenomes, 
     rarecurve(t(prokaryotes_metagenome_clean),
               col = vol_col, step = 100, label = FALSE,
               xlab = "Sequences",ylab  = "Species richness",lwd=2,
               main = "c \n MetaG, Prokaryotes"))
legend(x = "right",
       inset = c(0.1,0.1),
       legend = volumes_legend,
       col = colors_volume_rare, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 1,
       title = "Volume (L)")

# MetaG, Protists
with(env_all_metagenomes, 
     rarecurve(t(protists_metagenome_clean),
               col = vol_col, step = 100, label = FALSE,
               xlab = "Sequences", ylab  = "Species richness",lwd=2,
               main = "d \n MetaG, Protists")) 
legend(x = "right",
       inset = c(0.1,0.1),
       legend = volumes_legend,
       col = colors_volume_rare, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 1,
       title = "Volume (L)")

## Supplementary Figure S2 --> S4
tidy_prok_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>%
  group_by(FAKE_rank) %>% 
  complete(size_fraction,
           effected_volume,
           fill = list(Species_richness = 0)) %>% 
  filter(!is.na(effected_volume),
         !is.na(size_fraction)) %>%
  ggplot(aes(x = as.factor(effected_volume), 
             fill = size_fraction, y = Species_richness))+
  geom_col(position = position_dodge(preserve = "total"))+
  facet_wrap(~FAKE_rank, scales = "free")+
  scale_fill_manual(values = qualitative_colors[c(1,3,5,7)])+
  labs(fill = "Pore size (\U00B5m):",
       x = "Volume (L)", 
       y = "Species richness",
       title = "")+
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        axis.text = element_text(family = "Helvetica"),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 9),
        text = element_text(size = 12),
        axis.ticks.x = element_blank())




## Supplementary Figure S4 --> Sup Fig S5
tidy_prot_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>%
  mutate(FAKE_rank = case_when(
    FAKE_rank == "Other" ~ "Others",
    FAKE_rank == "Opalozoa" ~ "Others",
    FAKE_rank == "Sagenista" ~ "Others",
    FAKE_rank == "Apusozoa" ~ "Others",
    FAKE_rank == "Rhizaria" ~ "Others",
    FAKE_rank == "Alveolata" ~ "Others",
    FAKE_rank == "Amoebozoa" ~ "Others",
    FAKE_rank == FAKE_rank ~ FAKE_rank
  )) %>% 
  group_by(FAKE_rank) %>% 
  complete(size_fraction,
           effected_volume,
           fill = list(Species_richness = 0)) %>% 
  filter(!is.na(effected_volume),
         !is.na(size_fraction)) %>%
  ggplot(aes(x = as.factor(effected_volume), 
             fill = size_fraction, y = Species_richness))+
  geom_col(position = position_dodge(preserve = "total"))+
  facet_wrap(~FAKE_rank, scales = "free")+
  scale_fill_manual(values = qualitative_colors[c(1,3,5,7)])+
  labs(fill = "Pore size (\U00B5m):",
       x = "Volume (L)", 
       y = "Species richness",
       title = "")+
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        axis.text = element_text(angle = 90),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 9),
        text = element_text(size = 12))


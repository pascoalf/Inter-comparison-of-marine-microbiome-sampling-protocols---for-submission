## Data preparation ##

## Prokaryotes by 16S
#
prok_16S_clean_rarefied_taxa <- prokaryotes_16S_clean_rarefied
prok_16S_clean_rarefied_taxa$OTU <- rownames(prokaryotes_16S_clean_rarefied)
prok_16S_clean_rarefied_taxa <- 
  prok_16S_clean_rarefied_taxa %>% 
  left_join(prokaryotes_16S_OTUs,by="OTU")
#

# Tidy and add number of OTUs per taxonomy
tidy_prok_16S_clean_rarefied_taxa <-
  prok_16S_clean_rarefied_taxa %>% 
  gather(key="Sample",value="Reads",-OTU,-FAKE_rank) %>% 
  group_by(Sample,FAKE_rank) %>% 
  summarise(Species_richness = specnumber(Reads)) %>% 
  mutate(FAKE_rank = str_remove(FAKE_rank,"c__"),
         FAKE_rank = str_remove(FAKE_rank,"p__")) %>%
  ungroup() %>% 
  left_join(metadata,by=c("Sample"="run_accession")) %>%
  group_by(FAKE_rank) %>% 
  mutate(Total_species = sum(Species_richness)) %>% 
  ungroup() %>% 
  mutate(FAKE_rank = ifelse(Total_species > 100, FAKE_rank,"Others")) %>% 
  filter(Species_richness>0)

## Prokaryotes by Metagenomics
prok_metagenome_clean_rarefied_taxa <- prokaryotes_metagenome_clean_rarefied
prok_metagenome_clean_rarefied_taxa$OTU <- rownames(prokaryotes_metagenome_clean_rarefied)
prok_metagenome_clean_rarefied_taxa <- prok_metagenome_clean_rarefied_taxa %>% 
  left_join(prokaryotes_metagenome_OTUs,by="OTU")

#
# Tidy format
tidy_prok_metagenome_clean_rarefied_taxa <-
  prok_metagenome_clean_rarefied_taxa %>% 
  gather(key="Sample",value="Reads",-OTU,-FAKE_rank) %>% 
  group_by(Sample,FAKE_rank) %>% 
  summarise(Species_richness = specnumber(Reads)) %>% 
  mutate(FAKE_rank = str_remove(FAKE_rank,"c__"),
         FAKE_rank = str_remove(FAKE_rank,"p__")) %>%
  ungroup() %>% 
  left_join(metadata,by=c("Sample"="run_accession")) %>%
  group_by(FAKE_rank) %>% 
  mutate(Total_species = sum(Species_richness)) %>% 
  ungroup() %>% 
  mutate(FAKE_rank = ifelse(Total_species > 100, FAKE_rank,"Others"))%>% 
  filter(Species_richness>0)

## Protists
# MetaB18S
prot_18S_clean_rarefied_taxa <- protists_18S_clean_rarefied
prot_18S_clean_rarefied_taxa$OTU <- rownames(prot_18S_clean_rarefied_taxa)
prot_18S_clean_rarefied_taxa <- prot_18S_clean_rarefied_taxa %>% 
  left_join(protists_18S_OTUs,by="OTU")

# Tidy format
tidy_prot_18S_clean_rarefied_taxa <-
  prot_18S_clean_rarefied_taxa %>% 
  gather(key="Sample",value="Reads",-OTU,-FAKE_rank) %>% 
  group_by(Sample,FAKE_rank) %>% 
  summarise(Species_richness = specnumber(Reads)) %>% 
  mutate(FAKE_rank = str_remove(FAKE_rank,"c__"),
         FAKE_rank = str_remove(FAKE_rank,"p__")) %>% 
  ungroup() %>% 
  left_join(metadata,by=c("Sample"="run_accession"))%>% 
  filter(Species_richness>0)

# Protist, MetaG
# MetaG
prot_metagenome_clean_rarefied_taxa <- protists_metagenome_clean_rarefied
prot_metagenome_clean_rarefied_taxa$OTU <- rownames(prot_metagenome_clean_rarefied_taxa)
prot_metagenome_clean_rarefied_taxa <- prot_metagenome_clean_rarefied_taxa %>% 
  left_join(protists_metagenome_OTUs,by="OTU")

# Tidy format
tidy_prot_metagenome_clean_rarefied_taxa <-
  prot_metagenome_clean_rarefied_taxa %>% 
  gather(key="Sample",value="Reads",-OTU,-FAKE_rank) %>% 
  group_by(Sample,FAKE_rank) %>% 
  summarise(Species_richness = specnumber(Reads)) %>% 
  mutate(FAKE_rank = str_remove(FAKE_rank,"c__"),
         FAKE_rank = str_remove(FAKE_rank,"p__")) %>% 
  ungroup() %>% 
  left_join(metadata,by=c("Sample"="run_accession"))%>% 
  filter(Species_richness>0)

## Plot figures ##

## Functions to help taxonomy plots
points_taxa <- function(x, title){
  x %>% 
  ggplot(aes(effected_volume, Species_richness, col = size_fraction))+
    geom_jitter(position = position_jitter(height = 0, width = 2), alpha = 0.8, size = 2.5)+
    stat_summary(fun = mean, geom = "line", linewidth = 0.4, alpha = 0.9, position = position_dodge(width = 0.5))+
    facet_wrap(~reorder(FAKE_rank,-Species_richness),scale="free")+
    labs(x = "Volume (L)",
         y = "Number of OTU",
         colour = "Pore size (\U00B5m):",
         title = title)+
    scale_y_continuous(labels=label_number(accuracy = 1))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "top",
          axis.text = element_text(family = "Helvetica"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", size = 10),
          text = element_text(size = 12, family = "Helvetica"))+
    scale_color_manual(values = qualitative_colors[c(1,3,5,7)])
}
#
# taxonomy in collumns
col_tax_plot <- function(x, title){
  x %>% 
    ggplot(aes(x = as.factor(effected_volume), 
               fill = size_fraction, y = Species_richness))+
    geom_col(position = position_dodge(preserve = "single"))+
    facet_wrap(~FAKE_rank, scales = "free")+
    scale_fill_manual(values = qualitative_colors[c(1,3,5,7)])+
    labs(fill = "Pore size (\U00B5m):",
         x = "Volume (L)", 
         y = "Number of OTUs",
         title = title)+
    scale_y_continuous(labels=label_number(accuracy = 1))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "top",
          axis.text = element_text(family = "Helvetica"),
          strip.background = element_blank(),
          strip.text = element_text(color = "black", size = 9),
          text = element_text(size = 12),
          axis.ticks.x = element_blank())
}

# Prokaryotes, 16S
# Range 1L to 100L
# Plot number of OTUs per Fake rank
tidy_prok_16S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         effected_volume <= 100) %>% 
  mutate(FAKE_rank = ifelse(FAKE_rank == "Candidatus Marinimicrobia", 
                            "Candidatus \nMarinimicrobia", FAKE_rank),
         FAKE_rank = ifelse(FAKE_rank == "Candidate Phylum", 
                            "Candidate \nPhylum", FAKE_rank)) %>% 
  points_taxa(title = "Prokaryotes 16S | 2.5L to 100L") 

## Alternative
tidy_prok_16S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),
         !is.na(size_fraction),
         effected_volume <= 100) %>% 
  mutate(FAKE_rank = ifelse(FAKE_rank == "Candidatus Marinimicrobia", 
                            "Candidatus \nMarinimicrobia", FAKE_rank),
         FAKE_rank = ifelse(FAKE_rank == "Candidate Phylum", 
                            "Candidate \nPhylum", FAKE_rank)) %>% 
  col_tax_plot(title = "Prokaryotes 16S | 2.5L to 100L")

## Full range of volumes
# Plot number of OTUs per Fake rank
tidy_prok_16S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>% 
  mutate(FAKE_rank = ifelse(FAKE_rank == "Candidatus Marinimicrobia", "Candidatus \nMarinimicrobia", FAKE_rank), 
         FAKE_rank = ifelse(FAKE_rank == "Candidate Phylum", "Candidate \nPhylum", FAKE_rank)) %>% 
  points_taxa(title = "Prokaryotes 16S | 2.5L to 1000L")

# Alternative
tidy_prok_16S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>% 
  mutate(FAKE_rank = ifelse(FAKE_rank == "Candidatus Marinimicrobia", "Candidatus \nMarinimicrobia", FAKE_rank), 
         FAKE_rank = ifelse(FAKE_rank == "Candidate Phylum", "Candidate \nPhylum", FAKE_rank)) %>% 
  col_tax_plot(title = "Prokaryotes 16S | 2.5L to 1000L") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.95, size = 8))
  
## Prokaryotes, metagenome

# Range 1L to 100L
# Plot number of OTUs per Fake rank
tidy_prok_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         effected_volume <= 100) %>%
  points_taxa(title = "Prokaryotes MetaG | 1L to 100L") +
  geom_point(size = 3)

## Alternative
tidy_prok_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         effected_volume <= 100) %>%
  col_tax_plot(title = "Prokaryotes MetaG | 1L to 100L")

## Full range of volumes
# Plot number of OTUs per Fake rank
tidy_prok_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>%
  points_taxa(title = "Prokaryotes MetaG | 1L to 1000L") +
  geom_point(size = 3, alpha = 0.6)

## Alternative
tidy_prok_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>%
  col_tax_plot(title = "Prokaryotes MetaG | 1L to 1000L") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9))

# Protists, MetaB18S

## Range 1L to 100L
# Plot number of OTUs per Fake rank
tidy_prot_18S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         effected_volume <= 100) %>% 
  mutate(FAKE_rank = ifelse(FAKE_rank == "Archaeplastida", "Others", FAKE_rank)) %>% 
  points_taxa(title = "Protists 18S | 2.5L to 100L") + 
  geom_point(size = 3, alpha = 0.6)

## Alternative
tidy_prot_18S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         effected_volume <= 100) %>% 
  mutate(FAKE_rank = ifelse(FAKE_rank == "Archaeplastida", "Others", FAKE_rank)) %>% 
  col_tax_plot(title = "Protists 18S | 2.5L to 100L")

## Full range of volumes
# Plot number of OTUs per Fake rank
tidy_prot_18S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>%
  mutate(FAKE_rank = ifelse(FAKE_rank == "Archaeplastida", "Others", FAKE_rank)) %>% 
  points_taxa(title = "Protists 18S | 2.5L to 1000L") +
  geom_point(size = 2, alpha = 0.6)

## Alternative
tidy_prot_18S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction)) %>%
  mutate(FAKE_rank = ifelse(FAKE_rank == "Archaeplastida", "Others", FAKE_rank)) %>% 
  col_tax_plot(title = "Protists 18S | 2.5L to 1000L") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9))

# Protist, MetaG

## Range 1L to 100L
# Plot number of OTUs per Fake rank
tidy_prot_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         effected_volume <= 100) %>%
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
  points_taxa(title = "Protists MetaG | 1L to 100L") + 
  geom_point(size = 3, alpha = 0.6)

## Alternative
tidy_prot_metagenome_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         effected_volume <= 100) %>%
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
  col_tax_plot(title = "Protists MetaG | 1L to 100L")

## Full range of volumes
# Plot number of OTUs per Fake rank
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
  points_taxa(title = "Protists MetaG | 1L to 1000L") +
  geom_point(size = 2, alpha = 0.6)

## Alternative
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
  col_tax_plot(title = "Protists MetaG | 1L to 1000L") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9))

# Candidate marinimicrobia

grid.arrange(
tidy_prok_16S_clean_rarefied_taxa %>% 
  filter(!is.na(effected_volume),!is.na(size_fraction),
         FAKE_rank == "Candidatus Marinimicrobia",
         effected_volume %in% c(100)) %>% 
  ggplot(aes(size_fraction, Species_richness))+
  geom_col(aes(fill = size_fraction, col = size_fraction))+
  labs(y = "Number of OTU", 
       x = "",
       title = "16S | membrane 100L | Candidatus Marinimicrobia") + 
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
       y = "Relative abundance (%)")+
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

# end

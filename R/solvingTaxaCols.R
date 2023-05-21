tidy_prok_16S_clean_rarefied_taxa %>%
  group_by(FAKE_rank) %>% 
  complete(size_fraction,
           effected_volume,
           fill = list(Species_richness = 0)) %>% 
  filter(!is.na(effected_volume),
         !is.na(size_fraction),
         effected_volume <= 100) %>%
  ggplot(aes(x = as.factor(effected_volume), 
             fill = size_fraction, y = Species_richness))+
  geom_col(position = position_dodge(preserve = "total"))+
  facet_wrap(~FAKE_rank, scales = "free")+
  scale_fill_manual(values = qualitative_colors[c(1,3,5,7)])+
  labs(fill = "Pore size (\U00B5m):",
       x = "Volume (L)", 
       y = "Number of OTUs",
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



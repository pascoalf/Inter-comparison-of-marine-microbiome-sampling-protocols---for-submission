## mock community figures

# removed singletons, kept rare biosphere

phylumExpectedObserved_plot <- 
  phylumExpectedObserved %>% 
  filter(RelativeAbundance > 0) %>% 
  ggplot(aes(x = reorder(Phylum, RelativeAbundance), y = RelativeAbundance, fill = SequencingMachine)) + 
  geom_col(position = position_dodge()) + 
  coord_flip() + 
  theme_bw() +
  scale_fill_manual(values = qualitative_colors[c(2,3,4)]) + 
  labs(x = "Phylum",
       y = "Relative abundance (%)",
       fill = "Source",
       title = "Phylum level",
       subtitle = "Note: Removed singletons") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank(), text = element_text(size = 10)) + 
  scale_fill_brewer(palette = 2, type = "qual")

classExpectedObserved_plot <-
  classExpectedObserved %>% 
  filter(RelativeAbundance > 0) %>% 
  ggplot(aes(reorder(Class, RelativeAbundance), RelativeAbundance, fill = SequencingMachine))+
  geom_col(position = position_dodge()) + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Class",
       y = "Relative abundance (%)",
       fill = "Source",
       title = "Class level",
       subtitle = "Note: removed singletons") + 
  theme(legend.position = "top", 
        panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        axis.text.y = element_text(size = 7)) + 
  scale_fill_brewer(palette = 2, type = "qual")


percentage_classification <- stag_tidy %>% 
  pivot_longer(cols = contains("check"), 
               values_to = "Classification", 
               names_to = "TaxonomicLevel") %>%
  group_by(SequencingMachine, TaxonomicLevel) %>% 
  filter(Abundance > 1) %>% 
  count(Classification) %>% 
  mutate(Percentage = n*100/sum(n)) %>% 
  mutate(TaxonomicLevel = str_remove(TaxonomicLevel, "_check")) %>%
  mutate(TaxonomicLevel = factor(TaxonomicLevel, levels = c("Phylum", "Class", "Order"))) %>% 
  ggplot(aes(TaxonomicLevel, y = Percentage, fill = Classification))+
  geom_col() + 
  facet_grid(~SequencingMachine) +
  scale_fill_manual(values = c("dodgerblue4", "indianred3")) + 
  geom_hline(yintercept = 50, linetype = "dashed") + 
  theme_bw()+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 15), 
        text = element_text(size = 10)) +
  labs(x = "Taxonomic level",
       y = "Percentage",
       fill = "Classification was",
       subtitle = "Note: removed singletons",
       title = "Percentage of correct vs wrong classifications")


gridExtra::grid.arrange(phylumExpectedObserved_plot,
                        classExpectedObserved_plot, 
                        percentage_classification, ncol = 3)


### Remove rare biosphere
phylumExpectedObserved_plot_abundant <- 
  phylumExpectedObserved %>% 
  filter(RelativeAbundance > 0.1) %>% 
  ggplot(aes(x = reorder(Phylum, RelativeAbundance), y = RelativeAbundance, fill = SequencingMachine)) + 
  geom_col(position = position_dodge()) + 
  coord_flip() + 
  theme_bw() +
  scale_fill_manual(values = qualitative_colors[c(2,3,4)]) + 
  labs(x = "Phylum",
       y = "Relative abundance (%)",
       fill = "Source",
       title = "Phylum level",
       subtitle = "Note: relative abundance \U2265 0.1%") + 
  theme(legend.position = "top",
        panel.grid.minor = element_blank(), text = element_text(size = 10)) + 
  scale_fill_brewer(palette = 2, type = "qual")

classExpectedObserved_plot_abundant <-
  classExpectedObserved %>% 
  filter(RelativeAbundance > 0.1) %>% 
  ggplot(aes(reorder(Class, RelativeAbundance), RelativeAbundance, fill = SequencingMachine))+
  geom_col(position = position_dodge()) + 
  coord_flip() + 
  theme_bw() + 
  labs(x = "Class",
       y = "Relative abundance (%)",
       fill = "Source",
       title = "Class level",
       subtitle = "Note: relative abundance \U2265 0.1%") + 
  theme(legend.position = "top", 
        panel.grid.minor = element_blank(),
        text = element_text(size = 10)) + 
  scale_fill_brewer(palette = 2, type = "qual")


percentage_classification_abundant <- stag_tidy %>% 
  filter(RelativeAbundance > 0.1) %>% 
  pivot_longer(cols = contains("check"), 
               values_to = "Classification", 
               names_to = "TaxonomicLevel") %>%
  group_by(SequencingMachine, TaxonomicLevel) %>% 
  filter(Abundance > 1) %>% 
  count(Classification) %>% 
  mutate(Percentage = n*100/sum(n)) %>% 
  mutate(TaxonomicLevel = str_remove(TaxonomicLevel, "_check")) %>%
  mutate(TaxonomicLevel = factor(TaxonomicLevel, levels = c("Phylum", "Class", "Order"))) %>% 
  ggplot(aes(TaxonomicLevel, y = Percentage, fill = Classification))+
  geom_col() + 
  facet_grid(~SequencingMachine) +
  scale_fill_manual(values = c("dodgerblue4", "indianred3")) + 
  geom_hline(yintercept = 50, linetype = "dashed") + 
  theme_bw()+
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 15), text = element_text(size = 10)) +
  labs(x = "Taxonomic level",
       y = "Percentage",
       fill = "Classification was",
       subtitle = "Note: relative abundance \U2265 0.1%",
       title = "Percentage of correct vs wrong classifications")






gridExtra::grid.arrange(phylumExpectedObserved_plot,
                        classExpectedObserved_plot, 
                        percentage_classification,
                        phylumExpectedObserved_plot_abundant,
                        classExpectedObserved_plot_abundant, 
                        percentage_classification_abundant, ncol = 3)





## betadisper summary


bd_prok_16S_clean_rarefied_tukey <- 
  TukeyHSD(bd_prok_16S_clean_rarefied, type = "centroid")$group %>% 
  as.data.frame() %>% 
  mutate(Comparison = rownames(.),
         Sequencing_strategy = "MetaB16SV4V5",
         Taxonomic_group = "Prokaryotic")


bd_prok_meta_clean_rarefied_tukey <-   
  TukeyHSD(bd_prok_meta_clean_rarefied, type = "centroid")$group %>% 
  as.data.frame() %>% 
  mutate(Comparison = rownames(.),
         Sequencing_strategy = "MetaG",
         Taxonomic_group = "Prokaryotic")

betadisper_prokaryotes_tukey <- rbind(bd_prok_16S_clean_rarefied_tukey, bd_prok_meta_clean_rarefied_tukey)


bd_prot_18S_clean_rarefied_tukey <-    
  TukeyHSD(bd_prot_18S_clean_rarefied, type = "centroid")$group %>% 
  as.data.frame() %>% 
  mutate(Comparison = rownames(.),
         Sequencing_strategy = "MetaB18SV9",
         Taxonomic_group = "Protist")

bd_prot_meta_clean_rarefied_tukey <-    
  TukeyHSD(bd_prot_meta_clean_rarefied, type = "centroid")$group %>% 
  as.data.frame() %>% 
  mutate(Comparison = rownames(.),
         Sequencing_strategy = "MetaG",
         Taxonomic_group = "Protist")


betadisper_protist_tukey <-  rbind(bd_prot_18S_clean_rarefied_tukey, bd_prot_meta_clean_rarefied_tukey)



write.csv(betadisper_prokaryotes_tukey, "./Revised_results/betadisper_prokaryotes.csv")
write.csv(betadisper_protist_tukey, "./Revised_results/betadisper_protists.csv")

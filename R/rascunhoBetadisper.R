bd_prok_16S_clean_rarefied <- betadisper(
  vegdist(t(prokaryotes_16S_clean_rarefied)),
  group = pull(select(env_16S,all_of("size_fraction"))),
  type = "centroid")

## S3 method for class 'betadisper'
plot(bd_prok_16S_clean_rarefied, 
     main = "Prokaryotes| 16S")

## S3 method for class 'betadisper'
boxplot(bd_prok_16S_clean_rarefied)

## S3 method for class 'betadisper'
TukeyHSD(betadisper(
  vegdist(t(prokaryotes_16S_clean_rarefied)),
  group = pull(select(env_16S,all_of("size_fraction"))),
  type = "centroid"))


bd_prok_16S_clean_rarefied_df <- 
  bd_prok_16S_clean_rarefied$distances %>% 
  as.data.frame() %>% 
  rename(Distance = ".") %>% 
  mutate(Sample = rownames(.)) %>% 
  cbind(data.frame(group = bd_prok_16S_clean_rarefied$group))


bd_prok_16S_clean_rarefied_df %>% 
  ggplot(aes(group, Distance)) +
  geom_boxplot(outlier.shape = "cross") + 
  geom_jitter(alpha = 0.5, position = position_dodge2(width = 0.25)) + 
  labs(title = "C\nProkaryotes 16S",
       y = "Distance to centroid")

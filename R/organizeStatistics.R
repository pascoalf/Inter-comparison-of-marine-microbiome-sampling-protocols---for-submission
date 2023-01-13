## gather statistical tests 

# Prokaryotes
prokaryotes_statistical_tests <- 
# ww vs sf
prokaryotes_wilcox_test_ww_vs_sf_10L_membrane %>% 
  full_join(prokaryotes_kruskal_test_ww_vs_sf_10L_membrane) %>% 
  full_join(prokaryotes_dunn_ww_vs_sf_10L_membrane) %>% 
# sf 100L
  full_join(prokaryotes_size_fractions_100L_membrane_kruskal_test_16S) %>% 
  full_join(prokaryotes_size_fractions_100L_membrane_dunn_test_16S) %>% 
  full_join(prokaryotes_size_fractions_100L_membrane_wilcox_metagenomes) %>% 
# device
  full_join(prok_wilc_dev_10L_ww) %>% 
# 2.5L vs 10L
  full_join(prok_2.5L_vs_10L_wilcoxon)

## To export
#write.csv(sapply(prokaryotes_statistical_tests,as.character),"./output/prokaryotes_statistical_tests.csv")

# Protist
protists_statistical_tests <- 
# ww vs sf
protist_ww_sf_10L_membrane_kruskall_test %>% 
# sf 100L
full_join(protist_size_fractions_100L_membrane_kruskall) %>% 
full_join(protist_size_fractions_100L_membrane_dunn) %>% 
# device
full_join(prot_wilc_m_vs_st_10L_ww) %>% 
# 2.5L vs 10L
full_join(prot_2.5_vs_10_st_ww_wilcox)

## To export
#write.csv(sapply(protists_statistical_tests,as.character),"./output/protists_statistical_tests.csv")


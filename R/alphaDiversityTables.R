## Prokaryotes, 16S
# 16S samples
samples_16S <- 
  metadata %>%
  filter(Sequencing_strategy == "MetaB16SV4V5",
         !is.na(run_accession),
         !str_detect(sample_alias,"MOCK")) %>% 
  select(run_accession,Sequencing_strategy) 

# Samples for prokaryotes
prokaryotes_16S_clean <- 
  curated_taxonomy_prokaryotes_16S %>% 
  select(samples_16S$run_accession,FAKE_rank)

# Taxonomy of each OTU
prokaryotes_16S_OTUs <- tibble(FAKE_rank = prokaryotes_16S_clean$FAKE_rank,
                               OTU = paste0("OTU",1:nrow(prokaryotes_16S_clean)))

rownames(prokaryotes_16S_clean) <- prokaryotes_16S_OTUs$OTU
prokaryotes_16S_clean$FAKE_rank <- NULL

all_16S_samples <- prokaryotes_16S_clean %>% names()

prokaryotes_16S_diversity_not_rarefied <- 
  tibble(Sample = all_16S_samples,
         Species_richness = sapply(prokaryotes_16S_clean,specnumber),
         Shannon_index = sapply(prokaryotes_16S_clean,diversity,index="shannon"),
         Reads = sapply(prokaryotes_16S_clean,sum),
         Taxonomic_group = "Prokaryotic")

prokaryotes_16S_diversity_not_rarefied %>%
  left_join(metadata, by = c("Sample" = "run_accession")) %>% 
  ggplot(aes(Sample,Reads,col=device))+
  geom_point()+
  coord_flip()+
  geom_hline(yintercept = 250000)

# Rarefaction of prokaryotes, 16S
prok_16S_reads <- 
  tibble(Sample = samples_16S$run_accession,
         Reads = sapply(prokaryotes_16S_clean,sum))

# View total number of reads
prok_16S_reads %>% 
  left_join(metadata,by=c("Sample" = "run_accession")) %>% 
  ggplot(aes(Sample,Reads,fill=device))+
  geom_col()+
  coord_flip()+
  geom_hline(yintercept = 250000)+
  theme_bw()+
  theme(panel.grid.minor.x = element_line(color="red"))

# Rarefaction for prokaryotes 16S
rarefaction_prok_16S <- 250000

# Check proportion of lost samples for 250 000 reads
mean(prok_16S_reads <= rarefaction_prok_16S)*100
sum(prok_16S_reads <= rarefaction_prok_16S)

# Rarefy 16S to 250 000 reads
# select samples good to rarefy
samples_of_prokaryotes_16S_to_rarefy <- 
  prok_16S_reads %>%
  filter(Reads >= rarefaction_prok_16S) %>% 
  select(Sample)

prokaryotes_16S_clean_to_rarefy <- 
  prokaryotes_16S_clean %>% select(samples_of_prokaryotes_16S_to_rarefy$Sample)

# Remove singletons
prokaryotes_16S_clean_to_rarefy[prokaryotes_16S_clean_to_rarefy==1] <- 0

# prepare table for rarefaction
t_prok_16_to_rarefy <- t(prokaryotes_16S_clean_to_rarefy)
rownames(t_prok_16_to_rarefy) <- colnames(prokaryotes_16S_clean_to_rarefy)
colnames(t_prok_16_to_rarefy) <- rownames(prokaryotes_16S_clean_to_rarefy)

# Rarefaction
prokaryotes_16S_clean_rarefied <- 
  rrarefy(t_prok_16_to_rarefy, sample = rarefaction_prok_16S) %>% 
  t() %>% as.data.frame()

# Rarefied diversity
## note: added more alternative metrics
prokaryotic_16S_diversity <- 
  tibble(Sample = samples_of_prokaryotes_16S_to_rarefy$Sample,
         Species_richness = sapply(prokaryotes_16S_clean_rarefied,specnumber),
         Shannon_index = sapply(prokaryotes_16S_clean_rarefied,diversity,index="shannon"),
         Simpson_index = sapply(prokaryotes_16S_clean_rarefied,diversity,index="simpson"),
         Reads = sapply(prokaryotes_16S_clean_rarefied,sum))

# Sanity check on rarefaction
prokaryotic_16S_diversity %>% 
  ggplot(aes(Sample,Reads)) + 
  geom_col() + 
  coord_flip()

## Prokaryotes, Metagenomes
# select relevant samples for prokaryotes
metagenome_samples <- metadata %>% 
  filter(Sequencing_strategy == "MetaG") %>% 
  drop_na(Sequencing_strategy) %>% 
  filter(!str_detect(sample_alias,"MOCK"))

prokaryotes_metagenome_clean <- 
  curated_taxonomy_prokaryotes_metagenome %>% 
  select(metagenome_samples$run_accession,FAKE_rank)

prokaryotes_metagenomes_diversity_not_rarefied <- 
  tibble(Sample = metagenome_samples$run_accession,
         Species_richness = sapply(select(prokaryotes_metagenome_clean, !FAKE_rank), specnumber),
         Shannon_index = sapply(select(prokaryotes_metagenome_clean, !FAKE_rank), diversity, index = "shannon"),
         Reads = sapply(select(prokaryotes_metagenome_clean, !FAKE_rank), sum),
         Taxonomic_group = "Prokaryotic")

# Taxonomy of each OTU
prokaryotes_metagenome_OTUs <- tibble(FAKE_rank = prokaryotes_metagenome_clean$FAKE_rank,
                               OTU = paste0("OTU",1:nrow(prokaryotes_metagenome_clean)))

rownames(prokaryotes_metagenome_clean) <- prokaryotes_metagenome_OTUs$OTU
prokaryotes_metagenome_clean$FAKE_rank <- NULL

# Rarefaction of prokaryotes, metagenomes
prok_metagenome_reads <- 
  tibble(Sample = metagenome_samples$run_accession,
         Reads = sapply(prokaryotes_metagenome_clean, sum))

# View total number of reads
prok_metagenome_reads %>% 
  ggplot(aes(Sample,Reads))+
  geom_col()+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor.x = element_line(color="red"))

# Rarefaction for prokaryotes metagenome
rarefaction_prok_metagenomes <- 10000

# Check proportion of lost samples for 250 000 reads
mean(prok_metagenome_reads$Reads <= rarefaction_prok_metagenomes)*100
sum(prok_metagenome_reads$Reads <= rarefaction_prok_metagenomes)

# Rarefy prokaryotic metagenomes to 10 000 reads
# select samples good to rarefy
samples_of_prokaryotes_metagenomes_to_rarefy <- 
  prok_metagenome_reads %>%
  filter(Reads >= rarefaction_prok_metagenomes) %>% 
  select(Sample)

prokaryotes_metagenome_clean_to_rarefy <- 
  prokaryotes_metagenome_clean %>% 
  select(samples_of_prokaryotes_metagenomes_to_rarefy$Sample)

# Remove singletons
#prokaryotes_metagenome_clean_to_rarefy[prokaryotes_metagenome_clean_to_rarefy==1,] <- 0

t_prok_metagenome_to_rarefy <- t(prokaryotes_metagenome_clean_to_rarefy)
rownames(t_prok_metagenome_to_rarefy) <- colnames(prokaryotes_metagenome_clean_to_rarefy)
colnames(t_prok_metagenome_to_rarefy) <- rownames(prokaryotes_metagenome_clean_to_rarefy)

# Rarefaction
prokaryotes_metagenome_clean_rarefied <- 
  rrarefy(t_prok_metagenome_to_rarefy, sample = rarefaction_prok_metagenomes) %>% 
  t() %>% as.data.frame()

# Rarefied diversity
prokaryotic_metagenome_diversity <- 
  tibble(Sample = samples_of_prokaryotes_metagenomes_to_rarefy$Sample,
         Species_richness = sapply(prokaryotes_metagenome_clean_rarefied,specnumber),
         Shannon_index = sapply(prokaryotes_metagenome_clean_rarefied,diversity,index="shannon"),
         Simpson_index = sapply(prokaryotes_metagenome_clean_rarefied,diversity,index="simpson"),
         Reads = sapply(prokaryotes_metagenome_clean_rarefied,sum))

# Sanity check on rarefaction
prokaryotic_metagenome_diversity %>% 
  ggplot(aes(Sample,Reads)) + 
  geom_col() + 
  coord_flip()

## Add metadata to diversity tables
# 16S
prok_16S_diversity_rarefied_metadata <-
  prokaryotic_16S_diversity %>%
  left_join(metadata, by = c("Sample" = "run_accession"))

# prok, metagenome
prok_metagenome_diversity_rarefied_metadata <-
  prokaryotic_metagenome_diversity %>%
  left_join(metadata, by = c("Sample" = "run_accession"))

## Join 16S and metagenome for prokaryotes
prok_diversity_metadata <- 
  prok_16S_diversity_rarefied_metadata %>% 
  rbind(prok_metagenome_diversity_rarefied_metadata)

## To export
#write.csv(prok_diversity_metadata, "./output/otus_summary_prokaryotes.csv")

# Protists 

## protists, 18S
# 18S samples
samples_18S <- 
  metadata %>% 
  filter(Sequencing_strategy == "MetaB18SV9",
         !is.na(run_accession),
         !str_detect(sample_alias,"MOCK")) %>% 
  select(run_accession,Sequencing_strategy)

# Samples for protists
protists_18S_clean <- 
  curated_taxonomy_protists_18S %>% 
  select(samples_18S$run_accession,FAKE_rank)

# Taxonomy of each OTU
protists_18S_OTUs <- tibble(FAKE_rank = protists_18S_clean$FAKE_rank,
                            OTU = paste0("OTU",1:nrow(protists_18S_clean)))

protists_18S_clean$FAKE_rank <- NULL
rownames(protists_18S_clean) <- protists_18S_OTUs$OTU

protists_18S_diversity_not_rarefied <- 
  tibble(Sample = samples_18S$run_accession,
         Species_richness = sapply(protists_18S_clean,specnumber),
         Shannon_index = sapply(protists_18S_clean,diversity,index="shannon"),
         Reads = sapply(protists_18S_clean,sum),
         Taxonomic_group = "Protist")

protists_18S_diversity_not_rarefied %>% 
  ggplot(aes(Sample,Reads))+
  geom_point()+
  coord_flip()+
  geom_hline(yintercept = 250000)

# Rarefaction of protists, 18S
prot_18S_reads <- 
  tibble(Sample = samples_18S$run_accession,
         Reads = sapply(protists_18S_clean,sum))

# View total number of reads
prot_18S_reads %>% 
  ggplot(aes(Sample,Reads))+
  geom_col()+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor.x = element_line(color="red"))

# Rarefaction for protists 18S
rarefaction_prot_18S <- 100000

# Check proportion of lost samples for 250 000 reads
mean(prot_18S_reads <= rarefaction_prot_18S)*100
sum(prot_18S_reads <= rarefaction_prot_18S)

# Rarefy 18S to 100 000 reads
# select samples good to rarefy
samples_of_protists_18S_to_rarefy <- 
  prot_18S_reads %>%
  filter(Reads >= rarefaction_prot_18S) %>% 
  select(Sample)

protists_18S_clean_to_rarefy <- 
  as.matrix(protists_18S_clean)[,samples_of_protists_18S_to_rarefy$Sample]

# Remove singletons
#protists_18S_clean_to_rarefy[protists_18S_clean_to_rarefy==1] <- 0

# prepare table for rarefaction
t_prot_18S_to_rarefy <- t(protists_18S_clean_to_rarefy)
rownames(t_prot_18S_to_rarefy) <- colnames(protists_18S_clean_to_rarefy)
colnames(t_prot_18S_to_rarefy) <- rownames(protists_18S_clean_to_rarefy)

# Rarefaction
protists_18S_clean_rarefied <- 
  rrarefy(t_prot_18S_to_rarefy, sample = rarefaction_prot_18S) %>% 
  t() %>% as.data.frame()

# Rarefied diversity
protist_18S_diversity <- 
  tibble(Sample = samples_of_protists_18S_to_rarefy$Sample,
         Species_richness = sapply(protists_18S_clean_rarefied,specnumber),
         Shannon_index = sapply(protists_18S_clean_rarefied,diversity,index="shannon"),
         Simpson_index = sapply(protists_18S_clean_rarefied,diversity,index="simpson"),
         Reads = sapply(protists_18S_clean_rarefied,sum))

# Sanity check on rarefaction
protist_18S_diversity %>% 
  ggplot(aes(Sample,Reads)) + 
  geom_col() + 
  coord_flip()


## protists, Metagenomes
# select relevant samples for protists

#
protists_metagenome_clean <- 
  curated_taxonomy_protist_metagenome %>% 
  select(metagenome_samples$run_accession)

protist_metagenome_diversity_not_rarefied <- 
  tibble(Sample = metagenome_samples$run_accession,
         Species_richness = sapply(protists_metagenome_clean,specnumber),
         Shannon_index = sapply(protists_metagenome_clean,diversity,index="shannon"),
         Reads = sapply(protists_metagenome_clean,sum),
         Taxonomic_group = "Protist")


# Taxonomy of each OTU
protists_metagenome_OTUs <- 
  tibble(FAKE_rank = curated_taxonomy_protist_metagenome$FAKE_rank,
         OTU = paste0("OTU",1:nrow(protists_metagenome_clean)))

rownames(protists_metagenome_clean) <- protists_metagenome_OTUs$OTU

# Rarefaction of protists, metagenomes
prot_metagenome_reads <- 
  tibble(Sample = metagenome_samples$run_accession,
         Reads = sapply(protists_metagenome_clean,sum))

# View total number of reads
prot_metagenome_reads %>% 
  left_join(metadata,by=c("Sample"="run_accession")) %>% 
  ggplot(aes(Sample,Reads,fill=size_fraction))+
  geom_col(position = "dodge")+
  geom_hline(yintercept = c(1000,2000))+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor.x = element_line(color="red"))

# Rarefaction for protists metagenome
rarefaction_prot_metagenomes <- 1000

# Check proportion of lost samples for 2000 reads
mean(prot_metagenome_reads$Reads <= rarefaction_prot_metagenomes)*100
sum(prot_metagenome_reads$Reads <= rarefaction_prot_metagenomes)

# Rarefy protaryotic metagenomes to 2000 reads
# select samples good to rarefy
samples_of_protists_metagenomes_to_rarefy <- 
  prot_metagenome_reads %>%
  filter(Reads >= rarefaction_prot_metagenomes) %>% 
  select(Sample)

protists_metagenome_clean_to_rarefy <- 
  protists_metagenome_clean[samples_of_protists_metagenomes_to_rarefy$Sample]

# Remove singletons
#protists_metagenome_clean_to_rarefy[protists_metagenome_clean_to_rarefy==1,] <- 0

t_prot_metagenome_to_rarefy <- t(protists_metagenome_clean_to_rarefy)
rownames(t_prot_metagenome_to_rarefy) <- colnames(protists_metagenome_clean_to_rarefy)
colnames(t_prot_metagenome_to_rarefy) <- rownames(protists_metagenome_clean)

# Rarefaction
protists_metagenome_clean_rarefied <- 
  rrarefy(t_prot_metagenome_to_rarefy, sample = rarefaction_prot_metagenomes) %>% 
  t() %>% as.data.frame()

# Rarefied diversity
protist_metagenome_diversity <- 
  tibble(Sample = samples_of_protists_metagenomes_to_rarefy$Sample,
         Species_richness = sapply(protists_metagenome_clean_rarefied,specnumber),
         Shannon_index = sapply(protists_metagenome_clean_rarefied,diversity,index="shannon"),
         Simpson_index = sapply(protists_metagenome_clean_rarefied,diversity,index="simpson"),
         Reads = sapply(protists_metagenome_clean_rarefied,sum))

# Sanity check on rarefaction
protist_metagenome_diversity %>% 
  ggplot(aes(Sample,Reads)) + 
  geom_col() + 
  coord_flip()

## Add metadata to diversity tables
# 18S
prot_18S_diversity_rarefied_metadata <-
  protist_18S_diversity %>%
  left_join(metadata, by = c("Sample" = "run_accession"))

# prot, metagenome
prot_metagenome_diversity_rarefied_metadata <-
  protist_metagenome_diversity %>%
  left_join(metadata, by = c("Sample" = "run_accession"))

## Join 18S and metagenome for protists
prot_diversity_metadata <- 
  prot_18S_diversity_rarefied_metadata %>% 
  rbind(prot_metagenome_diversity_rarefied_metadata)

## To export 
#write.csv(prot_diversity_metadata, "./output/otus_summary_protists.csv")


## Summarise discarded samples, after rarefaction
discarded_16S <- 
  prok_16S_reads %>% 
  left_join(metadata,by=c("Sample" = "run_accession")) %>% 
  filter(Reads < rarefaction_prok_16S) %>% 
  select(Sample,Sequencing_strategy) %>% 
  mutate(Taxonomic_group = "Prokaryotic")

discarded_18S <- 
prot_18S_reads %>% 
  left_join(metadata,by=c("Sample" = "run_accession")) %>% 
  filter(Reads < rarefaction_prot_18S) %>% 
  select(Sample,Sequencing_strategy)%>% 
  mutate(Taxonomic_group = "Protist")

discarded_meta_prok <- 
prot_metagenome_reads %>% 
  left_join(metadata,by=c("Sample" = "run_accession")) %>% 
  filter(Reads < rarefaction_prok_metagenomes) %>% 
  select(Sample,Sequencing_strategy)%>% 
  mutate(Taxonomic_group = "Prokaryotic")

discarded_meta_prot <- 
prot_metagenome_reads %>% 
  left_join(metadata,by=c("Sample" = "run_accession")) %>% 
  filter(Reads < rarefaction_prot_metagenomes) %>% 
  select(Sample,Sequencing_strategy)%>% 
  mutate(Taxonomic_group = "Protist")

discarded_samples <- 
  rbind(discarded_16S,
      discarded_18S,
      discarded_meta_prok,
      discarded_meta_prot)

## To export 
#write.csv(discarded_samples,"./output/discarded_samples.csv")


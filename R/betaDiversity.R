# Beta diversity | Volumes and size fractions
  library(RColorBrewer)

## Data preparation ##
# Prokaryotes, 16S
env_16S <- metadata %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5",
         run_accession %in% samples_of_prokaryotes_16S_to_rarefy$Sample) %>% 
  select(size_fraction,effected_volume,sample = run_accession)

row.names(env_16S) <- env_16S$sample

# MDS calculation for 16S
mds_prok_16S <-
  metaMDS(t(prokaryotes_16S_clean_rarefied))

env_16S <- env_16S %>% 
  mutate(col_volume = 
           case_when(
             effected_volume == 2.5 ~ brewer.pal(n = 9,name = "Blues")[1],
             effected_volume == 10 ~ brewer.pal(n = 9,name = "Blues")[2],
             effected_volume == 30 ~ brewer.pal(n = 9,name = "Blues")[3],
             effected_volume == 60 ~ brewer.pal(n = 9,name = "Blues")[4],
             effected_volume == 100 ~ brewer.pal(n = 9,name = "Blues")[5],
             effected_volume == 120 ~ brewer.pal(n = 9,name = "Blues")[6],
             effected_volume == 496 ~ brewer.pal(n = 9,name = "Blues")[7],
             effected_volume == 716 ~ brewer.pal(n = 9,name = "Blues")[8],
             effected_volume == 776 ~ brewer.pal(n = 9,name = "Blues")[9],
             effected_volume == 1000 ~ "black")) 

# Prokaryotes, metagenome
env_metagenome_prok <- metadata %>% 
  filter(Sequencing_strategy == "MetaG", 
         run_accession %in% samples_of_prokaryotes_metagenomes_to_rarefy$Sample)  %>% 
  select(size_fraction,effected_volume,sample = run_accession)

row.names(env_metagenome_prok) <- env_metagenome_prok$sample

# MDS calculation for metagenome
mds_prok_metagenome <-
  metaMDS(t(prokaryotes_metagenome_clean_rarefied))

env_metagenome_prok <- env_metagenome_prok %>% 
  mutate(col_volume = 
           case_when(
             effected_volume <= 2.5 ~ brewer.pal(n = 9,name = "Blues")[1],
             effected_volume == 10 ~ brewer.pal(n = 9,name = "Blues")[2],
             effected_volume == 30 ~ brewer.pal(n = 9,name = "Blues")[3],
             effected_volume == 60 ~ brewer.pal(n = 9,name = "Blues")[4],
             effected_volume == 100 ~ brewer.pal(n = 9,name = "Blues")[5],
             effected_volume == 120 ~ brewer.pal(n = 9,name = "Blues")[6],
             effected_volume == 496 ~ brewer.pal(n = 9,name = "Blues")[7],
             effected_volume == 716 ~ brewer.pal(n = 9,name = "Blues")[8],
             effected_volume == 776 ~ brewer.pal(n = 9,name = "Blues")[9],
             effected_volume == 1000 ~ "black")) 


# Protists, 18S
env_18S <- metadata %>% 
  filter(Sequencing_strategy == "MetaB18SV9",
         run_accession %in%
           samples_of_protists_18S_to_rarefy$Sample) %>% 
  select(size_fraction,effected_volume,sample = run_accession)

row.names(env_18S) <- env_18S$sample

# MDS calculation for 18S
mds_prot_18S <-
  metaMDS(t(protists_18S_clean_rarefied))

env_18S <- env_18S %>% 
  mutate(col_volume = 
           case_when(
             effected_volume == 2.5 ~ brewer.pal(n = 9,name = "Blues")[1],
             effected_volume == 10 ~ brewer.pal(n = 9,name = "Blues")[2],
             effected_volume == 30 ~ brewer.pal(n = 9,name = "Blues")[3],
             effected_volume == 60 ~ brewer.pal(n = 9,name = "Blues")[4],
             effected_volume == 100 ~ brewer.pal(n = 9,name = "Blues")[5],
             effected_volume == 120 ~ brewer.pal(n = 9,name = "Blues")[6],
             effected_volume == 496 ~ brewer.pal(n = 9,name = "Blues")[7],
             effected_volume == 716 ~ brewer.pal(n = 9,name = "Blues")[8],
             effected_volume == 776 ~ brewer.pal(n = 9,name = "Blues")[9],
             effected_volume == 1000 ~ "black")) 

# Protist, metagenome
env_metagenome_prot <- metadata %>% 
  filter(Sequencing_strategy == "MetaG",
         run_accession %in% samples_of_protists_metagenomes_to_rarefy$Sample) %>% 
  select(size_fraction,effected_volume,sample = run_accession)

row.names(env_metagenome_prot) <- env_metagenome_prot$sample

# MDS calculation for metagenome
mds_prot_metagenome <-
  metaMDS(t(protists_metagenome_clean_rarefied))

env_metagenome_prot <- env_metagenome_prot %>% 
  mutate(col_volume = 
           case_when(
             effected_volume <= 2.5 ~ brewer.pal(n = 9,name = "Blues")[1],
             effected_volume == 10 ~ brewer.pal(n = 9,name = "Blues")[2],
             effected_volume == 30 ~ brewer.pal(n = 9,name = "Blues")[3],
             effected_volume == 60 ~ brewer.pal(n = 9,name = "Blues")[4],
             effected_volume == 100 ~ brewer.pal(n = 9,name = "Blues")[5],
             effected_volume == 120 ~ brewer.pal(n = 9,name = "Blues")[6],
             effected_volume == 496 ~ brewer.pal(n = 9,name = "Blues")[7],
             effected_volume == 716 ~ brewer.pal(n = 9,name = "Blues")[8],
             effected_volume == 776 ~ brewer.pal(n = 9,name = "Blues")[9],
             effected_volume == 1000 ~ "black")) 

## Plot Figures ##
# Prokaryotes, MetaB16S
plot(mds_prok_16S,
     display = "sites", type="p", main = "Prokaryotes, MetaB16SV4V5")
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

# Prokaryotes, metagenome
plot(mds_prok_metagenome,
     display = "sites", type="p", main = "Prokaryotes, MetaG")
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

## Protists, 18S
plot(mds_prot_18S,
     display = "sites", type="p", main = "Protists, MetaB18SV9")
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
     display = "sites", type="p", main = "Protists, MetaG")
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

## Statistical tests ##

# Function to run permanova and betadisper for 2 variables
full_permanova_summary <- function(y,env,model,a,b,alpha = 0.05,...){
  require(vegan)
  require(dplyr)
  # y is community table (species as columns!)
  # model is the model in formula nomenclature, see: help(formula)
  # a is the first variable used in the model (must be a column name of env)
  # b is the second variable used in the model (must be a column name of env)
  # alpha is the significance level for p-value

  # Make model
  model <- as.formula(model)

  # Permanova for model
  permanova <- 
    adonis2(model, data=env)
  permanova_clean <- 
    permanova %>% 
    as.data.frame() %>% 
    mutate(Test = "PERMANOVA",Term = rownames(permanova))

  # Betadisper helper function
  full_betadisper <- function(var){
   bd <- permutest(
      betadisper(
        vegdist(y),
        group = pull(select(env,all_of(var))),
        type = "centroid"),
      permutations = 999)
   #
   bd[[1]] %>% 
     as.data.frame() %>% 
     mutate(Term = all_of(var),Test = "Betadisper",Term2 = c("Groups","Residuals"))
  }

  # Betadisper for a
  betadisper_a <- full_betadisper(a)
  # Betadisper for b
  betadisper_b <- full_betadisper(b)
  
  # Merge in a single table and add more information
  permanova_clean %>% 
    full_join(betadisper_a) %>% 
    full_join(betadisper_b) %>% 
    mutate(Significance = ifelse(`Pr(>F)` < alpha,
                                 paste0("Significant (p<",alpha,")"),
                                 paste0("Not significant (p>=",alpha,")")))
}

# Permanova for 16S
full_permanova_16S <- 
  full_permanova_summary(y = t(prokaryotes_16S_clean_rarefied),
                       env = env_16S,
                       model = "y ~ effected_volume*size_fraction",
                       a = "size_fraction",
                       b = "effected_volume") %>% 
  mutate(Taxonomic_group = "Prokaryotic",Sequencing_strategy = "MetaB16SV4V5")

# Permanova for 18S
full_permanova_18S <- 
  full_permanova_summary(y = t(protists_18S_clean_rarefied),
                         env = env_18S,
                         model = "y ~ effected_volume*size_fraction",
                         a = "size_fraction",
                         b = "effected_volume") %>% 
  mutate(Taxonomic_group = "Protist",Sequencing_strategy = "MetaB18SV9")

# Permanova for prokaryotes, metagenomes
full_permanova_prok_metagenome <- 
  full_permanova_summary(y = t(prokaryotes_metagenome_clean_rarefied),
                         env = env_metagenome_prok,
                         model = "y ~ effected_volume*size_fraction",
                         a = "size_fraction",
                         b = "effected_volume") %>% 
  mutate(Taxonomic_group = "Prokaryotic",Sequencing_strategy = "MetaG")

# Permanova for protist, metagenomes
full_permanova_prot_metagenome <- 
  full_permanova_summary(y = t(protists_metagenome_clean_rarefied),
                         env = env_metagenome_prot,
                         model = "y ~ effected_volume*size_fraction",
                         a = "size_fraction",
                         b = "effected_volume") %>% 
  mutate(Taxonomic_group = "Protist",Sequencing_strategy = "MetaG")

# Merge all permanova tests
all_permanova_volume_size_fraction <- 
  full_permanova_16S %>% 
  full_join(full_permanova_18S) %>% 
  full_join(full_permanova_prok_metagenome) %>% 
  full_join(full_permanova_prot_metagenome)

## To export
#write.csv(all_permanova_volume_size_fraction, "./output/PERMANOVAs equivalent to betadiversity figures.csv",row.names = FALSE)



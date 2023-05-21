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

###









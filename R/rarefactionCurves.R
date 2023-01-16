rarefaction_lines <- 
  data.frame(Sequencing_strategy = c("MetaB16SV4V5","MetaB18SV9","MetaG","MetaG"),
             Taxonomic_group = c("Prokaryotic","Protist","Prokaryotic","Protist"),
             vline = c(rarefaction_prok_16S,rarefaction_prot_18S,rarefaction_prok_metagenomes,rarefaction_prot_metagenomes))


# Rarefaction figures
### Rarefaction curves by size fractioning
## Prepare data for rarefaction curves
colors_sf_rarecurve <- qualitative_colors[c(1,3,5,7)] #c("#E76BF3","#00B0F6","#00BF7D","yellow1")
colors_volume_rare <- RColorBrewer::brewer.pal(n=9,"PuBu")
#16S
env_all_prok_16S <- metadata %>% 
  filter(Sequencing_strategy == "MetaB16SV4V5",
         run_accession %in% samples_16S$run_accession,
         !is.na(size_fraction)) %>%
  select(size_fraction,effected_volume,device,sample = run_accession) %>% 
  mutate(sf_col = case_when(size_fraction == ">0.22 \U00B5m" ~ colors_sf_rarecurve[1],
                            size_fraction == "0.22-3 \U00B5m" ~ colors_sf_rarecurve[2],
                            size_fraction == "3-20 \U00B5m" ~ colors_sf_rarecurve[3],
                            size_fraction == ">20 \U00B5m" ~ colors_sf_rarecurve[4]),
         vol_col = case_when(effected_volume == 1 ~ colors_volume_rare[1],
                             effected_volume == 2.5 ~ colors_volume_rare[2],
                             effected_volume == 10 ~ colors_volume_rare[3],
                             effected_volume == 30 ~ colors_volume_rare[4],
                             effected_volume == 60 ~ colors_volume_rare[4],
                             effected_volume == 100 ~ colors_volume_rare[5],
                             effected_volume == 120 ~ colors_volume_rare[5],
                             effected_volume == 496 ~ colors_volume_rare[6],
                             effected_volume == 716 ~ colors_volume_rare[7],
                             effected_volume == 776 ~ colors_volume_rare[7],
                             effected_volume == 1000 ~ colors_volume_rare[8]),
         dev_col = case_when(device == "sterivex" ~ "red1",
                             device == "membrane" ~ "grey"),
         dev_lty = case_when(device == "sterivex" ~ "solid",
                             device != "sterivex" ~ "dashed"))
#18S
env_all_prot_18S <- metadata %>% 
  filter(Sequencing_strategy == "MetaB18SV9",
         run_accession %in% samples_18S$run_accession,
         !is.na(size_fraction)) %>%
  select(size_fraction,effected_volume,sample = run_accession,device) %>% 
  mutate(sf_col = case_when(size_fraction == ">0.22 \U00B5m" ~ colors_sf_rarecurve[1],
                            size_fraction == "0.22-3 \U00B5m" ~ colors_sf_rarecurve[2],
                            size_fraction == "3-20 \U00B5m" ~ colors_sf_rarecurve[3],
                            size_fraction == ">20 \U00B5m" ~ colors_sf_rarecurve[4]),
         vol_col = case_when(effected_volume == 1 ~ colors_volume_rare[1],
                             effected_volume == 2.5 ~ colors_volume_rare[2],
                             effected_volume == 10 ~ colors_volume_rare[3],
                             effected_volume == 30 ~ colors_volume_rare[4],
                             effected_volume == 60 ~ colors_volume_rare[4],
                             effected_volume == 100 ~ colors_volume_rare[5],
                             effected_volume == 120 ~ colors_volume_rare[5],
                             effected_volume == 496 ~ colors_volume_rare[6],
                             effected_volume == 716 ~ colors_volume_rare[7],
                             effected_volume == 776 ~ colors_volume_rare[7],
                             effected_volume == 1000 ~ colors_volume_rare[8]),
         dev_col = case_when(device == "sterivex" ~ "red2",
                             device == "membrane" ~ "grey"),
         dev_lty = case_when(device == "sterivex" ~ "solid",
                             device != "sterivex" ~ "dashed"))
#metagenomes - prokaryotes
env_all_metagenomes <- metadata %>% 
  filter(Sequencing_strategy == "MetaG",
         run_accession %in% metagenome_samples$run_accession,
         !is.na(size_fraction)) %>%
  select(size_fraction,effected_volume,sample = run_accession,device) %>% 
  mutate(sf_col = case_when(size_fraction == ">0.22 \U00B5m" ~ colors_sf_rarecurve[1],
                            size_fraction == "0.22-3 \U00B5m" ~ colors_sf_rarecurve[2],
                            size_fraction == "3-20 \U00B5m" ~ colors_sf_rarecurve[3],
                            size_fraction == ">20 \U00B5m" ~ colors_sf_rarecurve[4]),
         vol_col = case_when(effected_volume == 1 ~ colors_volume_rare[1],
                             effected_volume == 2.5 ~ colors_volume_rare[2],
                             effected_volume == 10 ~ colors_volume_rare[3],
                             effected_volume == 30 ~ colors_volume_rare[4],
                             effected_volume == 60 ~ colors_volume_rare[4],
                             effected_volume == 100 ~ colors_volume_rare[5],
                             effected_volume == 120 ~ colors_volume_rare[5],
                             effected_volume == 496 ~ colors_volume_rare[6],
                             effected_volume == 716 ~ colors_volume_rare[7],
                             effected_volume == 776 ~ colors_volume_rare[7],
                             effected_volume == 1000 ~ colors_volume_rare[8]),
         dev_col = case_when(device == "sterivex" ~ "red2",
                             device == "membrane" ~ "grey"),
         dev_lty = case_when(device == "sterivex" ~ "solid",
                             device != "sterivex" ~ "dashed"))

# Rarefaction curves per size fraction
#pdf("./output/Rarefaction curves per size fraction.pdf")
## Make label
plot(c(1), col = NA)
legend(x = "right",
       inset = c(0.1,0.1),
       legend = c(">0.22 µm","0.22-3 µm","3-20 µm",">20 µm"),
       col = colors_sf_rarecurve, 
       lwd = 2.5, 
       box.col = "white", 
       cex = 2,
       title = "Pore size (\U00B5m)", 
       text.font = "Helvetica")
#
par(mfrow = c(2,2))
# 16S
with(env_all_prok_16S, 
     rarecurve(t(prokaryotes_16S_clean),
               col = sf_col, 
               step = 1000, 
               label = FALSE,
               xlab = "Number of reads", 
               ylab  = "Nubmer of OTU", 
               lwd = 2,
               main = "MetaB16SV4V5"))
# 18S
with(env_all_prot_18S, 
     rarecurve(t(protists_18S_clean),
               col = sf_col, step = 1000, 
               label = FALSE,
               xlab = "Number of reads", 
               ylab  = "Number of OTU", 
               lwd = 2,
               main = "MetaB18SV9"))
# MetaG, Prokaryotes
with(env_all_metagenomes, 
     rarecurve(t(prokaryotes_metagenome_clean),
               col = sf_col, 
               step = 100, 
               label = FALSE,
               xlab = "Number of reads",
               ylab  = "Number of OTU",
               lwd = 2,
               main = "MetaG, Prokaryotes"))
# MetaG, Protists
with(env_all_metagenomes, 
     rarecurve(t(protists_metagenome_clean),
               col = sf_col, step = 100, label = FALSE,
               xlab = "Number of reads",
               ylab  = "Number of OTU",
               lwd = 2,
               main = "MetaG, Protists"))
## Volume
par(mfrow = c(2,2))
# 16S
with(env_all_prok_16S, 
     rarecurve(t(prokaryotes_16S_clean),
               col = vol_col, 
               step = 1000, 
               label = FALSE,
               xlab = "Sequences",
               ylab  = "OTUs",
               lwd=2,
               main = "A MetaB16SV4V5"))

# 18S
with(env_all_prot_18S, 
     rarecurve(t(protists_18S_clean),
               col = vol_col, step = 1000, label = FALSE,
               xlab = "Sequences",ylab  = "OTUs",lwd=2,
               main = "B MetaB18SV9"))

# MetaG, Prokaryotes
with(env_all_metagenomes, 
     rarecurve(t(prokaryotes_metagenome_clean),
               col = vol_col, step = 100, label = FALSE,
               xlab = "Sequences",ylab  = "OTUs",lwd=2,
               main = "C MetaG, Prokaryotes"))

# MetaG, Protists
with(env_all_metagenomes, 
     rarecurve(t(protists_metagenome_clean),
               col = vol_col, step = 100, label = FALSE,
               xlab = "Sequences",ylab  = "OTUs",lwd=2,
               main = "D MetaG, Protists")) 

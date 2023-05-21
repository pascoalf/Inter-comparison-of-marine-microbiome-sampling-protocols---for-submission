# Load necessary packages
packages <- c("data.table","readxl","tidyr","stringr","RColorBrewer",
              "vegan","ggplot2","ggpubr","gghalves","gridExtra",
              "patchwork","tibble","rstatix","scales","dplyr", 
              "ggrepel")
lapply(packages, library, character.only = TRUE)

# For reproducibility
set.seed(123)

# Set theme
theme_set(theme_bw())

## Get OTU tables
#load table with OTUs from MGnify v5
emose_mgnify <- fread("https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00001935/pipelines/5.0/file/ERP090011_taxonomy_abundances_SSU_v5.0.tsv") %>% 
  rename(taxa = '#SampleID')

#load table with OTUs from MGnify v4.1
emose_metagenomes_original <- fread("https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00001935/pipelines/4.1/file/ERP090011_taxonomy_abundances_SSU_v4.1.tsv") %>% 
  rename(taxa = '#SampleID')

## Get sample information
# laod ENA data (we need run accession and sample_alias)
ena_reference <- fread("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB87662&result=read_run&fields=study_accession,sample_accession,run_accession,fastq_ftp,sample_alias&format=tsv&download=true&limit=0")

ena_reference <- ena_reference %>% 
  select(run_accession,sample_alias)

## Get detailed information on each sample
## from Supplementary Table S1
sample_registry <- mapply(FUN = read_excel,
                          path="data/Supplementary Table S1.xlsx",
                          sheet = c("EMOSE_SAMPLES_METHODS","SEQUENCING RUNS"),
                          col_names=F)

## get information of sampling methodology of each sample
emose_samples_metadata <- sample_registry[[1]]

## get information related to the sequencing strategies for each sample
emose_seq_info <- sample_registry[[2]]

# Clean sample information table
emose_samples_metadata <- emose_samples_metadata[-c(1:5),-1]

# Name necessary columns
emose_samples_metadata <- emose_samples_metadata %>% 
  select(sample_alias = ...2,protocol = ...11)

## Combine samples methodologies with ENA run accessions
metadata <- ena_reference %>% 
  full_join(emose_samples_metadata, 
            by = "sample_alias")  # Note: there was mistyping error with some sample_alias ID's (that is why I used full_join instead of left_join, 
                                  #but this generated NA's that I dealt in a case by case basis :( ))

## Clean sequencing strategies sample information
emose_seq_info <- emose_seq_info[-c(1:5),-1]
emose_seq_info <- emose_seq_info %>% 
  dplyr::select(sample_alias = ...2,run_accession = ...7,Sequencing_strategy = ...9,
                planned_volume = ...14,effected_volume = ...15,device = ...17,
                method = ...18,replica = ...19,Sequencing_platform = ...25,
                Valid_sequences = ...28)

## Add sequencing strategies data to metadata table by sample_alias and run_accession
metadata <- metadata %>% 
  full_join(emose_seq_info,by=c("sample_alias","run_accession"))

## change obejct class where necessary
metadata$effected_volume <- as.double(metadata$effected_volume)

## Get a column with the size fractions
# based on protocol information
# Pre-cleaning step, because some protocol information is missing for sterivex filters
metadata <- metadata %>% 
  mutate(protocol = case_when(is.na(protocol) & device == "sterivex" ~ "_W>0.22",
                              !is.na(protocol) ~ protocol))
metadata <- separate(metadata, protocol, c("protocol","size_fraction"),"_W")
metadata$size_fraction <- factor(metadata$size_fraction,
                                 levels = c(">0.22","0.22-3","0.22-0.8",
                                            "0.8-20","3-20",">20"),ordered = T,
                                 exclude = NA)
metadata <- metadata %>% mutate(size_fraction = 
                                  factor(paste(size_fraction,"\U00B5m"),
                                         levels = c(">0.22 \U00B5m",
                                                    "0.22-3 \U00B5m",
                                                    "3-20 \U00B5m",
                                                    ">20 \U00B5m")))

## Clean OTU tables

## Remove unwanted taxa
# From MGnify version 5
emose_otus <- emose_mgnify %>% 
  filter(!str_detect(emose_mgnify$taxa,regex(c(".Mitochondria"),dotall=TRUE)),
         !str_detect(emose_mgnify$taxa,regex(c(".Chloroplast"),dotall=TRUE)),
         !str_detect(emose_mgnify$taxa,regex(c(".Metazoa"),dotall=TRUE)),
         !str_detect(emose_mgnify$taxa,regex(c(".Viridiplantae"),dotall=TRUE)), 
         !str_detect(emose_mgnify$taxa,regex(c(".Fungi"),dotall=TRUE)),
         taxa != "Unclassified")

# From MGnify version 4.1 (only metagenomes will be considered afterwards)
emose_metagenomes <- 
  emose_metagenomes_original %>% 
  filter(!str_detect(emose_metagenomes_original$taxa,regex(c(".Mitochondria"),dotall=TRUE)),
         !str_detect(emose_metagenomes_original$taxa,regex(c(".Chloroplast"),dotall=TRUE)),
         !str_detect(emose_metagenomes_original$taxa,regex(c(".Metazoa"),dotall=TRUE)),
         !str_detect(emose_metagenomes_original$taxa,regex(c(".Viridiplantae"),dotall=TRUE)), 
         !str_detect(emose_metagenomes_original$taxa,regex(c(".Fungi"),dotall=TRUE)),
         taxa != c("Unclassified"))

## More data cleaning
metadata$Sequencing_strategy <- 
  recode(metadata$Sequencing_strategy,
         "MetaB16SV4V5 No Sizing" = "MetaB16SV4V5")

## Separate taxonomic groups into different columns
# for MGnify version 5
emose_otus <- emose_otus %>% 
  separate(taxa, c("Super Kingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"),";")

# for MGnify version 4.1 (only metagenomes)
emose_metagenomes <- emose_metagenomes %>% 
  separate(taxa, c("Super Kingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"),";")

## Separate prokaryotes and protists into distinct sub-datasets
# prokaryotes from MGnify version 5
prok_otus <- emose_otus %>% 
  filter(`Super Kingdom` %in% c("sk__Archaea","sk__Bacteria"))
# protists from MGnify version 5
prot_otus <- emose_otus %>% 
  filter(`Super Kingdom` %in% c("sk__Eukaryota"))
# prokaryotes from Mgnify version 4.1 (only metagenomes)
prokaryotic_meta  <- emose_metagenomes %>%
  filter(`Super Kingdom` %in% c("sk__Archaea","sk__Bacteria"))
# protists from Mgnify version 4.1 (only metagenomes)
protist_meta <- emose_metagenomes %>%
  filter(`Super Kingdom` %in% c("sk__Eukaryota"))

## Taxonomy curation
## Protists
# 18S, curated by Roberta:
curated_taxonomy_protist_18S <- 
  read_xlsx("./data/euka_emose_roberta.xlsx",sheet = 2)

# Remove NAs from Fake rank
curated_taxonomy_protist_18S_remove_nas <-  
  curated_taxonomy_protist_18S %>% 
  filter(!is.na(FAKE_rank))

# Remove Rhodophyta (pluricelular)
curated_taxonomy_protist_18S_no_pluricelular <- 
  curated_taxonomy_protist_18S_remove_nas %>%  
  filter(Phylum != "p__Rhodophyta")

curated_taxonomy_protists_18S <- curated_taxonomy_protist_18S_no_pluricelular

# Metagenomes, curated by Roberta:
curated_taxonomy_protist_metagenome <- read_xlsx("./data/protist_meta_initial_Rob.xlsx",sheet = 3)
curated_taxonomy_protist_metagenome <- rename(curated_taxonomy_protist_metagenome, FAKE_rank = FAKE_rank_R)

## Prokaryotes
# Make fake rank with all phyla and Proteobacteria class
# 16S
curated_taxonomy_prokaryotes_16S <-
  prok_otus %>% 
  mutate(FAKE_rank =case_when(str_detect(Phylum,"p__Candidatus_Marinimicrobia") ~ "Candidatus Marinimicrobia",
                              str_detect(Phylum,"p__Candidatus") ~ "Candidate Phylum",
                              str_detect(Phylum,"andidate") ~ "Candidate Phylum",
                              str_detect(Phylum,"p__Proteobacteria") ~ Class,
                              !str_detect(Phylum,"p__Proteobacteria") ~ Phylum
                              )) %>% 
  filter(FAKE_rank != "c__", FAKE_rank != "p__",!is.na(FAKE_rank)) #remove NA assignments

# Metagenomes
curated_taxonomy_prokaryotes_metagenome <-
  prokaryotic_meta %>% 
  mutate(FAKE_rank =case_when(str_detect(Phylum,"p__Candidatus") ~ "Candidate Phylum",
                              str_detect(Phylum,"andidate") ~ "Candidate Phylum",
                              str_detect(Phylum,"p__Proteobacteria") ~ Class,
                              !str_detect(Phylum,"p__Proteobacteria") ~ Phylum)) %>% 
  filter(FAKE_rank != "c__", FAKE_rank != "p__",!is.na(FAKE_rank)) #remove NA assignments


## Brief metadata summary
brief_metadata_summary <- 
  metadata %>% 
  filter(Sequencing_strategy %in% c("MetaB16SV4V5","MetaB18SV9","MetaG"),!is.na(device)) %>%
  mutate(pre_filtration = ifelse(method == "filter>filter>filtrate>filter","yes","no")) %>% 
  group_by(Sequencing_strategy,device,effected_volume,method,pre_filtration,size_fraction,planned_volume) %>% 
  summarise(n())

## To export
#write.csv(brief_metadata_summary,"./output/brief_metadata_summary.csv")


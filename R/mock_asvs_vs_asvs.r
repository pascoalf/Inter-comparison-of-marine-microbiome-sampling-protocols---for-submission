### mock analysis --> ASV version
## just stag mock


## Get OTU tables
#load table with OTUs from MGnify v5
emose_mgnify <- fread("https://www.ebi.ac.uk/metagenomics/api/v1/studies/MGYS00001935/pipelines/5.0/file/ERP090011_taxonomy_abundances_SSU_v5.0.tsv") %>% 
  rename(taxa = '#SampleID')

## Get sample information
# laod ENA data (we need run accession and sample_alias)
ena_reference <- fread("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB87662&result=read_run&fields=study_accession,sample_accession,run_accession,fastq_ftp,sample_alias&format=tsv&download=true&limit=0")

ena_reference <- ena_reference %>% 
  select(run_accession,sample_alias)

## Get detailed information on each sample
## from Supplementary Table S1
sample_registry <- mapply(FUN = read_excel,
                          path="./data/Supplementary Table S1.xlsx",
                          sheet = c("EMOSE_SAMPLES_METHODS","SEQUENCING RUNS"),
                          col_names=F)


## get information of sampling methodology of each sample
emose_samples_metadata <- sample_registry[[1]]

## get information related to the sequencing strategies for each sample
emose_seq_info <- sample_registry[[2]]

# Clean sample information table
emose_samples_metadata <- emose_samples_metadata[-c(1:5),-1]

# Name necessary cols
emose_samples_metadata <- emose_samples_metadata %>% 
  select(sample_alias = ...2,protocol = ...11)

## Combine samples methodologies with ENA run accessions
metadata <- ena_reference %>% 
  full_join(emose_samples_metadata, 
            by = "sample_alias")

## Clean sequencing strategies sample information
emose_seq_info <- emose_seq_info[-c(1:5),-1]
emose_seq_info <- emose_seq_info %>% 
  dplyr::select(sample_alias = ...2,run_accession = ...7,Sequencing_strategy = ...9,
                planned_volume = ...14,effected_volume = ...15,device = ...17,
                method = ...18,replica = ...19,Sequencing_platform = ...25)

## Add sequencing strategies data to metadata table by sample_alias and run_accession
metadata <- metadata %>% 
  full_join(emose_seq_info,by=c("sample_alias","run_accession"))


#get stag mock samples ID (16S No Sizing)
#HiSeq 2500 Rapid Machine
mock_stag_16S <- metadata %>% 
  filter(sample_alias == "EMOSE_STAG-MOCK-16S",
         Sequencing_strategy == "MetaB16SV4V5 No Sizing",
         Sequencing_platform == "Hiseq 2500 Rapid") %>% 
  select(run_accession)
mock_stag_16S$run_accession #I got the fastq file with this ID mannually from ENA, but it could be done more automatically

#get path
stag_path <- "./sequences/stag"

# get file names
stag_fnFs <- paste("./data/sequences/stag/",list.files(stag_path)[1],sep="")
stag_fnRs <-  paste("./data/sequences/stag/",list.files(stag_path)[2],sep="")

# Sample name
stag_sample.names <- "stag_mock_observed"

# Make quality profile to check best trimming
purrr::map(list(stag_fnFs,stag_fnRs),plotQualityProfile)

# Place filtered files in filtered/ subdirectory
stag_filtFs <- file.path(stag_path, "filtered", paste0(stag_sample.names, "_F_filt.fastq.gz"))
stag_filtRs <- file.path(stag_path, "filtered", paste0(stag_sample.names, "_R_filt.fastq.gz"))
names(stag_filtFs) <- stag_sample.names
names(stag_filtRs) <- stag_sample.names

# filter and trim
stag_out <- filterAndTrim(stag_fnFs, stag_filtFs, stag_fnRs, stag_filtRs, 
                          truncLen=c(240,240), #trimming
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE,
                          trimLeft = c(19,20))

stag_error_1 <- learnErrors(stag_filtFs,multithread = TRUE)

stag_error_2 <- learnErrors(stag_filtRs,multithread = TRUE)


# get exact sequences
stag_dada_Fs <- dada(stag_filtFs,err = stag_error_1,multithread = TRUE)

stag_dada_Rs <- dada(stag_filtRs,err = stag_error_2,multithread = TRUE)

stag_mergers <- mergePairs(stag_dada_Fs,stag_filtFs,stag_dada_Rs,stag_filtRs, verbose=TRUE)

# make sequence table
stag_seqtab <- makeSequenceTable(stag_mergers)

# remove chimeras
stag_seqtab.nochim <- removeBimeraDenovo(stag_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


stag_mock_asvs <- colnames(stag_seqtab.nochim)

# Taxonomic assignment with exact matching against the custom database
stag_species <- assignSpecies(stag_seqtab.nochim, "16_custom_reference")

# get expected values from Parada et al.
#expected_abundance <- read.table("./mock_expected.txt",sep = "\t",header=T) #includes stag, even and Parada et al results.


#prepare abundance table

stag_asv <- stag_seqtab.nochim %>% t() %>% as.data.frame()
stag_asv$sequence <- row.names(stag_asv)
row.names(stag_asv) <- NULL
stag_asv$ASV <- paste0("ASV",seq_along(stag_asv[,1]))
stag_asv <- rename(stag_asv,Observed_abundance = "V1")

# Prepare species table

stag_species_observed <- stag_species %>% as.data.frame()
stag_species_observed$sequence <- row.names(stag_species_observed)
row.names(stag_species_observed) <- NULL

# join abundance table with species table

stag_asv_species <- stag_asv %>% 
  left_join(stag_species_observed,by="sequence")

# join expected values with observed values by colony names (Upper case!)

stag_asv_observed_and_expected <- stag_asv_species %>% 
  full_join(expected_abundance,by=c("Species"="Colony_name_upper_case"))

# get observed abundance in relative abundance

stag_asv_observed_and_expected <-
  stag_asv_observed_and_expected %>% 
  filter(!is.na(Observed_abundance)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance))

stag_asv_observed_and_expected %>% 
  ggplot(aes(x=Staggered_expected,y=Observed_abundance_ra))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  stat_regline_equation(label.y=14,aes(label = ..eq.label..))+
  stat_regline_equation(label.y=13,aes(label = ..rr.label..))+
  theme_pubclean()+
  labs(title="Observed abundance vs Staggered expected",
       x = "Staggered expected abundance (%)",
       y = "Observed abundance (%)")

stag_obs_vs_exp_lm <- lm(Observed_abundance_ra ~ Staggered_expected, stag_asv_observed_and_expected)
summary(stag_obs_vs_exp_lm)

stag_asv_observed_and_expected %>% 
  ggplot(aes(x=Staggered_expected,y=Observed_abundance_ra))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  stat_regline_equation(label.y=14,aes(label = ..eq.label..))+
  stat_regline_equation(label.y=13,aes(label = ..rr.label..))+
  scale_y_log10()+
  scale_x_log10()+
  theme_pubclean()+
  labs(title="Observed abundance vs Staggered expected (Log10 transformed)",
       x = "Staggered expected abundance (%)",
       y = "Observed abundance (%)")

#observed vs expected (same as Parada)
stag_asv_observed_and_expected %>% 
  gather(key = "Sample",value ="Relative_abundance",Observed_abundance_ra,Observed_Parada,Staggered_expected) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Expected vs Observed abundance, STAG mock \n with NA assignments",
       x= "Common name",
       y="Relative abundance",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")

## To export
#pdf("Taxonomy observed vs stag_mock vs parada.pdf")

## removed NA assignments
stag_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_STAG_MOCK,Parada_et_al_2015,Staggered_expected) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Expected vs Observed abundance, STAG mock \n removed NA assignments",
       x= "Colony name",
       y="Relative abundance (%)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")

##log transformed
stag_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_STAG_MOCK,Parada_et_al_2015,Staggered_expected) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Expected vs Observed abundance, STAG mock \n removed NA assignments",
       x= "Colony name",
       y="Relative abundance (%) (Log10 transformed)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")+
  scale_y_log10()
dev.off()

stag_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  ggplot(aes(x=Staggered_expected,y=EMOSE_STAG_MOCK))+
  geom_smooth(method="lm",se=F)+
  geom_point()+
  geom_abline(lty=2)+
  stat_regline_equation(label.x=27,label.y=25,aes(label = ..eq.label..))+
  stat_regline_equation(label.x=27,label.y=24,aes(label = ..rr.label..))+
  theme_pubclean()

## end ##


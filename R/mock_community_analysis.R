##Available positive control for 16S PCR

### 16S (No Sizing)
## prokaryotic stag mock
#HiSeq -> done
#MiSeq -> done
## prokaryotic even mock
#HiSeq-> removed
#MiSeq -> removed
## eukaryotic stag mock (???)
#HiSeq-> removed
#MiSeq-> removed
## eukaryotic even mock (???)
#HiSeq-> removed
#MiSeq -> removed

## Packages
packages <- c("tidyr","data.table","stringr","ggplot2","ggpubr","gridExtra",
              "tibble","readxl","dada2","dplyr")
lapply(packages,library,character.only=T)

##### Get data
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


### Prokaryotic stag mock HiSeq (main)

#get stag mock samples ID (16S No Sizing)
#HiSeq 2500 Rapid Machine
mock_stag_16S_prok_hiseq <- metadata %>% 
  filter(sample_alias == "EMOSE_STAG-MOCK-16S",
         Sequencing_strategy == "MetaB16SV4V5 No Sizing",
         Sequencing_platform == "Hiseq 2500 Rapid") %>% 
  select(run_accession)
mock_stag_16S_prok_hiseq$run_accession #I got the fastq file with this ID mannually from ENA, but it could be done more automatically

#get path
stag_prok_hiseq_path <- "./sequences/stag"

# get file names
stag_prok_hiseq_fnFs <- paste("./sequences/stag/",mock_stag_16S_prok_hiseq$run_accession,"_1.fastq",sep="")
stag_prok_hiseq_fnRs <-  paste("./sequences/stag/",mock_stag_16S_prok_hiseq$run_accession,"_2.fastq",sep="")

# Sample name
stag_prok_hiseq_sample.names <- "stag_mock_prok_hiseq_observed"

# Make quality profile to check best trimming
#purrr::map(list(stag_prok_hiseq_fnFs,stag_prok_hiseq_fnRs),plotQualityProfile)

# Place filtered files in filtered/ subdirectory
stag_prok_hiseq_filtFs <- file.path(stag_prok_hiseq_path, "filtered", paste0(stag_prok_hiseq_sample.names, "_F_filt.fastq.gz"))
stag_prok_hiseq_filtRs <- file.path(stag_prok_hiseq_path, "filtered", paste0(stag_prok_hiseq_sample.names, "_R_filt.fastq.gz"))
names(stag_prok_hiseq_filtFs) <- stag_prok_hiseq_sample.names
names(stag_prok_hiseq_filtRs) <- stag_prok_hiseq_sample.names

# filter and trim
stag_prok_hiseq_out <- filterAndTrim(stag_prok_hiseq_fnFs,
                                     stag_prok_hiseq_filtFs, 
                                     stag_prok_hiseq_fnRs, 
                                     stag_prok_hiseq_filtRs, 
                          truncLen=c(240,240), #trimming
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE,
                          trimLeft = c(19,20))

stag_prok_hiseq_error_1 <- learnErrors(stag_prok_hiseq_filtFs,multithread = TRUE)

stag_prok_hiseq_error_2 <- learnErrors(stag_prok_hiseq_filtRs,multithread = TRUE)


# get exact sequences
stag_prok_hiseq_dada_Fs <- dada(stag_prok_hiseq_filtFs,err = stag_prok_hiseq_error_1,multithread = TRUE)

stag_prok_hiseq_dada_Rs <- dada(stag_prok_hiseq_filtRs,err = stag_prok_hiseq_error_2,multithread = TRUE)

stag_prok_hiseq_mergers <- 
  mergePairs(stag_prok_hiseq_dada_Fs,stag_prok_hiseq_filtFs,
             stag_prok_hiseq_dada_Rs,stag_prok_hiseq_filtRs, verbose=TRUE)

# make sequence table
stag_prok_hiseq_seqtab <- makeSequenceTable(stag_prok_hiseq_mergers)

# remove chimeras
stag_prok_hiseq_seqtab.nochim <- 
  removeBimeraDenovo(stag_prok_hiseq_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


stag_prok_hiseq_mock_asvs <- colnames(stag_prok_hiseq_seqtab.nochim)

# Taxonomic assignment with exact matching against the custom database
stag_prok_hiseq_species <- assignSpecies(stag_prok_hiseq_seqtab.nochim, "./prok_stag_custom_reference")

# get expected values from Parada et al.
expected_abundance <- read.table("./mock_expected.txt",sep = "\t",header=T) #includes stag, even and Parada et al results.
expected_abundance$expected_even <- 9.091
#prepare abundance table

stag_prok_hiseq_asv <- stag_prok_hiseq_seqtab.nochim %>% t() %>% as.data.frame()
stag_prok_hiseq_asv$sequence <- row.names(stag_prok_hiseq_asv)
row.names(stag_prok_hiseq_asv) <- NULL
stag_prok_hiseq_asv$ASV <- paste0("ASV",seq_along(stag_prok_hiseq_asv[,1]))
stag_prok_hiseq_asv <- rename(stag_prok_hiseq_asv,Observed_abundance = "V1")

# Prepare species table

stag_prok_hiseq_species_observed <- stag_prok_hiseq_species %>% as.data.frame()
stag_prok_hiseq_species_observed$sequence <- row.names(stag_prok_hiseq_species_observed)
row.names(stag_prok_hiseq_species_observed) <- NULL

# join abundance table with species table

stag_prok_hiseq_asv_species <- stag_prok_hiseq_asv %>% 
  left_join(stag_prok_hiseq_species_observed,by="sequence")

# join expected values with observed values by colony names (Upper case!)

stag_prok_hiseq_asv_observed_and_expected <- stag_prok_hiseq_asv_species %>% 
  full_join(expected_abundance,by=c("Species"="Colony_name_upper_case"))

stag_prok_hiseq_asv_observed_and_expected <-
  stag_prok_hiseq_asv_observed_and_expected %>% 
  filter(!is.na(Observed_abundance)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance))

############################################################

# get observed abundance in relative abundance

pdf("Positive control results: Prokaryotes STAG Mock HiSeq.pdf")

stag_prok_hiseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_STAG_MOCK,Staggered_expected) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Taxonomy: Expected vs Observed abundance \nProkaryotes stag mock community \n(removed NA assignments) \n(HiSeq 2500 Rapid)",
       x= "Colony name",
       y="Relative abundance (%)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")

##log transformed
stag_prok_hiseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_STAG_MOCK,Staggered_expected) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Taxonomy: Expected vs Observed abundance \nProkaryotes stag mock community \n(removed NA assignments) \n(Log10 transformed to see lower abundance)\n(HiSeq 2500 Rapid)",
       x= "Colony name",
       y="Relative abundance (%) (Log10 transformed)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")+
  scale_y_log10()

stag_prok_hiseq_asv_observed_and_expected %>% 
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
  theme_pubclean()+
  labs(title="Relative abundance: Expected vs Observed \nProkaryotes stag mock community \n(removed NA assignments) \n(HiSeq 2500 Rapid)")

dev.off()

###### Prokaryotes STAG MiSeq

#get stag mock samples ID (16S No Sizing)
#miseq 2500 Rapid Machine
mock_stag_16S_prok_miseq <- metadata %>% 
  filter(sample_alias == "EMOSE_STAG-MOCK-16S",
         Sequencing_strategy == "MetaB16SV4V5 No Sizing",
         Sequencing_platform == "MiSeq") %>% 
  select(run_accession)
mock_stag_16S_prok_miseq$run_accession #I got the fastq file with this ID mannually from ENA, but it could be done more automatically

#get path
stag_prok_miseq_path <- "./sequences/stag"

# get file names
stag_prok_miseq_fnFs <- paste("./sequences/stag/",mock_stag_16S_prok_miseq$run_accession,"_1.fastq",sep="")
stag_prok_miseq_fnRs <-  paste("./sequences/stag/",mock_stag_16S_prok_miseq$run_accession,"_2.fastq",sep="")

# Sample name
stag_prok_miseq_sample.names <- "stag_mock_prok_miseq_observed"

# Make quality profile to check best trimming
#purrr::map(list(stag_prok_miseq_fnFs,stag_prok_miseq_fnRs),plotQualityProfile)

# Place filtered files in filtered/ subdirectory
stag_prok_miseq_filtFs <- file.path(stag_prok_miseq_path, "filtered", paste0(stag_prok_miseq_sample.names, "_F_filt.fastq.gz"))
stag_prok_miseq_filtRs <- file.path(stag_prok_miseq_path, "filtered", paste0(stag_prok_miseq_sample.names, "_R_filt.fastq.gz"))
names(stag_prok_miseq_filtFs) <- stag_prok_miseq_sample.names
names(stag_prok_miseq_filtRs) <- stag_prok_miseq_sample.names

# filter and trim
stag_prok_miseq_out <- filterAndTrim(stag_prok_miseq_fnFs, stag_prok_miseq_filtFs, stag_prok_miseq_fnRs, stag_prok_miseq_filtRs, 
                                     truncLen=c(250,225), #trimming
                                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                     compress=TRUE, multithread=TRUE,
                                     trimLeft = c(19,20))

stag_prok_miseq_error_1 <- learnErrors(stag_prok_miseq_filtFs,multithread = TRUE)

stag_prok_miseq_error_2 <- learnErrors(stag_prok_miseq_filtRs,multithread = TRUE)


# get exact sequences
stag_prok_miseq_dada_Fs <- dada(stag_prok_miseq_filtFs,err = stag_prok_miseq_error_1,multithread = TRUE)

stag_prok_miseq_dada_Rs <- dada(stag_prok_miseq_filtRs,err = stag_prok_miseq_error_2,multithread = TRUE)

stag_prok_miseq_mergers <- 
  mergePairs(stag_prok_miseq_dada_Fs,stag_prok_miseq_filtFs,
             stag_prok_miseq_dada_Rs,stag_prok_miseq_filtRs, verbose=TRUE)

# make sequence table
stag_prok_miseq_seqtab <- makeSequenceTable(stag_prok_miseq_mergers)

# remove chimeras
stag_prok_miseq_seqtab.nochim <- 
  removeBimeraDenovo(stag_prok_miseq_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


stag_prok_miseq_mock_asvs <- colnames(stag_prok_miseq_seqtab.nochim)

# Taxonomic assignment with exact matching against the custom database
stag_prok_miseq_species <- assignSpecies(stag_prok_miseq_seqtab.nochim, "./prok_stag_custom_reference")

#prepare abundance table
stag_prok_miseq_asv <- stag_prok_miseq_seqtab.nochim %>% t() %>% as.data.frame()
stag_prok_miseq_asv$sequence <- row.names(stag_prok_miseq_asv)
row.names(stag_prok_miseq_asv) <- NULL
stag_prok_miseq_asv$ASV <- paste0("ASV",seq_along(stag_prok_miseq_asv[,1]))
stag_prok_miseq_asv <- rename(stag_prok_miseq_asv,Observed_abundance = "V1")

# Prepare species table

stag_prok_miseq_species_observed <- stag_prok_miseq_species %>% as.data.frame()
stag_prok_miseq_species_observed$sequence <- row.names(stag_prok_miseq_species_observed)
row.names(stag_prok_miseq_species_observed) <- NULL

# join abundance table with species table

stag_prok_miseq_asv_species <- stag_prok_miseq_asv %>% 
  left_join(stag_prok_miseq_species_observed,by="sequence")

# join expected values with observed values by colony names (Upper case!)

stag_prok_miseq_asv_observed_and_expected <- stag_prok_miseq_asv_species %>% 
  full_join(expected_abundance,by=c("Species"="Colony_name_upper_case"))

stag_prok_miseq_asv_observed_and_expected <-
  stag_prok_miseq_asv_observed_and_expected %>% 
  filter(!is.na(Observed_abundance)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance))
############################################################

# get observed abundance in relative abundance

pdf("Positive control results: Prokaryotes STAG Mock MiSeq.pdf")

stag_prok_miseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_STAG_MOCK,Staggered_expected) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Taxonomy: Expected vs Observed abundance \nProkaryotes stag mock community \n(removed NA assignments) \n(miseq 2500 Rapid)",
       x= "Colony name",
       y="Relative abundance (%)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")

##log transformed
stag_prok_miseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_STAG_MOCK,Staggered_expected) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Taxonomy: Expected vs Observed abundance \nProkaryotes stag mock community \n(removed NA assignments) \n(Log10 transformed to see lower abundance)\n(miseq 2500 Rapid)",
       x= "Colony name",
       y="Relative abundance (%) (Log10 transformed)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")+
  scale_y_log10()

stag_prok_miseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra") %>% 
  ggplot(aes(x=Staggered_expected,y=EMOSE_STAG_MOCK))+
  geom_smooth(method="lm",se=F)+
  geom_point()+
  geom_abline(lty=2)+
  stat_regline_equation(label.x=27,label.y=25,aes(label = ..eq.label..))+
  stat_regline_equation(label.x=27,label.y=24,aes(label = ..rr.label..))+
  theme_pubclean()+
  labs(title="Relative abundance: Expected vs Observed \nProkaryotes stag mock community \n(removed NA assignments) \n(miseq 2500 Rapid)")

dev.off()

#### Join stag on single figure

stag_prok_hiseq_merged <- stag_prok_hiseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra")

stag_prok_miseq_merged <- stag_prok_miseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_STAG_MOCK = "Observed_abundance_ra")

stag_all_combined <- 
  stag_prok_hiseq_merged %>% full_join(stag_prok_miseq_merged,by=c("Common_name","Staggered_expected"))

pdf("Positive control results for 16S HiSeq and MiSeq")

stag_all_combined %>%
  rename(HiSeq = "EMOSE_STAG_MOCK.x",MiSeq = "EMOSE_STAG_MOCK.y") %>% 
  gather(key="Sequencing_machine",
         value="EMOSE_STAG_MOCK",
         HiSeq,MiSeq) %>% 
  ggplot(aes(Staggered_expected,EMOSE_STAG_MOCK,col=Sequencing_machine))+
  geom_point()+
  geom_smooth(se=FALSE,method="lm")+
  geom_abline(lty=2)+ 
  stat_regline_equation(label.x=27,label.y=c(16,19),aes(label = ..eq.label..))+
  stat_regline_equation(label.x=27,label.y=c(15,18),aes(label = ..rr.label..))+
  labs(title="Relative abundance: Expected vs Observed \nProkaryotes stag mock community \n(removed NA assignments)",
       x="Staggered expected relative abundance (%)",
       y="Observed Staggered Mock relative abundance (%)",
       col="Sequencing machine")+
  scale_color_brewer(palette=2,type="qual")+
  theme_pubclean()

#Taxonomy together

stag_all_combined %>%
  rename(EMOSE_HiSeq = "EMOSE_STAG_MOCK.x",EMOSE_MiSeq = "EMOSE_STAG_MOCK.y") %>% 
  gather(key="Sample",
         value="Relative_abundance",
         EMOSE_HiSeq,EMOSE_MiSeq,Staggered_expected) %>%
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  scale_fill_brewer(palette=3,type="qual")+
  theme_pubclean()+
  labs(title="Relative abundance: Expected vs Observed \nProkaryotes stag mock community \n(removed NA assignments)",
       y="Relative abundance (%)",
       x="Colony name",
       fill="")
  
stag_all_combined %>%
  rename(EMOSE_HiSeq = "EMOSE_STAG_MOCK.x",EMOSE_MiSeq = "EMOSE_STAG_MOCK.y") %>% 
  gather(key="Sample",
         value="Relative_abundance",
         EMOSE_HiSeq,EMOSE_MiSeq,Staggered_expected) %>%
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  scale_fill_brewer(palette=3,type="qual")+
  theme_pubclean()+
  labs(title="Relative abundance: Expected vs Observed \nProkaryotes stag mock community \n(removed NA assignments)",
       y="Relative abundance (%)",
       x="Colony name",
       fill="")+
  scale_y_log10()

stag_all_combined %>% 
  ggplot(aes(Observed_abundance.x*100/sum(Observed_abundance.x),Observed_abundance.y*100/sum(Observed_abundance.y)))+
  geom_point()+
  geom_smooth(se=FALSE,method="lm")+
  theme_pubclean()+  
  geom_abline(lty=2)+ 
  stat_regline_equation(label.x=25,label.y=c(16,19),aes(label = ..eq.label..))+
  stat_regline_equation(label.x=25,label.y=c(15,18),aes(label = ..rr.label..))+
  labs(x="HiSeq 2500 Rapid (Relative abundance %)",
       y="MiSeq (Relative abundance %)",
       title="HiSeq 2500 Rapid vs MiSeq machine \n(NA assignments removed)")

dev.off()
######################################
####below I think will not be used####

#get even mock samples ID (16S No Sizing)
#hiseq 2500 Rapid Machine
mock_even_16S_prok_hiseq <- metadata %>% 
  filter(sample_alias == "EMOSE_EVEN-MOCK-16S",
         Sequencing_strategy == "MetaB16SV4V5 No Sizing",
         Sequencing_platform == "Hiseq 2500 Rapid") %>% 
  select(run_accession)
mock_even_16S_prok_hiseq$run_accession #I got the fastq file with this ID mannually from ENA, but it could be done more automatically

#get path
even_prok_hiseq_path <- "./sequences/even"

# get file names
even_prok_hiseq_fnFs <- paste("./sequences/even/",mock_even_16S_prok_hiseq$run_accession,"_1.fastq",sep="")
even_prok_hiseq_fnRs <-  paste("./sequences/even/",mock_even_16S_prok_hiseq$run_accession,"_2.fastq",sep="")

# Sample name
even_prok_hiseq_sample.names <- "even_mock_prok_hiseq_observed"

# Make quality profile to check best trimming
purrr::map(list(even_prok_hiseq_fnFs,even_prok_hiseq_fnRs),plotQualityProfile)

# Place filtered files in filtered/ subdirectory
even_prok_hiseq_filtFs <- file.path(even_prok_hiseq_path, "filtered", paste0(even_prok_hiseq_sample.names, "_F_filt.fastq.gz"))
even_prok_hiseq_filtRs <- file.path(even_prok_hiseq_path, "filtered", paste0(even_prok_hiseq_sample.names, "_R_filt.fastq.gz"))
names(even_prok_hiseq_filtFs) <- even_prok_hiseq_sample.names
names(even_prok_hiseq_filtRs) <- even_prok_hiseq_sample.names

# filter and trim
even_prok_hiseq_out <- filterAndTrim(even_prok_hiseq_fnFs, even_prok_hiseq_filtFs, even_prok_hiseq_fnRs, even_prok_hiseq_filtRs, 
                                     truncLen=c(240,240), #trimming
                                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                     compress=TRUE, multithread=TRUE,
                                     trimLeft = c(19,20))

even_prok_hiseq_error_1 <- learnErrors(even_prok_hiseq_filtFs,multithread = TRUE)

even_prok_hiseq_error_2 <- learnErrors(even_prok_hiseq_filtRs,multithread = TRUE)


# get exact sequences
even_prok_hiseq_dada_Fs <- dada(even_prok_hiseq_filtFs,err = even_prok_hiseq_error_1,multithread = TRUE)

even_prok_hiseq_dada_Rs <- dada(even_prok_hiseq_filtRs,err = even_prok_hiseq_error_2,multithread = TRUE)

even_prok_hiseq_mergers <- 
  mergePairs(even_prok_hiseq_dada_Fs,even_prok_hiseq_filtFs,
             even_prok_hiseq_dada_Rs,even_prok_hiseq_filtRs, verbose=TRUE)

# make sequence table
even_prok_hiseq_seqtab <- makeSequenceTable(even_prok_hiseq_mergers)

# remove chimeras
even_prok_hiseq_seqtab.nochim <- 
  removeBimeraDenovo(even_prok_hiseq_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


even_prok_hiseq_mock_asvs <- colnames(even_prok_hiseq_seqtab.nochim)

# Taxonomic assignment with exact matching against the custom database
even_prok_hiseq_species <- assignSpecies(even_prok_hiseq_seqtab.nochim, "./prok_even_custom_reference")

# get expected values from Parada et al.
expected_abundance <- read.table("./mock_expected.txt",sep = "\t",header=T) #includes even, even and Parada et al results.

#prepare abundance table

even_prok_hiseq_asv <- even_prok_hiseq_seqtab.nochim %>% t() %>% as.data.frame()
even_prok_hiseq_asv$sequence <- row.names(even_prok_hiseq_asv)
row.names(even_prok_hiseq_asv) <- NULL
even_prok_hiseq_asv$ASV <- paste0("ASV",seq_along(even_prok_hiseq_asv[,1]))
even_prok_hiseq_asv <- rename(even_prok_hiseq_asv,Observed_abundance = "V1")

# Prepare species table

even_prok_hiseq_species_observed <- even_prok_hiseq_species %>% as.data.frame()
even_prok_hiseq_species_observed$sequence <- row.names(even_prok_hiseq_species_observed)
row.names(even_prok_hiseq_species_observed) <- NULL

# join abundance table with species table

even_prok_hiseq_asv_species <- even_prok_hiseq_asv %>% 
  left_join(even_prok_hiseq_species_observed,by="sequence")

# join expected values with observed values by colony names (Upper case!)

even_prok_hiseq_asv_observed_and_expected <- even_prok_hiseq_asv_species %>% 
  full_join(expected_abundance,by=c("Species"="Colony_name_upper_case"))

even_prok_hiseq_asv_observed_and_expected <-
  even_prok_hiseq_asv_observed_and_expected %>% 
  filter(!is.na(Observed_abundance)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance))
############################################################

# get observed abundance in relative abundance

pdf("Positive control results: Prokaryotes even Mock HiSeq.pdf")

even_prok_hiseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_even_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_even_MOCK,expected_even) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Taxonomy: Expected vs Observed abundance \nProkaryotes even mock community \n(removed NA assignments) \n(hiseq 2500 Rapid)",
       x= "Colony name",
       y="Relative abundance (%)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")

dev.off()

#### Porkaryotes even mock MiSeq
#get even mock samples ID (16S No Sizing)
#miseq 2500 Rapid Machine
mock_even_16S_prok_miseq <- metadata %>% 
  filter(sample_alias == "EMOSE_EVEN-MOCK-16S",
         Sequencing_strategy == "MetaB16SV4V5 No Sizing",
         Sequencing_platform == "MiSeq") %>% 
  select(run_accession)
mock_even_16S_prok_miseq$run_accession #I got the fastq file with this ID mannually from ENA, but it could be done more automatically

#get path
even_prok_miseq_path <- "./sequences/even"

# get file names
even_prok_miseq_fnFs <- paste("./sequences/even/",mock_even_16S_prok_miseq$run_accession,"_1.fastq",sep="")
even_prok_miseq_fnRs <-  paste("./sequences/even/",mock_even_16S_prok_miseq$run_accession,"_2.fastq",sep="")

# Sample name
even_prok_miseq_sample.names <- "even_mock_prok_miseq_observed"

# Make quality profile to check best trimming
purrr::map(list(even_prok_miseq_fnFs,even_prok_miseq_fnRs),plotQualityProfile)

# Place filtered files in filtered/ subdirectory
even_prok_miseq_filtFs <- file.path(even_prok_miseq_path, "filtered", paste0(even_prok_miseq_sample.names, "_F_filt.fastq.gz"))
even_prok_miseq_filtRs <- file.path(even_prok_miseq_path, "filtered", paste0(even_prok_miseq_sample.names, "_R_filt.fastq.gz"))
names(even_prok_miseq_filtFs) <- even_prok_miseq_sample.names
names(even_prok_miseq_filtRs) <- even_prok_miseq_sample.names

# filter and trim
even_prok_miseq_out <- filterAndTrim(even_prok_miseq_fnFs, even_prok_miseq_filtFs, even_prok_miseq_fnRs, even_prok_miseq_filtRs, 
                                     truncLen=c(250,225), #trimming
                                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                                     compress=TRUE, multithread=TRUE,
                                     trimLeft = c(19,20))

even_prok_miseq_error_1 <- learnErrors(even_prok_miseq_filtFs,multithread = TRUE)

even_prok_miseq_error_2 <- learnErrors(even_prok_miseq_filtRs,multithread = TRUE)


# get exact sequences
even_prok_miseq_dada_Fs <- dada(even_prok_miseq_filtFs,err = even_prok_miseq_error_1,multithread = TRUE)

even_prok_miseq_dada_Rs <- dada(even_prok_miseq_filtRs,err = even_prok_miseq_error_2,multithread = TRUE)

even_prok_miseq_mergers <- 
  mergePairs(even_prok_miseq_dada_Fs,even_prok_miseq_filtFs,
             even_prok_miseq_dada_Rs,even_prok_miseq_filtRs, verbose=TRUE)

# make sequence table
even_prok_miseq_seqtab <- makeSequenceTable(even_prok_miseq_mergers)

# remove chimeras
even_prok_miseq_seqtab.nochim <- 
  removeBimeraDenovo(even_prok_miseq_seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


even_prok_miseq_mock_asvs <- colnames(even_prok_miseq_seqtab.nochim)

# Taxonomic assignment with exact matching against the custom database
even_prok_miseq_species <- assignSpecies(even_prok_miseq_seqtab.nochim, "./prok_even_custom_reference")

#prepare abundance table

even_prok_miseq_asv <- even_prok_miseq_seqtab.nochim %>% t() %>% as.data.frame()
even_prok_miseq_asv$sequence <- row.names(even_prok_miseq_asv)
row.names(even_prok_miseq_asv) <- NULL
even_prok_miseq_asv$ASV <- paste0("ASV",seq_along(even_prok_miseq_asv[,1]))
even_prok_miseq_asv <- rename(even_prok_miseq_asv,Observed_abundance = "V1")

# Prepare species table

even_prok_miseq_species_observed <- even_prok_miseq_species %>% as.data.frame()
even_prok_miseq_species_observed$sequence <- row.names(even_prok_miseq_species_observed)
row.names(even_prok_miseq_species_observed) <- NULL

# join abundance table with species table

even_prok_miseq_asv_species <- even_prok_miseq_asv %>% 
  left_join(even_prok_miseq_species_observed,by="sequence")

# join expected values with observed values by colony names (Upper case!)

even_prok_miseq_asv_observed_and_expected <- even_prok_miseq_asv_species %>% 
  full_join(expected_abundance,by=c("Species"="Colony_name_upper_case"))

even_prok_miseq_asv_observed_and_expected <-
  even_prok_miseq_asv_observed_and_expected %>% 
  filter(!is.na(Observed_abundance)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance))
############################################################

# get observed abundance in relative abundance

pdf("Positive control results: Prokaryotes even Mock MiSeq.pdf")

even_prok_miseq_asv_observed_and_expected %>% 
  filter(!is.na(Species)) %>% 
  mutate(Observed_abundance_ra = Observed_abundance*100/sum(Observed_abundance)) %>% 
  rename(EMOSE_even_MOCK = "Observed_abundance_ra",
         Parada_et_al_2015 = "Observed_Parada") %>% 
  gather(key = "Sample",value ="Relative_abundance",EMOSE_even_MOCK,expected_even) %>% 
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge",col="black")+
  coord_flip()+
  theme_pubclean()+
  labs(title="Taxonomy: Expected vs Observed abundance \nProkaryotes even mock community \n(removed NA assignments) \n(miseq 2500 Rapid)",
       x= "Colony name",
       y="Relative abundance (%)",
       fill = "From:")+
  scale_fill_brewer(palette=4,type="div")

dev.off()

###
### Arrange

#
qualitative_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

gridExtra::grid.arrange(
stag_all_combined %>%
  rename(HiSeq = "EMOSE_STAG_MOCK.x",
         MiSeq = "EMOSE_STAG_MOCK.y") %>%
  gather(key="Sequencing_machine",
         value="EMOSE_STAG_MOCK",
         HiSeq,MiSeq) %>%
  ggplot(aes(Staggered_expected,EMOSE_STAG_MOCK,col=Sequencing_machine))+
  geom_point()+
  geom_smooth(se=FALSE,method="lm")+
  geom_abline(lty=2)+
  stat_regline_equation(label.x=26,label.y=c(14,18),aes(label = ..eq.label..))+
  stat_regline_equation(label.x=26,label.y=c(12,16),aes(label = ..rr.label..))+
  labs(title="a",
       x="Expected relative abundance (%)",
       y="Observed relative abundance (%)",
       col="Sequencing machine:")+
  scale_color_manual(values = qualitative_colors[c(1,2)])+
  theme_pubclean()+
  theme(text = element_text(size = 12, family = "Helvetica"))

,

stag_all_combined %>%
  rename(HiSeq = "EMOSE_STAG_MOCK.x",
         MiSeq = "EMOSE_STAG_MOCK.y",
         Expected = "Staggered_expected") %>%
  gather(key="Sample",
         value="Relative_abundance",
         HiSeq, MiSeq, Expected) %>%
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge") + #,col="black")+
  coord_flip()+
  scale_fill_manual(values = qualitative_colors[c(3, 1, 2)])+
  theme_pubclean()+
  labs(title="b",
       y="Relative abundance (%)",
       x="Colony name",
       fill="")+
  theme(text = element_text(size = 12, family = "Helvetica"))

,

stag_all_combined %>%
  rename(HiSeq = "EMOSE_STAG_MOCK.x",
         MiSeq = "EMOSE_STAG_MOCK.y",
         Expected = "Staggered_expected") %>%
  gather(key="Sample",
         value="Relative_abundance",
         HiSeq,MiSeq,Expected) %>%
  ggplot(aes(x=reorder(Common_name,-Relative_abundance),y=Relative_abundance,fill=Sample))+
  geom_col(position="dodge") + #,col="black")+
  coord_flip()+
  scale_fill_manual(values = qualitative_colors[c(3, 1, 2)])+
  theme_pubclean()+
  labs(title="c",
       y="Relative abundance (%) in Log10 scale",
       x="Colony name",
       fill="")+
  scale_y_log10() +
  theme(text = element_text(size = 12, family = "Helvetica"))

,

stag_all_combined %>%
  ggplot(aes(Observed_abundance.x*100/sum(Observed_abundance.x),Observed_abundance.y*100/sum(Observed_abundance.y)))+
  geom_point()+
  geom_smooth(se=FALSE,method="lm")+
  theme_pubclean()+
  geom_abline(lty=2)+
  stat_regline_equation(label.x=24,label.y=c(16,19),aes(label = ..eq.label..))+
  stat_regline_equation(label.x=24,label.y=c(14,17),aes(label = ..rr.label..))+
  labs(x="HiSeq 2500 Rapid (relative abundance %)",
       y="MiSeq (relative abundance %)",
       title="d")+
  theme(text = element_text(size = 12, family = "Helvetiva"))
)


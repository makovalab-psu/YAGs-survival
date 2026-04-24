library(tidyr)
library(ggplot2)
library(ggpubr)
library(microseq)
library(purrr)

#from 'Attributes' column of a gff/gtf file create separate columns; returns a dataframe
gff_attr <- function(gff_file) {
  
  gff_file_df <- data.frame(gff_file)
  gff_file_inter <- strsplit(gsub("\"", "", gff_file_df$Attributes), split=";")
  gff_file_inter_no_space <- lapply(gff_file_inter, trimws)
  gff_file_pairs_list <- lapply(gff_file_inter_no_space, strsplit, split=' ')
  
  df_attributes <- data.frame()
  
  for (pair in gff_file_pairs_list){
    named_list_iter <- lapply(pair, function(sub) setNames(list(sub[[2]]), sub[[1]]))
    df_iter <- data.frame(t(data.frame(unlist(named_list_iter))))
    df_attributes <- bind_rows(df_attributes, df_iter)
  }
  
  rownames(df_attributes) <- NULL
  
  gff_file_attributes <- cbind(gff_file_df[, 1:8], df_attributes)
  return(gff_file_attributes)
  
}

#read gff files
gtf_annot_final_bonobo <- readGFF("C:/Users/gresh/Downloads/temp/source/bonobo_techrep1.stringtie.annotation.hc_cORF.Ychr.source.gff")
gtf_annot_final_chimp <- readGFF("C:/Users/gresh/Downloads/temp/source/chimp_techrep1.stringtie.annotation.hc_cORF.Ychr.source.gff")
gtf_annot_final_gor <- readGFF("C:/Users/gresh/Downloads/temp/source/gor_techrep1.stringtie.annotation.hc_cORF.Ychr.source.gff")
gtf_annot_final_bor_orang <- readGFF("C:/Users/gresh/Downloads/temp/source/bor_orang_techrep1.stringtie.annotation.hc_cORF.Ychr.source.gff")
gtf_annot_final_sum_orang <- readGFF("C:/Users/gresh/Downloads/temp/source/sum_orang_techrep1.stringtie.annotation.hc_cORF.Ychr.source.gff")
gtf_annot_final_human <- readGFF("C:/Users/gresh/Downloads/temp/source/human_techrep1.stringtie.annotation.hc_cORF.Ychr.source.gff")

#remove empty rows
gtf_annot_final_bonobo <- gtf_annot_final_bonobo %>% filter(Start != 'NA')
gtf_annot_final_chimp <- gtf_annot_final_chimp %>% filter(Seqid != 'NA')
gtf_annot_final_gor <- gtf_annot_final_gor %>% filter(Seqid != 'NA')
gtf_annot_final_bor_orang <- gtf_annot_final_bor_orang %>% filter(Seqid != 'NA')
gtf_annot_final_sum_orang <- gtf_annot_final_sum_orang %>% filter(Seqid != 'NA')
gtf_annot_final_human <- gtf_annot_final_human %>% filter(Seqid != 'NA')

#make each attribute into a separate column
gtf_annot_final_bonobo <- gff_attr(gtf_annot_final_bonobo)
gtf_annot_final_chimp <- gff_attr(gtf_annot_final_chimp)
gtf_annot_final_gor <- gff_attr(gtf_annot_final_gor)
gtf_annot_final_bor_orang <- gff_attr(gtf_annot_final_bor_orang)
gtf_annot_final_sum_orang <- gff_attr(gtf_annot_final_sum_orang)
gtf_annot_final_human <- gff_attr(gtf_annot_final_human)

#add species column
gtf_annot_final_bonobo$species <- 'bonobo'
gtf_annot_final_chimp$species <- 'chimp'
gtf_annot_final_gor$species <- 'gor'
gtf_annot_final_bor_orang$species <- 'bor_orang'
gtf_annot_final_sum_orang$species <- 'sum_orang'
gtf_annot_final_human$species <- 'human'

################################################################################
#remove duplicate transcripts - currently the column sseqid has unique inputs
#but column transcript_id has duplicated inputs

# unique_transcripts_hcORF_bonobo <- gtf_annot_final_bonobo %>% 
#   filter(Type == 'transcript') %>% 
#   filter(sseqid != 'NA') %>% 
#   group_by(transcript_id) %>% 
#   filter(length == max(length)) %>% 
#   filter(pident == max(pident)) %>% 
#   filter(qcovs == max(qcovs)) %>% 
#   as.data.frame() 
# 
# gtf_annot_final_bonobo_with_hcORF <-
#   gtf_annot_final_bonobo %>% 
#   filter(sseqid %in% unique_transcripts_hcORF_bonobo$qseqid)
# 
# gtf_annot_final_bonobo_no_hcORF <-
#   gtf_annot_final_bonobo %>% 
#   filter(sseqid == 'NA')
# 
# gtf_annot_final_bonobo <-
#   rbind(gtf_annot_final_bonobo_with_hcORF, gtf_annot_final_bonobo_no_hcORF)
# 
# unique_transcripts_hcORF_chimp <- gtf_annot_final_chimp %>% 
#   filter(Type == 'transcript') %>% 
#   filter(sseqid != 'NA') %>% 
#   group_by(transcript_id) %>% 
#   filter(length == max(length)) %>% 
#   filter(pident == max(pident)) %>% 
#   filter(qcovs == max(qcovs)) %>% 
#   as.data.frame() 
# 
# gtf_annot_final_chimp_with_hcORF <-
#   gtf_annot_final_chimp %>% 
#   filter(sseqid %in% unique_transcripts_hcORF_chimp$qseqid)
# 
# gtf_annot_final_chimp_no_hcORF <-
#   gtf_annot_final_chimp %>% 
#   filter(sseqid == 'NA')
# 
# gtf_annot_final_chimp <-
#   rbind(gtf_annot_final_chimp_with_hcORF, gtf_annot_final_chimp_no_hcORF)
# 
# unique_transcripts_hcORF_gor <- gtf_annot_final_gor %>% 
#   filter(Type == 'transcript') %>% 
#   filter(sseqid != 'NA') %>% 
#   group_by(transcript_id) %>% 
#   filter(length == max(length)) %>% 
#   filter(pident == max(pident)) %>% 
#   filter(qcovs == max(qcovs)) %>% 
#   as.data.frame() 
# 
# gtf_annot_final_gor_with_hcORF <-
#   gtf_annot_final_gor %>% 
#   filter(sseqid %in% unique_transcripts_hcORF_gor$qseqid)
# 
# gtf_annot_final_gor_no_hcORF <-
#   gtf_annot_final_gor %>% 
#   filter(sseqid == 'NA')
# 
# gtf_annot_final_gor <-
#   rbind(gtf_annot_final_gor_with_hcORF, gtf_annot_final_gor_no_hcORF)
# 
# unique_transcripts_hcORF_bor_orang <- gtf_annot_final_bor_orang %>% 
#   filter(Type == 'transcript') %>% 
#   filter(sseqid != 'NA') %>% 
#   group_by(transcript_id) %>% 
#   filter(length == max(length)) %>% 
#   filter(pident == max(pident)) %>% 
#   filter(qcovs == max(qcovs)) %>% 
#   as.data.frame() 
# 
# gtf_annot_final_bor_orang_with_hcORF <-
#   gtf_annot_final_bor_orang %>% 
#   filter(sseqid %in% unique_transcripts_hcORF_bor_orang$qseqid)
# 
# gtf_annot_final_bor_orang_no_hcORF <-
#   gtf_annot_final_bor_orang %>% 
#   filter(sseqid == 'NA')
# 
# gtf_annot_final_bor_orang <-
#   rbind(gtf_annot_final_bor_orang_with_hcORF, gtf_annot_final_bor_orang_no_hcORF)
# 
# unique_transcripts_hcORF_sum_orang <- gtf_annot_final_sum_orang %>% 
#   filter(Type == 'transcript') %>% 
#   filter(sseqid != 'NA') %>% 
#   group_by(transcript_id) %>% 
#   filter(length == max(length)) %>% 
#   filter(pident == max(pident)) %>% 
#   filter(qcovs == max(qcovs)) %>% 
#   as.data.frame() 
# 
# gtf_annot_final_sum_orang_with_hcORF <-
#   gtf_annot_final_sum_orang %>% 
#   filter(sseqid %in% unique_transcripts_hcORF_sum_orang$qseqid)
# 
# gtf_annot_final_sum_orang_no_hcORF <-
#   gtf_annot_final_sum_orang %>% 
#   filter(sseqid == 'NA')
# 
# gtf_annot_final_sum_orang <-
#   rbind(gtf_annot_final_sum_orang_with_hcORF, gtf_annot_final_sum_orang_no_hcORF)
# 
# unique_transcripts_hcORF_human <- gtf_annot_final_human %>% 
#   filter(Type == 'transcript') %>% 
#   filter(sseqid != 'NA') %>% 
#   group_by(transcript_id) %>% 
#   filter(length == max(length)) %>% 
#   filter(pident == max(pident)) %>% 
#   filter(qcovs == max(qcovs)) %>% 
#   as.data.frame() 
# 
# gtf_annot_final_human_with_hcORF <-
#   gtf_annot_final_human %>% 
#   filter(sseqid %in% unique_transcripts_hcORF_human$qseqid)
# 
# gtf_annot_final_human_no_hcORF <-
#   gtf_annot_final_human %>% 
#   filter(sseqid == 'NA')
# 
# gtf_annot_final_human <-
#   rbind(gtf_annot_final_human_with_hcORF, gtf_annot_final_human_no_hcORF)

################################################################################

#recode factors
gtf_annot_final_bonobo <-gtf_annot_final_bonobo %>%
  mutate(across(c(untargeted, targeted1, targeted2), ~ ifelse(. == "-", 0, 1)))
gtf_annot_final_chimp <-gtf_annot_final_chimp %>%
  mutate(across(c(untargeted, targeted1, targeted2), ~ ifelse(. == "-", 0, 1)))
gtf_annot_final_gor <-gtf_annot_final_gor %>%
  mutate(across(c(untargeted, targeted1, targeted2), ~ ifelse(. == "-", 0, 1)))
gtf_annot_final_bor_orang <-gtf_annot_final_bor_orang %>%
  mutate(across(c(untargeted, targeted1, targeted2), ~ ifelse(. == "-", 0, 1)))
gtf_annot_final_sum_orang <-gtf_annot_final_sum_orang %>%
  mutate(across(c(untargeted, targeted1, targeted2), ~ ifelse(. == "-", 0, 1)))
gtf_annot_final_human <-gtf_annot_final_human %>%
  mutate(across(c(untargeted, targeted1, targeted2), ~ ifelse(. == "-", 0, 1)))

#subset replicate-supported transcripts
gtf_annot_final_bonobo_rep <- gtf_annot_final_bonobo %>%
  mutate(sum_col = rowSums(across(c(untargeted, targeted1, targeted2)))) %>%
  filter(sum_col > 1)
gtf_annot_final_chimp_rep <- gtf_annot_final_chimp %>%
  mutate(sum_col = rowSums(across(c(untargeted, targeted1, targeted2)))) %>%
  filter(sum_col > 1)
gtf_annot_final_gor_rep <- gtf_annot_final_gor %>%
  mutate(sum_col = rowSums(across(c(untargeted, targeted1, targeted2)))) %>%
  filter(sum_col > 1)
gtf_annot_final_bor_orang_rep <- gtf_annot_final_bor_orang %>%
  mutate(sum_col = rowSums(across(c(untargeted, targeted1, targeted2)))) %>%
  filter(sum_col > 1)
gtf_annot_final_sum_orang_rep <- gtf_annot_final_sum_orang %>%
  mutate(sum_col = rowSums(across(c(untargeted, targeted1, targeted2)))) %>%
  filter(sum_col > 1)
gtf_annot_final_human_rep <- gtf_annot_final_human %>%
  mutate(sum_col = rowSums(across(c(untargeted, targeted1, targeted2)))) %>%
  filter(sum_col > 1)

#merge all gff files into one table
#replicated-supported
gff_annot_final_species_rep <- rbind(gtf_annot_final_bonobo_rep,
                                 gtf_annot_final_chimp_rep, 
                                 gtf_annot_final_gor_rep,
                                 gtf_annot_final_bor_orang_rep, 
                                 gtf_annot_final_sum_orang_rep,
                                 gtf_annot_final_human_rep)

#including non replicate-supported
gff_annot_final_species_all <- rbind(gtf_annot_final_bonobo,
                                     gtf_annot_final_chimp, 
                                     gtf_annot_final_gor,
                                     gtf_annot_final_bor_orang, 
                                     gtf_annot_final_sum_orang,
                                     gtf_annot_final_human)

protein_list <- read.delim("C:/Users/gresh/Downloads/temp/human_apes_Y_NCBI.tsv", header=TRUE)

high_conf_cORF_bonobo <- read.table("C:/Users/gresh/Downloads/temp/bonobo_techrep1_Y_homologs.Ychr.txt", header=FALSE)
high_conf_cORF_chimp <- read.table("C:/Users/gresh/Downloads/temp/chimp_techrep1_Y_homologs.Ychr.txt", header=FALSE)
high_conf_cORF_gor <- read.table("C:/Users/gresh/Downloads/temp/gor_techrep1_Y_homologs.Ychr.txt", header=FALSE)
high_conf_cORF_bor_orang <- read.table("C:/Users/gresh/Downloads/temp/bor_orang_techrep1_Y_homologs.Ychr.txt", header=FALSE)
high_conf_cORF_sum_orang <- read.table("C:/Users/gresh/Downloads/temp/sum_orang_techrep1_Y_homologs.Ychr.txt", header=FALSE)
high_conf_cORF_human <- read.table("C:/Users/gresh/Downloads/temp/human_techrep1_Y_homologs.Ychr.txt", header=FALSE)

#find new gene copies

new_gene_copies <- function(gff_file_rep, species){
  
  new_gene_copies_species <-
    gff_file_rep %>% 
    filter(ref_gene_id == 'NA') %>% 
    filter(lastz_annot != 'NA') %>% 
    dplyr::select(c(gene_id, annot_final)) %>% 
    distinct() %>% 
    group_by(annot_final) %>% 
    summarise(n = n()) %>% 
    mutate(species = species)
  
  return(new_gene_copies_species)
  
}

new_gene_copies_rep_bonobo <- new_gene_copies(gtf_annot_final_bonobo_rep, 'bonobo')
new_gene_copies_rep_chimp <- new_gene_copies(gtf_annot_final_chimp_rep, 'chimp')
new_gene_copies_rep_gor <- new_gene_copies(gtf_annot_final_gor_rep, 'gor')
new_gene_copies_rep_bor_orang <- new_gene_copies(gtf_annot_final_bor_orang_rep, 'bor_orang')
new_gene_copies_rep_sum_orang <- new_gene_copies(gtf_annot_final_sum_orang_rep, 'sum_orang')
new_gene_copies_rep_human <- new_gene_copies(gtf_annot_final_human_rep, 'human')

new_gene_copies_rep_all_species <- rbind(new_gene_copies_rep_bonobo,
                                         new_gene_copies_rep_chimp,
                                         new_gene_copies_rep_gor,
                                         new_gene_copies_rep_bor_orang,
                                         new_gene_copies_rep_sum_orang,
                                         new_gene_copies_rep_human)

new_gene_copies_rep_all_species$type <- 'rep'

new_gene_copies_bonobo <- new_gene_copies(gtf_annot_final_bonobo, 'bonobo')
new_gene_copies_chimp <- new_gene_copies(gtf_annot_final_chimp, 'chimp')
new_gene_copies_gor <- new_gene_copies(gtf_annot_final_gor, 'gor')
new_gene_copies_bor_orang <- new_gene_copies(gtf_annot_final_bor_orang, 'bor_orang')
new_gene_copies_sum_orang <- new_gene_copies(gtf_annot_final_sum_orang, 'sum_orang')
new_gene_copies_human <- new_gene_copies(gtf_annot_final_human, 'human')

################################################################################

potentaily_new_transcripts_bonobo <- 
  gtf_annot_final_bonobo %>% 
  filter(ref_gene_id == 'NA') %>% 
  filter(lastz_annot != 'NA') %>% 
  filter(Type == 'transcript')

potentaily_new_transcripts_bor_orang <- 
  gtf_annot_final_bor_orang %>% 
  filter(ref_gene_id == 'NA') %>% 
  filter(lastz_annot != 'NA') %>% 
  filter(Type == 'transcript')

potentaily_new_transcripts_sum_orang <- 
  gtf_annot_final_sum_orang %>% 
  filter(ref_gene_id == 'NA') %>% 
  filter(lastz_annot != 'NA') %>% 
  filter(Type == 'transcript')

potentaily_new_transcripts_gor <- 
  gtf_annot_final_gor %>% 
  filter(ref_gene_id == 'NA') %>% 
  filter(lastz_annot != 'NA') %>% 
  filter(Type == 'transcript')

potentaily_new_transcripts_chimp <- 
  gtf_annot_final_chimp %>% 
  filter(ref_gene_id == 'NA') %>% 
  filter(lastz_annot != 'NA') %>% 
  filter(Type == 'transcript')

new_gene_copies_all_species <- rbind(new_gene_copies_bonobo,
                                         new_gene_copies_chimp,
                                         new_gene_copies_gor,
                                         new_gene_copies_bor_orang,
                                         new_gene_copies_sum_orang,
                                         new_gene_copies_human)

new_gene_copies_all_species$type <- 'non_rep'

new_gene_copies_summary <- rbind(new_gene_copies_rep_all_species, 
                                 new_gene_copies_all_species)

new_gene_copies_summary$annot_final <- factor(new_gene_copies_summary$annot_final,
                                      levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                                 'RBMY', 'TSPY', 'VCY', 'XKRY'))

new_gene_copies_summary$species <- factor(new_gene_copies_summary$species,
                                  levels = c('chimp','bonobo', 'gor', 
                                             'bor_orang', 'sum_orang', 'human'))

new_gene_copies_summary <-
  new_gene_copies_summary %>% 
  mutate(species = recode(species, "chimp" = "chimpanzee", "bonobo" = "bonobo", 
                          "bor_orang" = "Bornean orangutan", "sum_orang" = "Sumatran orangutan",
                          "human" =  "human", "gor" = "gorilla"))

plot_new_gene_copies <-
  ggplot(new_gene_copies_summary,
         aes(x=annot_final, y=n, fill=species, alpha=type))+
  geom_col(stat = "identity", position = "stack", color='black')+
  facet_wrap(vars(species))+
  xlab('Ampliconic gene family')+
  ylab('New genes copies')+
  theme_pubclean()+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

################################################################################
#new gene copies - human

# library(writexl)
# 
# gtf_annot_final_human %>% 
#   filter(ref_gene_id == 'NA') %>% 
#   filter(lastz_annot != 'NA') %>% 
#   filter(Type == 'transcript') %>% 
#   write_xlsx("C:/Users/gresh/Downloads/temp/new_copies_human.xlsx")

################################################################################
#count transcripts

count_transcripts <- function(gff_file){
  
  sum_transcripts <- 
    gff_file %>% 
    dplyr::select(c('transcript_id', 'annot_final')) %>% 
    distinct() %>% 
    group_by(annot_final) %>% 
    summarise(n = n())
  
  return(sum_transcripts)
  
}

n_transcripts_rep_all_species <-
  gff_annot_final_species_rep %>%
  group_split(species) %>%
  map_dfr(~ count_transcripts(.x) %>% mutate(species = unique(.x$species)), .id = NULL)

n_transcripts_rep_all_species <-
  n_transcripts_rep_all_species %>% 
  mutate(species = recode(species, "chimp" = "chimpanzee", "bonobo" = "bonobo", 
                          "bor_orang" = "Bornean orangutan", "sum_orang" = "Sumatran orangutan",
                          "human" =  "human", "gor" = "gorilla"))

n_transcripts_rep_all_species$annot_final = factor(n_transcripts_rep_all_species$annot_final, 
                           levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                      'RBMY', 'TSPY', 'VCY', 'XKRY'))

n_transcripts_rep_all_species$species = factor(n_transcripts_rep_all_species$species, 
                       levels = c('chimpanzee', 'bonobo', 'gorilla', 
                                  'Bornean orangutan', 'Sumatran orangutan', 'human'))

n_transcripts_rep <-
  ggplot(n_transcripts_rep_all_species,
         aes(x=annot_final, y=n, fill=species))+
  geom_col(color='black',
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label = n), 
            position = position_dodge2(width=.9, preserve='single'), 
            vjust=-0.5, size=3)+
  xlab('Ampliconic gene family')+
  ylab('Transcript count')+
  theme_pubclean()+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

n_transcripts_all_all_species <-
  gff_annot_final_species_all %>%
  group_split(species) %>%
  map_dfr(~ count_transcripts(.x) %>% mutate(species = unique(.x$species)), .id = NULL)

n_transcripts_all_all_species <-
  n_transcripts_all_all_species %>% 
  mutate(species = recode(species, "chimp" = "chimpanzee", "bonobo" = "bonobo", 
                          "bor_orang" = "Bornean orangutan", "sum_orang" = "Sumatran orangutan",
                          "human" =  "human", "gor" = "gorilla"))

n_transcripts_all_all_species$annot_final = factor(n_transcripts_all_all_species$annot_final, 
                                                   levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                                              'RBMY', 'TSPY', 'VCY', 'XKRY'))

n_transcripts_all_all_species$species = factor(n_transcripts_all_all_species$species, 
                                               levels = c('chimpanzee', 'bonobo', 'gorilla', 
                                                          'Bornean orangutan', 'Sumatran orangutan', 'human'))

n_transcripts_all <-
  ggplot(n_transcripts_all_all_species,
         aes(x=annot_final, y=n, fill=species))+
  geom_col(color='black',
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label = n), 
            position = position_dodge2(width=.9, preserve='single'), 
            vjust=-0.5, size=3)+
  xlab('Ampliconic gene family')+
  ylab('Transcript count')+
  theme_pubclean()+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

################################################################################
#transcripts length

gff_annot_final_species_rep$transcript_len <- as.numeric(gff_annot_final_species_rep$transcript_len)

gff_annot_final_species_rep$width_NA <- as.numeric(gff_annot_final_species_rep$width_NA)

gff_annot_final_species_rep$annot_final = factor(gff_annot_final_species_rep$annot_final, 
                                             levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                                        'RBMY', 'TSPY', 'VCY', 'XKRY'))

gff_annot_final_species_rep$species = factor(gff_annot_final_species_rep$species, 
                                         levels = c('chimp', 'bonobo', 
                                                    'gor', 'bor_orang', 'sum_orang', 'human'))

len_transcript_rep <-
  ggpubr::ggboxplot(gff_annot_final_species_rep[which(gff_annot_final_species_rep$Type == 'transcript'),],
                    x= "annot_final", y = "transcript_len",
                    color = "species", size = 0.3,
                    add = "point", add.params = list(size = 1, alpha = 0.25),
                    xlab = 'Ampliconic gene family', ylab = 'Transcript length')+
  scale_color_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

################################################################################

gff_annot_final_species_all$transcript_len <- as.numeric(gff_annot_final_species_all$transcript_len)

gff_annot_final_species_all$annot_final = factor(gff_annot_final_species_all$annot_final, 
                                                 levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                                            'RBMY', 'TSPY', 'VCY', 'XKRY'))

gff_annot_final_species_all$species = factor(gff_annot_final_species_all$species, 
                                             levels = c('chimp', 'bonobo', 
                                                        'gor', 'bor_orang', 'sum_orang', 'human'))

len_transcript_all <-
  ggpubr::ggboxplot(gff_annot_final_species_all,
                    x= "annot_final", y = "transcript_len",
                    color = "species", size = 0.3,
                    add = "point", add.params = list(size = 1, alpha = 0.25),
                    xlab = 'Ampliconic gene family', ylab = 'Transcript length')+
  scale_color_manual(values = c("bonobo" = "#AA4499",
                                "chimp" = "#F0E442",
                                "gor" = "#0072B2",
                                "bor_orang" = "#D55E00",
                                "sum_orang" = "#888888",
                                "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

gff_annot_final_species_rep$transcript_len <- as.numeric(gff_annot_final_species_rep$transcript_len)

gff_annot_final_species_rep$annot_final = factor(gff_annot_final_species_rep$annot_final, 
                                                 levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                                            'RBMY', 'TSPY', 'VCY', 'XKRY'))

gff_annot_final_species_rep$species = factor(gff_annot_final_species_rep$species, 
                                             levels = c('chimp', 'bonobo', 
                                                        'gor', 'bor_orang', 'sum_orang', 'human'))

len_transcript_rep <-
  ggpubr::ggboxplot(gff_annot_final_species_rep,
                    x= "annot_final", y = "transcript_len",
                    color = "species", size = 0.3,
                    add = "point", add.params = list(size = 1, alpha = 0.25),
                    xlab = 'Ampliconic gene family', ylab = 'Transcript length')+
  scale_color_manual(values = c("bonobo" = "#AA4499",
                                "chimp" = "#F0E442",
                                "gor" = "#0072B2",
                                "bor_orang" = "#D55E00",
                                "sum_orang" = "#888888",
                                "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

################################################################################
#n_transcripts_with_cORF 

n_transcripts_with_cORF_all <-
  gff_annot_final_species_all %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  filter(Type == 'transcript') %>% 
  dplyr::select(c(transcript_id, annot_final, species)) %>% 
  distinct() %>% 
  group_by(species, annot_final) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

n_transcripts_with_cORF_all <-
  n_transcripts_with_cORF_all %>% 
  mutate(species = recode(species, "chimp" = "chimpanzee", "bonobo" = "bonobo", 
                          "bor_orang" = "Bornean orangutan", "sum_orang" = "Sumatran orangutan",
                          "human" =  "human", "gor" = "gorilla"))

n_transcripts_with_cORF_all$annot_final = factor(n_transcripts_with_cORF_all$annot_final, 
                                                 levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                                            'RBMY', 'TSPY', 'VCY', 'XKRY'))

n_transcripts_with_cORF_all$species = factor(n_transcripts_with_cORF_all$species, 
                                             levels = c('bonobo', 'chimpanzee', 
                                                        'Bornean orangutan', 'Sumatran orangutan', 
                                                        'gorilla', 'human'))

plot_n_transcripts_with_cORF_all <-
  ggplot(n_transcripts_with_cORF_all,
         aes(x=annot_final, y=n, fill=species))+
  geom_col(color='black',
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label = n), 
            position = position_dodge2(width=.9, preserve='single'), 
            vjust=-0.5, size=3)+
  xlab('Ampliconic gene family')+
  ylab('Number of transcripts with high-confidence cORF')+
  theme_pubclean()+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

#replicate-supported only
n_transcripts_with_cORF_rep <-
  gff_annot_final_species_rep %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  filter(Type == 'transcript') %>% 
  dplyr::select(c(transcript_id, annot_final, species)) %>% 
  distinct() %>% 
  group_by(species, annot_final) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

n_transcripts_with_cORF_rep <-
  n_transcripts_with_cORF_rep %>% 
  mutate(species = recode(species, "chimp" = "chimpanzee", "bonobo" = "bonobo", 
                          "bor_orang" = "Bornean orangutan", "sum_orang" = "Sumatran orangutan",
                          "human" =  "human", "gor" = "gorilla"))

n_transcripts_with_cORF_rep$annot_final = factor(n_transcripts_with_cORF_rep$annot_final, 
                                                 levels = c('BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 
                                                            'RBMY', 'TSPY', 'VCY', 'XKRY'))

n_transcripts_with_cORF_rep$species = factor(n_transcripts_with_cORF_rep$species, 
                                             levels = c('bonobo', 'chimpanzee', 
                                                        'Bornean orangutan', 'Sumatran orangutan', 
                                                        'gorilla', 'human'))

plot_n_transcripts_with_cORF_rep <-
  ggplot(n_transcripts_with_cORF_rep,
         aes(x=annot_final, y=n, fill=species))+
  geom_col(color='black',
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label = n), 
            position = position_dodge2(width=.9, preserve='single'), 
            vjust=-0.5, size=3)+
  xlab('Ampliconic gene family')+
  ylab('Number of replicate-supported transcripts with high-confidence cORF')+
  theme_pubclean()+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

################################################################################

recode_cORF <- function(hig_conf_cORF){
  
  colnames(hig_conf_cORF) <- c('qseqid', 'sseqid', 'qcovs', 'length', 'pident', 
                               'evalue', 'bitscore', 'mismatch', 'gaps', 'qstart', 
                               'qend', 'sstart', 'send', 'qseq', 'sseq', 'qlen')
  
  hig_conf_cORF$transcript_id <- sub("_\\d+$", "", hig_conf_cORF$qseqid)
  
  hig_conf_cORF <-
    hig_conf_cORF %>% 
    filter(sseqid %in% c(protein_list$id))
  
  recoding_df <- data.frame(
    original_values = protein_list$id,
    new_values = protein_list$gene_fam
  )
  
  recode_dict <- recoding_df %>%
    pull(new_values) %>%
    setNames(recoding_df$original_values)
  
  hig_conf_cORF_recoded <- hig_conf_cORF %>%
    mutate(sseqid = recode(sseqid, !!!recode_dict))
  
  return(hig_conf_cORF_recoded)
  
}

high_conf_cORF_bonobo_recoded <- recode_cORF(high_conf_cORF_bonobo)
high_conf_cORF_chimp_recoded <- recode_cORF(high_conf_cORF_chimp)
high_conf_cORF_gor_recoded <- recode_cORF(high_conf_cORF_gor)
high_conf_cORF_bor_orang_recoded <- recode_cORF(high_conf_cORF_bor_orang)
high_conf_cORF_sum_orang_recoded <- recode_cORF(high_conf_cORF_sum_orang)
high_conf_cORF_human_recoded <- recode_cORF(high_conf_cORF_human)

bonobo_transcripts_with_cORF <-
  gff_annot_final_species_all %>% 
  filter(species == 'bonobo') %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_bonobo_recoded <-
  high_conf_cORF_bonobo_recoded %>% 
  filter(transcript_id %in% unique(bonobo_transcripts_with_cORF$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'bonobo')

chimp_transcripts_with_cORF <-
  gff_annot_final_species_all %>% 
  filter(species == 'chimp') %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_chimp_recoded <-
  high_conf_cORF_chimp_recoded %>% 
  filter(transcript_id %in% unique(chimp_transcripts_with_cORF$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'chimp')

gor_transcripts_with_cORF <-
  gff_annot_final_species_all %>% 
  filter(species == 'gor') %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_gor_recoded <-
  high_conf_cORF_gor_recoded %>% 
  filter(transcript_id %in% unique(gor_transcripts_with_cORF$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'gor')

bor_orang_transcripts_with_cORF <-
  gff_annot_final_species_all %>% 
  filter(species == 'bor_orang') %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_bor_orang_recoded <-
  high_conf_cORF_bor_orang_recoded %>% 
  filter(transcript_id %in% unique(bor_orang_transcripts_with_cORF$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'bor_orang')

sum_orang_transcripts_with_cORF <-
  gff_annot_final_species_all %>% 
  filter(species == 'sum_orang') %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_sum_orang_recoded <-
  high_conf_cORF_sum_orang_recoded %>% 
  filter(transcript_id %in% unique(sum_orang_transcripts_with_cORF$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'sum_orang')

human_transcripts_with_cORF <-
  gff_annot_final_species_all %>% 
  filter(species == 'human') %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_human_recoded <-
  high_conf_cORF_human_recoded %>% 
  filter(transcript_id %in% unique(human_transcripts_with_cORF$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'human')

high_conf_cORF_species_recoded_all <- rbind(high_conf_cORF_sum_orang_recoded,
                                       high_conf_cORF_bor_orang_recoded,
                                       high_conf_cORF_bonobo_recoded,
                                       high_conf_cORF_chimp_recoded,
                                       high_conf_cORF_gor_recoded,
                                       high_conf_cORF_human_recoded)

bonobo_transcripts_with_cORF_rep <-
  gtf_annot_final_bonobo_rep %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_bonobo_recoded_rep <-
  high_conf_cORF_bonobo_recoded %>% 
  filter(transcript_id %in% unique(bonobo_transcripts_with_cORF_rep$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'bonobo')

chimp_transcripts_with_cORF_rep <-
  gtf_annot_final_chimp_rep %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_chimp_recoded_rep <-
  high_conf_cORF_chimp_recoded %>% 
  filter(transcript_id %in% unique(chimp_transcripts_with_cORF_rep$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'chimp')

gor_transcripts_with_cORF_rep <-
  gtf_annot_final_gor_rep %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_gor_recoded_rep <-
  high_conf_cORF_gor_recoded %>% 
  filter(transcript_id %in% unique(gor_transcripts_with_cORF_rep$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'gor')

bor_orang_transcripts_with_cORF_rep <-
  gtf_annot_final_bor_orang_rep %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_bor_orang_recoded_rep <-
  high_conf_cORF_bor_orang_recoded %>% 
  filter(transcript_id %in% unique(bor_orang_transcripts_with_cORF_rep$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'bor_orang')

sum_orang_transcripts_with_cORF_rep <-
  gtf_annot_final_sum_orang_rep %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_sum_orang_recoded_rep <-
  high_conf_cORF_sum_orang_recoded %>% 
  filter(transcript_id %in% unique(sum_orang_transcripts_with_cORF_rep$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'sum_orang')

human_transcripts_with_cORF_rep <-
  gtf_annot_final_human_rep %>% 
  filter(sseqid != 'NA') %>% 
  filter(sseqid == annot_final) %>% 
  select(transcript_id) %>% 
  distinct()

high_conf_cORF_human_recoded_rep <-
  high_conf_cORF_human_recoded %>% 
  filter(transcript_id %in% unique(human_transcripts_with_cORF_rep$transcript_id)) %>% 
  distinct() %>% 
  mutate(species = 'human')

high_conf_cORF_species_recoded_rep <- rbind(high_conf_cORF_sum_orang_recoded_rep,
                                           high_conf_cORF_bor_orang_recoded_rep,
                                           high_conf_cORF_bonobo_recoded_rep,
                                           high_conf_cORF_chimp_recoded_rep,
                                           high_conf_cORF_gor_recoded_rep,
                                           high_conf_cORF_human_recoded_rep)

################################################################################
#total number of high-confidence cORFs identified in transcripts

high_conf_cORF_distinct_all <-
  high_conf_cORF_species_recoded_all %>% 
  filter(pident >= 80) %>% 
  filter(qcovs >= 75) %>% 
  group_by(species, transcript_id) %>% 
  select(species, sseqid, transcript_id) %>% 
  distinct() %>% 
  group_by(species, sseqid) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

high_conf_cORF_distinct_all$sseqid <- factor(high_conf_cORF_distinct_all$sseqid,
                                                  levels = c('BPY2', 'CDY', 'DAZ',
                                                             'HSFY', 'PRY', 'RBMY',
                                                             'TSPY', 'VCY', 'XKRY'))

high_conf_cORF_distinct_all$species <- factor(high_conf_cORF_distinct_all$species,
                                                   levels = c('bonobo', 
                                                              'chimp',
                                                              'bor_orang',
                                                              'sum_orang',
                                                              'gor',
                                                              'human'))

transcripts_with_high_conf_cORF_all <-
  ggplot(high_conf_cORF_distinct_all,
         aes(x=sseqid, y=n, fill = species))+
  geom_col(color = 'black', 
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label=n), vjust=-0.5,
            position = position_dodge2(width=.9, preserve='single'), size=3)+
  theme_pubclean()+
  xlab('Ampliconic gene family')+
  ylab('Number of transcripts with high-confidence cORFs')+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

#replicate-supported
high_conf_cORF_distinct_rep <-
  high_conf_cORF_species_recoded_rep %>% 
  filter(pident >= 80) %>% 
  filter(qcovs >= 75) %>% 
  select(species, sseqid, transcript_id) %>% 
  distinct() %>% 
  group_by(species, sseqid) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

high_conf_cORF_distinct_rep$sseqid <- factor(high_conf_cORF_distinct_rep$sseqid,
                                             levels = c('BPY2', 'CDY', 'DAZ',
                                                        'HSFY', 'PRY', 'RBMY',
                                                        'TSPY', 'VCY', 'XKRY'))

high_conf_cORF_distinct_rep$species <- factor(high_conf_cORF_distinct_rep$species,
                                              levels = c('bonobo', 
                                                         'chimp',
                                                         'bor_orang',
                                                         'sum_orang',
                                                         'gor',
                                                         'human'))

transcripts_with_high_conf_cORF_rep <-
  ggplot(high_conf_cORF_distinct_rep,
         aes(x=sseqid, y=n, fill = species))+
  geom_col(color = 'black', 
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label=n), vjust=-0.5,
            position = position_dodge2(width=.9, preserve='single'), size=3)+
  theme_pubclean()+
  xlab('Ampliconic gene family')+
  ylab('Number of transcripts with high-confidence cORFs')+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

################################################################################

#total number of high-confidence cORF 
high_conf_cORF_all <-
  high_conf_cORF_species_recoded_all %>% 
  filter(pident >= 80) %>% 
  filter(qcovs >= 75) %>% 
  select(species, sseqid, qseqid) %>% 
  distinct() %>% 
  group_by(species, sseqid) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

high_conf_cORF_all$sseqid <- factor(high_conf_cORF_all$sseqid,
                                       levels = c('BPY2', 'CDY', 'DAZ',
                                                  'HSFY', 'PRY', 'RBMY',
                                                  'TSPY', 'VCY', 'XKRY'))

high_conf_cORF_all$species <- factor(high_conf_cORF_all$species,
                                     levels = c('bonobo', 
                                                'chimp',
                                                'bor_orang',
                                                'sum_orang',
                                                'gor',
                                                'human'))

plot_high_conf_cORF_all <-
  ggplot(high_conf_cORF_all, 
         aes(x=sseqid, y=n, fill = species))+
  geom_col(color = 'black',
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label=n), 
            vjust=-0.5,
            stat = "identity", position = position_dodge2(width=.9, preserve='single'),
            size = 3)+
  theme_pubclean()+
  xlab('Gene family')+
  ylab('Total number of cORFs homologs')+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

#replicate-supported

high_conf_cORF_rep <-
  high_conf_cORF_species_recoded_rep %>% 
  filter(pident >= 80) %>% 
  filter(qcovs >= 75) %>% 
  select(species, sseqid, qseqid) %>% 
  distinct() %>% 
  group_by(species, sseqid) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

high_conf_cORF_rep$sseqid <- factor(high_conf_cORF_rep$sseqid,
                                       levels = c('BPY2', 'CDY', 'DAZ',
                                                  'HSFY', 'PRY', 'RBMY',
                                                  'TSPY', 'VCY', 'XKRY'))

high_conf_cORF_rep$species <- factor(high_conf_cORF_rep$species,
                                        levels = c('bonobo', 
                                                   'chimp',
                                                   'bor_orang',
                                                   'sum_orang',
                                                   'gor', 
                                                   'human'))

plot_high_conf_cORF_rep <-
  ggplot(high_conf_cORF_rep, 
         aes(x=sseqid, y=n, fill = species))+
  geom_col(color = 'black',
           position = position_dodge2(width=.9, preserve='single'))+
  geom_text(aes(label=n), 
            vjust=-0.5,
            stat = "identity", position = position_dodge2(width=.9, preserve='single'),
            size = 3)+
  theme_pubclean()+
  xlab('Gene family')+
  ylab('Number of ORFs homologs')+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))


################################################################################
#length of high-confidence cORF

high_conf_cORF_all_distinct <-
  high_conf_cORF_species_recoded_all %>% 
  dplyr::select(c(qseqid, sseqid, species, qlen, 
                  length, pident, qcovs, transcript_id)) %>% 
  group_by(species, sseqid, transcript_id) %>% 
  filter(pident >= 80) %>% 
  filter(qcovs >= 75) %>% 
  filter(length == max(length)) %>% 
  filter(pident == max(pident)) %>% 
  filter(qcovs == max(qcovs)) %>% 
  ungroup() %>% 
  distinct() %>% 
  data.frame()

high_conf_cORF_all_distinct$sseqid <- factor(high_conf_cORF_all_distinct$sseqid,
                                            levels = c('BPY2', 'CDY', 'DAZ',
                                                       'HSFY', 'PRY', 'RBMY',
                                                       'TSPY', 'VCY', 'XKRY'))

high_conf_cORF_all_distinct$species <- factor(high_conf_cORF_all_distinct$species,
                                             levels = c('bonobo', 
                                                        'chimp',
                                                        'bor_orang',
                                                        'sum_orang',
                                                        'gor',
                                                        'human'))

len_high_conf_cORF_all <-
  ggpubr::ggboxplot(high_conf_cORF_all_distinct,
                    "sseqid", "qlen",
                    color = "species",
                    add = "point", add.params = list(size = 1, alpha = 0.25),
                    xlab = 'Ampliconic gene family', ylab = 'High-confidence cORF length (aa)')+
  scale_color_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

high_conf_cORF_rep_distinct <-
  high_conf_cORF_species_recoded_rep %>% 
  dplyr::select(c(qseqid, sseqid, species, qlen, 
                  length, pident, qcovs, transcript_id)) %>% 
  group_by(species, sseqid, transcript_id) %>% 
  filter(pident >= 80) %>% 
  filter(qcovs >= 75) %>% 
  filter(length == max(length)) %>% 
  filter(pident == max(pident)) %>% 
  filter(qcovs == max(qcovs)) %>% 
  ungroup() %>% 
  distinct() %>% 
  data.frame()

high_conf_cORF_rep_distinct$sseqid <- factor(high_conf_cORF_rep_distinct$sseqid,
                                            levels = c('BPY2', 'CDY', 'DAZ',
                                                       'HSFY', 'PRY', 'RBMY',
                                                       'TSPY', 'VCY', 'XKRY'))

high_conf_cORF_rep_distinct$species <- factor(high_conf_cORF_rep_distinct$species,
                                             levels = c('bonobo', 
                                                        'chimp',
                                                        'bor_orang',
                                                        'sum_orang',
                                                        'gor',
                                                        'human'))

len_high_conf_cORF_rep <-
  ggpubr::ggboxplot(high_conf_cORF_rep_distinct,
                    "sseqid", "qlen",
                    color = "species",
                    add = "point", add.params = list(size = 1, alpha = 0.25),
                    xlab = 'Ampliconic gene family', ylab = 'High-confidence replicate-supported cORF length (aa)')+
  scale_color_manual(values = c("bonobo" = "#AA4499",
                                "chimp" = "#F0E442",
                                "gor" = "#0072B2",
                                "bor_orang" = "#D55E00",
                                "sum_orang" = "#888888",
                                "human" = "#009E73"))+
  theme(legend.position = "none",                     
        axis.text.x = element_text(face = "bold.italic"))

################################################################################

BPY2_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'BPY2',],
            x="length", y="count", 
            fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                                "chimp" = "#F0E442",
                                "gor" = "#0072B2",
                                "bor_orang" = "#D55E00",
                                "sum_orang" = "#888888",
                                "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="BPY2")+
  ylim(c(0, 50))+
  xlim(c(0, 120))+
  geom_vline(xintercept = 106, linetype = "dashed", size = 0.8)+
  theme(legend.position = "none",
        aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

CDY_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'CDY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="CDY")+
  ylim(c(0, 50))+
  xlim(c(0, 650))+
  geom_vline(xintercept = 470, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 540, linetype = "dashed", size = 0.8)+
  theme(legend.position = "none",
        aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

DAZ_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'DAZ',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="DAZ")+
  geom_vline(xintercept = 390, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 744, linetype = "dashed", size = 0.8)+
  xlim(c(0, 850))+
  ylim(c(0, 50))+
  theme(legend.position = "none",
        aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

HSFY_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'HSFY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="HSFY")+
  geom_vline(xintercept = 203, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 323, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 401, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 450))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

RBMY_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'RBMY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='RBMY')+
  geom_vline(xintercept = 356, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 422, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 496, linetype = "dashed", size = 0.8)+
  ylim(c(0, 125))+
  xlim(c(0, 550))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

PRY_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'PRY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  geom_vline(xintercept = 147, linetype = "dashed", size = 0.8)+
  labs(x=NULL, y=NULL, title="PRY")+
  ylim(c(0, 50))+
  xlim(c(0, 160))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

TSPY_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'TSPY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='TSPY')+
  geom_vline(xintercept = 288, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 308, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 350))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

VCY_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'VCY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='VCY')+
  geom_vline(xintercept = 125, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 150))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

XKRY_plot <-
  gghistogram(high_conf_cORF_all_distinct[high_conf_cORF_all_distinct$sseqid == 'XKRY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='XKRY')+
  geom_vline(xintercept = 117, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 130))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

ggarrange(BPY2_plot, CDY_plot, DAZ_plot,
          HSFY_plot, RBMY_plot, PRY_plot,
          TSPY_plot, VCY_plot, XKRY_plot)

################################################################################

BPY2_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'BPY2',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="BPY2")+
  ylim(c(0, 50))+
  xlim(c(0, 150))+
  geom_vline(xintercept = 106, linetype = "dashed", size = 0.8)+
  theme(legend.position = "none",
        aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

CDY_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'CDY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="CDY")+
  ylim(c(0, 50))+
  xlim(c(0, 650))+
  geom_vline(xintercept = 470, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 540, linetype = "dashed", size = 0.8)+
  theme(legend.position = "none",
        aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

DAZ_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'DAZ',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="DAZ")+
  geom_vline(xintercept = 390, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 744, linetype = "dashed", size = 0.8)+
  theme(legend.position = "none",
        aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))+
  xlim(c(0, 800))+
  ylim(c(0, 50))

HSFY_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'HSFY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title="HSFY")+
  geom_vline(xintercept = 203, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 323, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 401, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 450))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

RBMY_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'RBMY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='RBMY')+
  geom_vline(xintercept = 356, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 422, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 496, linetype = "dashed", size = 0.8)+
  ylim(c(0, 125))+
  xlim(c(0, 550))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

PRY_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'PRY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  geom_vline(xintercept = 147, linetype = "dashed", size = 0.8)+
  labs(x=NULL, y=NULL, title="PRY")+
  ylim(c(0, 50))+
  xlim(c(0, 160))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

TSPY_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'TSPY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='TSPY')+
  geom_vline(xintercept = 288, linetype = "dashed", size = 0.8)+
  geom_vline(xintercept = 308, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 400))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

VCY_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'VCY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='VCY')+
  geom_vline(xintercept = 125, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 150))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

XKRY_plot <-
  gghistogram(high_conf_cORF_rep_distinct[high_conf_cORF_rep_distinct$sseqid == 'XKRY',],
              x="length", y="count", 
              fill="species", alpha=1, position = "stack")+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  labs(x=NULL, y=NULL, title='XKRY')+
  geom_vline(xintercept = 117, linetype = "dashed", size = 0.8)+
  ylim(c(0, 50))+
  xlim(c(0, 150))+
  theme(legend.position = "none", aspect.ratio = .5,
        plot.title = element_text(hjust = 0, face = "bold.italic"))

ggarrange(BPY2_plot, CDY_plot, DAZ_plot,
          HSFY_plot, RBMY_plot,
          TSPY_plot, VCY_plot)

###############################################################################

high_conf_cORF_all_distinct$qcovs <- 
  as.numeric(high_conf_cORF_all_distinct$qcovs)
high_conf_cORF_all_distinct$pident <- 
  as.numeric(high_conf_cORF_all_distinct$pident)

coverage_vs_seqid_high_conf_cORF_all <-
  ggplot(high_conf_cORF_all_distinct)+
  geom_point(aes(x=qcovs, y=pident, color=species))+
  facet_wrap(vars(sseqid))+
  theme_pubclean()+
  xlab('Coverage of predicted YAGs proteins to NCBI YAGs proteins')+
  ylab('Sequence identity between predicted and NCBI YAGs proteins')+
  scale_color_manual(values = c("bonobo" = "#AA4499",
                                "chimp" = "#F0E442",
                                "gor" = "#0072B2",
                                "bor_orang" = "#D55E00",
                                "sum_orang" = "#888888",
                                "human" = "#009E73"))+
  theme(strip.text = element_text(face = "bold.italic"),
        legend.position = "none")

coverage_vs_seqid_high_conf_cORF_rep <-
  ggplot(high_conf_cORF_rep_distinct)+
  geom_point(aes(x=qcovs, y=pident, color=species))+
  facet_wrap(vars(sseqid))+
  theme_pubclean()+
  xlab('Coverage of predicted YAGs proteins to NCBI YAGs proteins')+
  ylab('Sequence identity between predicted and NCBI YAGs proteins')+
  scale_color_manual(values = c("bonobo" = "#AA4499",
                                "chimp" = "#F0E442",
                                "gor" = "#0072B2",
                                "bor_orang" = "#D55E00",
                                "sum_orang" = "#888888",
                                "human" = "#009E73"))+
  theme(strip.text = element_text(face = "bold.italic"),
        legend.position = "none")

number_of_homologs_vs_algn_length_all <-
  ggplot(high_conf_cORF_all_distinct) +
  geom_histogram(aes(x=length, fill = species))+
  facet_wrap(vars(sseqid))+
  theme_pubclean()+
  xlab('Alignment length (aa) of predicted YAGs proteins homologus to NCBI YAGs proteins')+
  ylab('Number of predicted YAGs proteins homologus to NCBI YAGs proteins')+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(strip.text = element_text(face = "bold.italic"),
        legend.position = "none")

number_of_homologs_vs_algn_length_all_rep <-
  ggplot(high_conf_cORF_rep_distinct) +
  geom_histogram(aes(x=length, fill = species))+
  facet_wrap(vars(sseqid))+
  theme_pubclean()+
  xlab('Alignment length (aa) of predicted YAGs proteins homologus to NCBI YAGs proteins')+
  ylab('Number of predicted YAGs proteins homologus to NCBI YAGs proteins')+
  scale_fill_manual(values = c("bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888",
                               "human" = "#009E73"))+
  theme(strip.text = element_text(face = "bold.italic"),
        legend.position = "none")

################################################################################

#count categories - all transcripts

count_cat <- function(gff_file){
  
  count_transcripts_per_gene_family <-
    gff_file %>% 
    filter(Type == 'transcript') %>% 
    group_by(species, annot_final, gene_id) %>% 
    summarise(n = n()) %>% 
    data.frame()
  
  cat_1 <-
    count_transcripts_per_gene_family %>%
    filter(n == 1) %>%
    group_by(species, annot_final) %>%
    summarise(n = n()) %>%
    mutate(cat = '1')
  
  cat_2 <-
    count_transcripts_per_gene_family %>%
    filter(n > 1) %>%
    filter(n < 4) %>%
    group_by(species, annot_final) %>%
    summarise(n = n()) %>%
    mutate(cat = '2-3')
  
  cat_3 <-
    count_transcripts_per_gene_family %>%
    filter(n > 3) %>%
    filter(n < 6) %>%
    group_by(species, annot_final) %>%
    summarise(n = n()) %>%
    mutate(cat = '4-5')
  
  cat_4 <-
    count_transcripts_per_gene_family %>%
    filter(n > 5) %>%
    filter(n < 8) %>%
    group_by(species, annot_final) %>%
    summarise(n = n()) %>%
    mutate(cat = '6-7')
  
  cat_5 <-
    count_transcripts_per_gene_family %>%
    filter(n > 7) %>%
    filter(n < 10) %>%
    group_by(species, annot_final) %>%
    summarise(n = n()) %>%
    mutate(cat = '8-9')
  
  cat_6 <-
    count_transcripts_per_gene_family %>%
    filter(n > 9) %>%
    group_by(species, annot_final) %>%
    summarise(n = n()) %>%
    mutate(cat = '>10')
  
  cat_all <- rbind(cat_1, cat_2, cat_3, cat_4, cat_5, cat_6)
  cat_all <- data.frame(cat_all)
  
  cat_all$annot_final <- factor(cat_all$annot_final,
                                levels = c('BPY2', 'CDY', 'DAZ',
                                           'HSFY', 'PRY', 'RBMY',
                                           'TSPY', 'VCY', 'XKRY'))
  
  cat_all$species <- factor(cat_all$species,
                            levels = c('bonobo', 
                                       'chimp',
                                       'bor_orang',
                                       'sum_orang',
                                       'gor',
                                       'human'))
  
  cat_all <-
    cat_all %>% 
    mutate(species = recode(species, "chimp" = "chimpanzee", "bonobo" = "bonobo", 
                            "bor_orang" = "Bornean orangutan", "sum_orang" = "Sumatran orangutan",
                            "human" =  "human", "gor" = "gorilla"))
  
  return(cat_all)
  
}

categories_gff_annot_final_species_all <- count_cat(gff_annot_final_species_all)

n_isoforms_per_gene_copy_all <-
  ggplot(categories_gff_annot_final_species_all)+
  geom_col(aes(x=cat, y=n, fill=species),
           position = position_dodge2(.9,  preserve = "single"),
           color='black', width=.8)+
  geom_text(aes(x=cat, y=n, label=n), vjust=-0.5,
            position = position_dodge2(.9,  preserve = "single"), size=3)+
  facet_grid(rows=vars(annot_final),
             cols=vars(species))+
  scale_x_discrete(limits=c("1","2-3","4-5","6-7","8-9",">10"))+
  xlab("Number of Isoforms per Gene Copy")+
  ylab("Number of Gene Copies")+
  ylim(0, 50)+
  theme_pubclean()+
  theme(panel.border=element_rect(colour="black", fill=NA, size=1),
        strip.text.y = element_text(face = "bold.italic"),
        legend.position = "none")+
  scale_fill_manual(values = c("human" = "#009E73",
                               "bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888"))

categories_gff_annot_final_species_rep <- count_cat(gff_annot_final_species_rep)

n_isoforms_per_gene_copy_rep <-
  ggplot(categories_gff_annot_final_species_rep)+
  geom_col(aes(x=cat, y=n, fill=species),
           position = position_dodge2(.9,  preserve = "single"),
           color='black', width=.8)+
  geom_text(aes(x=cat, y=n, label=n), vjust=-0.5,
            position = position_dodge2(.9,  preserve = "single"), size=3)+
  facet_grid(rows=vars(annot_final),
             cols=vars(species))+
  scale_x_discrete(limits=c("1","2-3","4-5","6-7","8-9",">10"))+
  xlab("Number of Isoforms per Gene Copy, Replicate-supported")+
  ylab("Number of Gene Copies")+
  ylim(0, 50)+
  theme_pubclean()+
  theme(panel.border=element_rect(colour="black", fill=NA, size=1),
        strip.text.y = element_text(face = "bold.italic"),
        legend.position = "none")+
  scale_fill_manual(values = c("human" = "#009E73",
                               "bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888"))

gff_annot_final_species_all_hc_ORF <-
  gff_annot_final_species_all %>% 
  filter(annot_final == sseqid) 

categories_gff_annot_final_species_all_hc_ORF <- count_cat(gff_annot_final_species_all_hc_ORF)

n_isoforms_per_gene_copy_hc_ORF_all <-
  ggplot(categories_gff_annot_final_species_all_hc_ORF)+
  geom_col(aes(x=cat, y=n, fill=species),
           position = position_dodge2(.9,  preserve = "single"),
           color='black', width=.8)+
  geom_text(aes(x=cat, y=n, label=n), vjust=-0.5,
            position = position_dodge2(.9,  preserve = "single"), size=3)+
  facet_grid(rows=vars(annot_final),
             cols=vars(species))+
  scale_x_discrete(limits=c("1","2-3","4-5","6-7","8-9",">10"))+
  xlab("Number of Isoforms per Gene Copy, All transcripts with high-confidence cORFs")+
  ylab("Number of Gene Copies")+
  ylim(0, 50)+
  theme_pubclean()+
  theme(panel.border=element_rect(colour="black", fill=NA, size=1),
        strip.text.y = element_text(face = "bold.italic"),
        legend.position = "none")+
  scale_fill_manual(values = c("human" = "#009E73",
                               "bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888"))

gff_annot_final_species_rep_hc_ORF <-
  gff_annot_final_species_rep %>% 
  filter(annot_final == sseqid) 

categories_gff_annot_final_species_rep_hc_ORF <- count_cat(gff_annot_final_species_rep_hc_ORF)

n_isoforms_per_gene_copy_hc_ORF_rep <-
  ggplot(categories_gff_annot_final_species_rep_hc_ORF)+
  geom_col(aes(x=cat, y=n, fill=species),
           position = position_dodge2(.9,  preserve = "single"),
           color='black', width=.8)+
  geom_text(aes(x=cat, y=n, label=n), vjust=-0.5,
            position = position_dodge2(.9,  preserve = "single"), size=3)+
  facet_grid(rows=vars(annot_final),
             cols=vars(species))+
  scale_x_discrete(limits=c("1","2-3","4-5","6-7","8-9",">10"))+
  xlab("Number of Isoforms per Gene Copy, Replicate-supported transcripts with high-confidence cORFs")+
  ylab("Number of Gene Copies")+
  ylim(0, 50)+
  theme_pubclean()+
  theme(panel.border=element_rect(colour="black", fill=NA, size=1),
        strip.text.y = element_text(face = "bold.italic"),
        legend.position = "none")+
  scale_fill_manual(values = c("human" = "#009E73",
                               "bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888"))

################################################################################

#number of identified gene copies and pseudo genes 

#number of gene copies
number_of_gene_copies <-
  gff_annot_final_species_all %>% 
  select(species, annot_final, gene_id, sseqid) %>% 
  filter(annot_final == sseqid) %>% 
  distinct() %>% 
  group_by(species, annot_final) %>% 
  summarise(n_gene_copies = n()) %>% 
  data.frame()

#number of pseudo genes 
number_of_pseudogenes <-
  gff_annot_final_species_all %>% 
  select(species, annot_final, gene_id, sseqid) %>% 
  filter(sseqid == 'NA') %>% 
  distinct() %>% 
  group_by(species, annot_final) %>% 
  summarise(n_pseudogenes = n()) %>% 
  data.frame()

number_genes_pseudogenes <- 
  merge(number_of_gene_copies, number_of_pseudogenes, all=TRUE) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  mutate(n_pseudogenes_negative = -n_pseudogenes) %>% 
  select(c(species, annot_final, n_gene_copies, n_pseudogenes_negative)) %>% 
  pivot_longer(cols = n_gene_copies:n_pseudogenes_negative) %>% 
  data.frame()

number_genes_pseudogenes$annot_final <- factor(number_genes_pseudogenes$annot_final,
                              levels = c('BPY2', 'CDY', 'DAZ',
                                         'HSFY', 'PRY', 'RBMY',
                                         'TSPY', 'VCY', 'XKRY'))

number_genes_pseudogenes$species <- factor(number_genes_pseudogenes$species,
                          levels = c('bonobo', 
                                     'chimp',
                                     'bor_orang',
                                     'sum_orang',
                                     'gor',
                                     'human'))

number_genes_pseudogenes <-
  number_genes_pseudogenes %>% 
  mutate(species = recode(species, "chimp" = "chimpanzee", "bonobo" = "bonobo", 
                          "bor_orang" = "Bornean orangutan", "sum_orang" = "Sumatran orangutan",
                          "human" =  "human", "gor" = "gorilla"))

genes_pseudogenes <-
  ggplot(number_genes_pseudogenes)+
  geom_col(aes(x=value, y=species, alpha=name, fill=species), color='black')+
  facet_wrap(vars(annot_final))+
  theme_pubclean()+
  coord_flip()+
  scale_alpha_discrete(range = c(1, 0.35))+
  theme(strip.text = element_text(face = "bold.italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none")+
  scale_fill_manual(values = c("human" = "#009E73",
                               "bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888"))+
  ylab('Species')+
  xlab('Number of copies')+
  xlim(c(-45, 45))


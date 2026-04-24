library(tidyr)
library(ggplot2)
library(ggpubr)
library(microseq)

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

#convert df to gff format (column 'Attributes')
df_to_gff <- function(df){
  
  if (nrow(df) > 0){
    
    attr_only_df <-
      df %>%
      dplyr::select(-c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase"))
    
    no_attr_df <-
      df %>%
      dplyr::select(c("Seqid", "Source", "Type", "Start", "End", "Score", "Strand", "Phase"))
    
    attr_only_bq <- data.frame(lapply(attr_only_df,
                                      function(x) paste('"', paste(x, '"', sep=''), sep='')))
    attr_only_merged <- apply(attr_only_bq, 1, function(row) paste(names(row), row, sep = " ", collapse = ";"))
    attr_only_merged <- data.frame(Attributes = attr_only_merged)
    
    gff_attr_from_ref <- cbind(no_attr_df, attr_only_merged)
    
    return(gff_attr_from_ref)
    
  } else {
    
    return(0)
    
  }
  
}

################################################################################

species_list <- c('bonobo', 'bor_orang', 'sum_orang', 'gor', 'chimp', 'human')
for (species in species_list){
  
  ################################################################################
  #read gff file
  
  read_species_ref <- readGFF(paste("/galaxy/home/ajg7274/isoseq_apes/mix_union_targeted_untargeted/compare_with_original/reference/", species, "_techrep1.stringtie.annotation.hc_cORF.Ychr.gff", sep=''))
  
  read_species_ref <-
    read_species_ref %>% 
    filter(Start != 'NA')
  
  read_species_ref_attr <- gff_attr(read_species_ref)
  
  ################################################################################
  #read tmap file untargeted
  
  ref_untargeted <-
    read.table(paste("/galaxy/home/ajg7274/isoseq_apes/mix_union_targeted_untargeted/compare_with_original/untargeted/", species, "_untargeted.", species, "_techrep1.stringtie.Ychr.gtf.tmap", sep=''), 
               sep = "\t", header = TRUE)
  
  ref_untargeted_eq <-
    ref_untargeted %>% 
    filter(class_code == '=')
  
  untargeted_list <- ref_untargeted_eq$ref_id
  
  ################################################################################
  #read tmap file targeted 1
  
  ref_targeted_1 <-
    read.table(paste("/galaxy/home/ajg7274/isoseq_apes/mix_union_targeted_untargeted/compare_with_original/targeted/", species, "_targeted1.", species, "_techrep1.stringtie.Ychr.gtf.tmap",
                     sep=''), sep = "\t", header = TRUE)
  
  ref_targeted_1_eq <-
    ref_targeted_1 %>% 
    filter(class_code == '=')
  
  targeted_1_list <- ref_targeted_1_eq$ref_id
  
  ################################################################################
  #read tmap file targeted 2
  
  ref_targeted_2 <-
    read.table(paste("/galaxy/home/ajg7274/isoseq_apes/mix_union_targeted_untargeted/compare_with_original/targeted/", species, "_targeted2.", species, "_techrep2.stringtie.Ychr.gtf.tmap", 
                     sep=''), sep = "\t", header = TRUE)
  
  ref_targeted_2_eq <-
    ref_targeted_2 %>% 
    filter(class_code == '=')
  
  targeted_2_list <- ref_targeted_2_eq$ref_id
  
  ################################################################################
  
  #add three columns "targeted1", "targeted2", "untargeted" reflecting the origin of the transcripts in the merged file
  read_species_ref_attr$untargeted <- ifelse(read_species_ref_attr$transcript_id %in% untargeted_list, "=", "-")
  read_species_ref_attr$targeted1 <- ifelse(read_species_ref_attr$transcript_id %in% targeted_1_list, "=", "-")
  read_species_ref_attr$targeted2 <- ifelse(read_species_ref_attr$transcript_id %in% targeted_2_list, "=", "-")
  
  ################################################################################
  #write to gff file
  
  read_species_ref_attr <- read_species_ref_attr %>% filter(annot_final != 'NA')

  ouput_name_gff <- paste("/galaxy/home/ajg7274/isoseq_apes/mix_union_targeted_untargeted/compare_with_original/source/", species, "_techrep1.stringtie.annotation.hc_cORF.Ychr.source.gff", sep='')
  
  writeGFF(df_to_gff(read_species_ref_attr), ouput_name_gff)

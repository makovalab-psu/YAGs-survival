library(stringr)
library(dplyr)

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

recode_cORF <- function(high_conf_cORF_input){
  
  high_conf_cORF_input$transcript_id <- sub("_\\d+$", "", high_conf_cORF_input$qseqid)
  
  high_conf_cORF_input <-
    high_conf_cORF_input %>% 
    filter(sseqid %in% c(protein_list$id))
  
  recoding_df <- data.frame(
    original_values = protein_list$id,
    new_values = protein_list$gene_fam
  )
  
  recode_dict <- recoding_df %>%
    pull(new_values) %>%
    setNames(recoding_df$original_values)
  
  high_conf_cORF_input$sseqid_original <- high_conf_cORF_input$sseqid
  
  high_conf_cORF_input_recoded <- high_conf_cORF_input %>%
    mutate(sseqid = recode(sseqid, !!!recode_dict))
  
  return(high_conf_cORF_input_recoded)
  
}

################################################################################

# protein_lengths <-
#    read.table("C:/Users/gresh/Downloads/temp/protein_lengths.tsv", header = FALSE,
#            sep = "\t", stringsAsFactors = FALSE)
# 
# protein_list <- read.delim("C:/Users/gresh/Downloads/temp/human_apes_Y_NCBI.tsv", header=TRUE)
# 
# gtf_annot <- microseq::readGFF("C:/Users/gresh/Downloads/temp/bor_orang_techrep1.stringtie.annotation.final.Ychr.gff")
# 
# high_conf_cORF <- read.table("C:/Users/gresh/Downloads/temp/bor_orang_techrep1_Y_homologs.Ychr.txt", header=FALSE)

protein_lengths <-
  read.table(snakemake@input[[1]], header = FALSE,
             sep = "\t",
             stringsAsFactors = FALSE)
protein_list <- read.delim(snakemake@input[[2]], header=TRUE)
gtf_annot <- microseq::readGFF(snakemake@input[[3]])
high_conf_cORF <- read.table(snakemake@input[[4]], header=FALSE)

################################################################################

protein_lengths$V1 <- str_extract(protein_lengths$V1, "^[^ ]+")

colnames(protein_lengths) <- c('sseqid_original', 'slen')

gtf_annot <- gff_attr(gtf_annot)

################################################################################

colnames(high_conf_cORF) <- c("qseqid", "sseqid", "qcovs", "length", 
                                        "pident", "evalue", "bitscore", "mismatch", 
                                        "gaps", "qstart", "qend", "sstart", "send", 
                                        "qseq", "sseq", "qlen")

high_conf_cORF_recoded <- recode_cORF(high_conf_cORF)

high_conf_cORF_recoded <- left_join(high_conf_cORF_recoded, protein_lengths)

high_conf_cORF_recoded <- 
  high_conf_cORF_recoded %>% 
  mutate(scovs = (length - gaps)/slen*100) %>% 
  dplyr::select(qseqid, sseqid, qcovs, scovs, length, pident, qlen, transcript_id)

high_conf_cORF_recoded_unique <-
  high_conf_cORF_recoded %>% 
  filter(scovs >= 50) %>% 
  filter(qcovs >= 80) %>% 
  filter(pident >= 80) %>% 
  group_by(transcript_id) %>% 
  filter(scovs == max(scovs)) %>% 
  filter(qcovs == max(qcovs)) %>% 
  filter(pident == max(pident)) %>% 
  ungroup() %>% 
  select(-c(qseqid)) %>% 
  distinct() %>% 
  data.frame()

distinct_high_conf_cORF_merged <-
  left_join(gtf_annot, high_conf_cORF_recoded_unique)

merged_gff_df <- df_to_gff(distinct_high_conf_cORF_merged)

writeGFF(merged_gff_df, snakemake@output[[1]])

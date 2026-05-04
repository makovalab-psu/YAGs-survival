library(dplyr)

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

# protein_list <- read.delim("C:/Users/gresh/Downloads/temp/human_apes_Y_NCBI.tsv", header=TRUE)
# high_conf_cORF <- read.table("C:/Users/gresh/Downloads/temp/gor_techrep1_Y_homologs.Ychr.txt", header=FALSE)
# output_file <- "C:/Users/gresh/Downloads/temp/gor_gene_fam_ids.txt"
# gene_fam <- 'BPY2'

protein_list <- read.delim(commandArgs(trailingOnly = TRUE)[1], header=TRUE)
high_conf_cORF <- read.table(commandArgs(trailingOnly = TRUE)[2], header=FALSE) 
output_file <- commandArgs(trailingOnly = TRUE)[3]
gene_fam <- commandArgs(trailingOnly = TRUE)[4]

high_conf_cORF_recoded <- recode_cORF(high_conf_cORF)

high_conf_cORF_recoded_gene_fam <-
  high_conf_cORF_recoded %>% 
  filter(pident >= 80) %>% 
  filter(qcovs >= 75) %>% 
  group_by(sseqid, transcript_id) %>% 
  filter(length == max(length)) %>% 
  filter(pident == max(pident)) %>% 
  filter(qcovs == max(qcovs)) %>% 
  ungroup() %>% 
  data.frame() %>% 
  filter(sseqid == gene_fam) %>% 
  dplyr::select(c(transcript_id, qseqid, sseqid)) %>% 
  distinct()

write.table(high_conf_cORF_recoded_gene_fam$qseqid, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)

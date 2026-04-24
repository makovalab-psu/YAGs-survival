library(microseq)
library(dplyr)
library(tidyr)
library(stringr)

#from 'Attributes' column of a gff/gtf file create separate columns; returns a dataframe
gff_attr <- function(gff_file) {
  
  gff_file_df <- data.frame(gff_file)
  gff_file_inter <- strsplit(gff_file_df$Attributes, split=";")
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

#from columns of dataframe create a gff file: take columns c(1:8) without changes,
#values in the rest of the columns merge using ';' and set column name as 'Attributes'
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

parse_attributes <- function(data) {
  
  # Assume the attributes column is named "attributes"
  attributes_col <- "Attributes"
  
  data <-
    data %>% 
    filter(!Type %in% c('region', 'CDS', 'V_gene_segment', 'cDNA_match')) 
  
  # Parse the attributes column into key-value pairs and reshape into columns
  parsed_attributes <- data %>% 
    select(!!sym(attributes_col)) %>%
    mutate(temp_attributes = str_split(!!sym(attributes_col), ";")) %>%
    unnest(temp_attributes) %>%
    mutate(
      key = str_extract(temp_attributes, "^[^=]+"),
      value = str_extract(temp_attributes, "(?<=\\=).*")
    ) %>%
    select(-temp_attributes) %>%
    pivot_wider(names_from = key, values_from = value)
  
  # Combine original first 8 columns with parsed attributes
  result <- data %>%
    select(1:8) %>%
    bind_cols(parsed_attributes)
  
  return(result)
  
}

################################################################################

gff_annot <- microseq::readGFF(snakemake@input[[1]])
gff_annot_grange <- microseq::readGFF(snakemake@input[[2]])
lastz <- read.table(snakemake@input[[3]], sep="\t")

# gff_annot <- microseq::readGFF("C:/Users/gresh/Downloads/temp/gor_techrep1.stringtie.Ychr.gtf")
# gff_annot_grange <- microseq::readGFF("C:/Users/gresh/Downloads/temp/gor_techrep1.stringtie.transcripts.annotation.Ychr.gtf")
# lastz <- read.table("C:/Users/gresh/Downloads/temp/gor_lastz.upd.Ychr.lastz", sep="\t")

gff_annot$Attributes <- gsub("\"", "", gff_annot$Attributes)
gff_annot <- gff_attr(gff_annot)

gff_annot_grange$Attributes <- gsub("\"", "", gff_annot_grange$Attributes)
gff_annot_grange <- gff_attr(gff_annot_grange)

################################################################################

# human <- read.csv2("C:/Users/gresh/Downloads/temp/human_gene_list.csv", sep=',')
# bonobo <- read.csv2("C:/Users/gresh/Downloads/temp/bonobo_gene_list.csv", sep=',')
# chimp <- read.csv2("C:/Users/gresh/Downloads/temp/chimp_gene_list.csv", sep=',')
# sum_orang <- read.csv2("C:/Users/gresh/Downloads/temp/sum_orang_gene_list.csv", sep=',')
# bor_orang <- read.csv2("C:/Users/gresh/Downloads/temp/bor_orang_gene_list.csv", sep=',')
# gor <- read.csv2("C:/Users/gresh/Downloads/temp/gor_gene_list.csv", sep=',')

human <- read.csv2(snakemake@input[[4]], sep=',')
bonobo <- read.csv2(snakemake@input[[5]], sep=',')
chimp <- read.csv2(snakemake@input[[6]], sep=',')
sum_orang <- read.csv2(snakemake@input[[7]], sep=',')
bor_orang <- read.csv2(snakemake@input[[8]], sep=',')
gor <- read.csv2(snakemake@input[[9]], sep=',')

all_species <- rbind(human, bonobo, chimp, sum_orang, bor_orang, gor)

################################################################################

gff_annot_grange_transcript_list <- unique(gff_annot_grange$transcript_id)

gff_annot_full <- data.frame()

for (i in 1:length(gff_annot_grange_transcript_list)){
  
  gff_annot_grange_transcript_list <- gff_annot_grange$transcript_id
  
  gff_annot_inter <-
    gff_annot %>% 
    filter(transcript_id == gff_annot_grange_transcript_list[i])
  
  gff_annot_grange_inter <-
    gff_annot_grange %>% 
    filter(transcript_id == gff_annot_grange_transcript_list[i])
  
  gff_annot_inter$ref_gene_id <- gff_annot_grange_inter$gene_ids
  gff_annot_inter$ref_biotype <- gff_annot_grange_inter$ref_gene_biotype
  gff_annot_inter$ref_family <- gff_annot_grange_inter$ref_gene_family
  
  gff_annot_full <- rbind(gff_annot_full, gff_annot_inter)
  
}

gff_annot_full_gff <- df_to_gff(gff_annot_full)

################################################################################

recode_dict <- setNames(all_species$gene_family, all_species$gene.name)

lastz <- 
  lastz %>%
  mutate(lastz_family = recode(lastz$V4, !!!recode_dict))

################################################################################

lastz$V8 <- as.numeric(sub("%", "", lastz$V8))
lastz$V11 <- as.numeric(sub("%", "", lastz$V11))

lastz_full <-
  lastz %>% 
  dplyr::select(c(V1, V4, V8, V11, lastz_family)) %>% 
  distinct()

lastz_fltr_full_max <-
  lastz %>% 
  filter(V11 > 80) %>% 
  filter(V8 > 80) %>% 
  dplyr::select(c(V1, V4, V8, V11, lastz_family)) %>% 
  group_by(V1) %>% 
  slice_max(order_by = V11) %>%
  slice_max(order_by = V8) %>%
  ungroup() %>% 
  data.frame() %>% 
  select(c(V1, V8, V11, lastz_family)) %>% 
  distinct() %>% 
  filter(lastz_family %in% c('HSFY', 'TSPY', 'RBMY', 'VCY', 'CDY', 'PRY', 'DAZ', 'BPY2', 'XKRY'))

colnames(lastz_fltr_full_max) <- c('transcript_id', 'lastz_prc', 'lastz_cov', 'lastz_annot')
  
lastz_fltr_full_max <- 
  lastz_fltr_full_max %>% 
  group_by(transcript_id) %>%
  filter(n() == 1) %>% 
  ungroup() %>% 
  data.frame()
  
################################################################################

gff_annot_full_lastz <-
  left_join(gff_annot_full, lastz_fltr_full_max)

################################################################################

# gtf_bonobo_Ychr <- readGFF("C:/Users/gresh/Downloads/chrY_bonobo.gff")
# gtf_human_Ychr <- readGFF("C:/Users/gresh/Downloads/chrY_human.gff")
# gtf_gor_Ychr <- readGFF("C:/Users/gresh/Downloads/chrY_gor.gff")
# gtf_chimp_Ychr <- readGFF("C:/Users/gresh/Downloads/chrY_chimp.gff")
# gtf_bor_orang_Ychr <- readGFF("C:/Users/gresh/Downloads/chrY_bor_orang.gff")
# gtf_sum_orang_Ychr <- readGFF("C:/Users/gresh/Downloads/chrY_sum_orang.gff")

gtf_bonobo_Ychr <- readGFF(snakemake@input[[10]])
gtf_human_Ychr <- readGFF(snakemake@input[[11]])
gtf_gor_Ychr <- readGFF(snakemake@input[[12]])
gtf_chimp_Ychr <- readGFF(snakemake@input[[13]])
gtf_bor_orang_Ychr <- readGFF(snakemake@input[[14]])
gtf_sum_orang_Ychr <- readGFF(snakemake@input[[15]])

################################################################################
#CREATE DICTIONARY WITH ALL GENES AND TRANSCRIPTS

gtf_bonobo_Ychr_Attributes <- parse_attributes(gtf_bonobo_Ychr)

gtf_bonobo_Ychr_Attributes_YAGs <- 
  gtf_bonobo_Ychr_Attributes %>% 
  filter(gene %in% bonobo$gene.name)

gtf_bonobo_Ychr_Attributes_YAGs_gene_transcripts <- 
  gtf_bonobo_Ychr_Attributes_YAGs %>% 
  filter(Type != 'exon')

bonobo_dictionary <-
  gtf_bonobo_Ychr_Attributes_YAGs_gene_transcripts %>% 
  data.frame() %>% 
  select(ID, gene) %>% 
  mutate(lastz_family = recode(gtf_bonobo_Ychr_Attributes_YAGs_gene_transcripts$gene, !!!recode_dict)) %>% 
  select(-c(gene))

################################################################################

gtf_human_Ychr_Attributes <- parse_attributes(gtf_human_Ychr)

gtf_human_Ychr_Attributes_YAGs <- 
  gtf_human_Ychr_Attributes %>% 
  filter(gene %in% human$gene.name)

gtf_human_Ychr_Attributes_YAGs_gene_transcripts <- 
  gtf_human_Ychr_Attributes_YAGs %>% 
  filter(Type != 'exon')

human_dictionary <-
  gtf_human_Ychr_Attributes_YAGs_gene_transcripts %>% 
  data.frame() %>% 
  select(ID, gene) %>% 
  mutate(lastz_family = recode(gtf_human_Ychr_Attributes_YAGs_gene_transcripts$gene, !!!recode_dict)) %>% 
  select(-c(gene))

################################################################################

gtf_gor_Ychr_Attributes <- parse_attributes(gtf_gor_Ychr)

gtf_gor_Ychr_Attributes_YAGs <- 
  gtf_gor_Ychr_Attributes %>% 
  filter(gene %in% gor$gene.name)

gtf_gor_Ychr_Attributes_YAGs_gene_transcripts <- 
  gtf_gor_Ychr_Attributes_YAGs %>% 
  filter(Type != 'exon')

gor_dictionary <-
  gtf_gor_Ychr_Attributes_YAGs_gene_transcripts %>% 
  data.frame() %>% 
  select(ID, gene) %>% 
  mutate(lastz_family = recode(gtf_gor_Ychr_Attributes_YAGs_gene_transcripts$gene, !!!recode_dict)) %>% 
  select(-c(gene))

################################################################################

gtf_chimp_Ychr_Attributes <- parse_attributes(gtf_chimp_Ychr)

gtf_chimp_Ychr_Attributes_YAGs <- 
  gtf_chimp_Ychr_Attributes %>% 
  filter(gene %in% chimp$gene.name)

gtf_chimp_Ychr_Attributes_YAGs_gene_transcripts <- 
  gtf_chimp_Ychr_Attributes_YAGs %>% 
  filter(Type != 'exon')

chimp_dictionary <-
  gtf_chimp_Ychr_Attributes_YAGs_gene_transcripts %>% 
  data.frame() %>% 
  select(ID, gene) %>% 
  mutate(lastz_family = recode(gtf_chimp_Ychr_Attributes_YAGs_gene_transcripts$gene, !!!recode_dict)) %>% 
  select(-c(gene))

################################################################################

gtf_bor_orang_Ychr_Attributes <- parse_attributes(gtf_bor_orang_Ychr)

gtf_bor_orang_Ychr_Attributes_YAGs <- 
  gtf_bor_orang_Ychr_Attributes %>% 
  filter(gene %in% bor_orang$gene.name)

gtf_bor_orang_Ychr_Attributes_YAGs_gene_transcripts <- 
  gtf_bor_orang_Ychr_Attributes_YAGs %>% 
  filter(Type != 'exon')

bor_orang_dictionary <-
  gtf_bor_orang_Ychr_Attributes_YAGs_gene_transcripts %>% 
  data.frame() %>% 
  select(ID, gene) %>% 
  mutate(lastz_family = recode(gtf_bor_orang_Ychr_Attributes_YAGs_gene_transcripts$gene, !!!recode_dict)) %>% 
  select(-c(gene))

################################################################################

gtf_sum_orang_Ychr_Attributes <- parse_attributes(gtf_sum_orang_Ychr)

gtf_sum_orang_Ychr_Attributes_YAGs <- 
  gtf_sum_orang_Ychr_Attributes %>% 
  filter(gene %in% sum_orang$gene.name)

gtf_sum_orang_Ychr_Attributes_YAGs_gene_transcripts <- 
  gtf_sum_orang_Ychr_Attributes_YAGs %>% 
  filter(Type != 'exon')

sum_orang_dictionary <-
  gtf_sum_orang_Ychr_Attributes_YAGs_gene_transcripts %>% 
  data.frame() %>% 
  select(ID, gene) %>% 
  mutate(lastz_family = recode(gtf_sum_orang_Ychr_Attributes_YAGs_gene_transcripts$gene, !!!recode_dict)) %>% 
  select(-c(gene))

################################################################################

all_species_dictionary <-
  rbind(sum_orang_dictionary, bor_orang_dictionary, gor_dictionary, chimp_dictionary, human_dictionary, bonobo_dictionary)

all_species_dictionary_dict <- setNames(all_species_dictionary$lastz_family, all_species_dictionary$ID)

################################################################################
# potential_novel_gene_copies

gff_annot_full_lastz_potentially_novel <-
  gff_annot_full_lastz %>% 
  filter(!is.na(lastz_annot)) %>% 
  filter(is.na(ref_gene_id))

# additional check - lastz against all Y chr transcripts
# lastz_all_Y_transcripts <- read.table("C:/Users/gresh/Downloads/temp/lastz_all_Y_transcripts/gor_lastz.Ychr.lastz", sep="\t")
lastz_all_Y_transcripts <- read.table(snakemake@input[[16]], sep="\t")

lastz_all_Y_transcripts$V8 <- as.numeric(gsub("%", "", lastz_all_Y_transcripts$V8))
lastz_all_Y_transcripts$V11 <- as.numeric(gsub("%", "", lastz_all_Y_transcripts$V11))

lastz_all_Y_transcripts_filter <-
  lastz_all_Y_transcripts %>% 
  filter(V1 %in% unique(gff_annot_full_lastz_potentially_novel$transcript_id)) %>% 
  filter(V4 %in% all_species_dictionary$ID) %>% 
  filter(V11 > 80) %>% 
  filter(V8 > 80) %>% 
  group_by(V1) %>% 
  filter(V11 == max(V11)) %>% 
  filter(V8 == max(V8)) %>% 
  ungroup()

lastz_all_Y_transcripts_filter$lastz_fam <- recode(lastz_all_Y_transcripts_filter$V4, !!!all_species_dictionary_dict)

Y_chr_novel_YAG_copies <-
  lastz_all_Y_transcripts_filter %>% 
  select(c(V1, V8, V11, lastz_fam)) %>% 
  distinct() %>% 
  data.frame()

################################################################################

#identify novel_gene_copies that are not verified by lastz against all Y transcripts
gff_annot_full_lastz_novel_upd <-
  gff_annot_full_lastz_potentially_novel %>% 
  filter(!transcript_id %in% Y_chr_novel_YAG_copies$V1)

#identify transcripts with inconsistent annotation
gff_annot_full_mltp <-
  gff_annot_full_lastz %>% 
  filter(lastz_annot != ref_family)

YAGs <- c('BPY2', 'PRY', 'CDY', 'XKRY', 'DAZ', 'HSFY', 'RBMY', 'TSPY', 'VCY')

gff_annot_full_YAG_clear <- 
  gff_annot_full_lastz %>% 
  filter(lastz_annot %in% YAGs | ref_family %in% YAGs) %>% 
  filter(!transcript_id %in% gff_annot_full_lastz_novel_upd$transcript_id) %>% 
  filter(!transcript_id %in% gff_annot_full_mltp$transcript_id)

gff_annot_full_lastz_novel_vf <-
  gff_annot_full_lastz_potentially_novel %>% 
  filter(transcript_id %in% Y_chr_novel_YAG_copies$V1)

microseq::writeGFF(tibble(df_to_gff(gff_annot_full_YAG_clear)), snakemake@output[[1]])

microseq::writeGFF(tibble(df_to_gff(gff_annot_full_mltp)), snakemake@output[[2]])

microseq::writeGFF(tibble(df_to_gff(gff_annot_full_lastz_novel_vf)), snakemake@output[[3]])

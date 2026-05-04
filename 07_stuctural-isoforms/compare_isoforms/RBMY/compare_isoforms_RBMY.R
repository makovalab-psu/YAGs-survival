library(tidyverse)
library(Biostrings)
library(ggpubr)
library(gplots)
library(microseq)

################################################################################

#create summary table
tmap <- data.frame(matrix(ncol=3, nrow=0))
colnames(tmap) <- c('ref_id', 'class_code', 'qry_gene_id')

for (i in 1:167){
  
  filename <- paste("C:/Users/gresh/Downloads/temp/RBMY/transcript_", i, 
                    ".RBMY_transcripts.gff.tmap", sep='')
  
  tmap_inter <- read.table(filename, header = TRUE)
  
  tmap_inter$class_code <- ifelse(tmap_inter$class_code == '=', 1, 0)
  
  tmap_inter_sub <-
    tmap_inter %>% 
    select(c(ref_id , class_code, qry_gene_id))
  
  tmap <- rbind(tmap, tmap_inter_sub)
  
}

#make a copy of a column with original transcript id
tmap$qry_original <- tmap$qry_gene_id

#add column with species
your_list = tmap$qry_gene_id

your_list <- gsub('bor_orang', 'bororang', your_list)
your_list <- gsub('sum_orang', 'sumorang', your_list)

tmap$qry_species <- sub(".*_(\\w+)_techrep1_.*", "\\1", your_list)

tmap$qry_species <- gsub('bororang', 'bor_orang', tmap$qry_species)
tmap$qry_species <- gsub('sumorang', 'sum_orang', tmap$qry_species)

#change transcript names to isoform names
transcript_names <- unique(tmap$ref_id)

tmap %>%
  group_by(ref_id) %>%
  mutate(isoform = cumsum(class_code == 1 & lag(class_code, default = 0) == 0)) %>% 
  as.data.frame()

for (i in 1:length(transcript_names)){
  
  isoform_name <- paste('isoform', i, sep='_')
  
  transcript_name <- transcript_names[i]
  
  same_isoform_list <- 
    tmap[which(tmap$ref_id == transcript_name & 
                 tmap$class_code == '1'), 3]
  
  tmap$ref_id <- ifelse(tmap$ref_id %in% same_isoform_list, 
                        isoform_name, tmap$ref_id)
  
  tmap$qry_gene_id <- ifelse(tmap$qry_gene_id %in% same_isoform_list, 
                             isoform_name, tmap$qry_gene_id)
  
}

################################################################################
#sequence isoforms

# Load your multiple sequence alignment (MSA)
alignment <- 
  readAAStringSet("C:/Users/gresh/Downloads/temp/RBMY/RBMY_algn.fa")

names(alignment) <- gsub("_extracted.*$", "", names(alignment))

# Create a function to calculate pairwise coverage
calculate_pairwise_coverage <- function(seq1, seq2) {
  # Convert sequences to character vectors
  seq1_chars <- unlist(strsplit(as.character(seq1), ""))
  seq2_chars <- unlist(strsplit(as.character(seq2), ""))
  
  # Calculate the number of positions with any differences
  different_positions <- sum(seq1_chars != seq2_chars)
  
  # Calculate coverage as a fraction
  coverage <- different_positions / length(seq1_chars)
  
  return(coverage)
}

# Create a matrix to store pairwise coverage values
num_seqs <- length(alignment)
coverage_matrix <- matrix(0, nrow = num_seqs, ncol = num_seqs)

# Assign sequence IDs as row and column names
row.names(coverage_matrix) <- names(alignment)
colnames(coverage_matrix) <- names(alignment)

# Calculate pairwise coverage for all pairs of sequences
for (i in 1:num_seqs) {
  for (j in 1:num_seqs) {
    coverage_matrix[i, j] <- calculate_pairwise_coverage(alignment[i], alignment[j])
  }
}

sorted_indices <- order(names(alignment))
coverage_matrix <- coverage_matrix[sorted_indices, sorted_indices]

# Create a lower diagonal matrix
lower_diagonal_matrix <- coverage_matrix

# Replace elements above the main diagonal with zeros
# lower_diagonal_matrix[upper.tri(coverage_matrix)] <- NA

seq_comparison_df <- data.frame(lower_diagonal_matrix)

seq_comparison_df <- rownames_to_column(seq_comparison_df, var = "ref_orf")

seq_comparison_df_long <-
  pivot_longer(seq_comparison_df,
               cols = -ref_orf, 
               names_to = "qry_orf", 
               values_to = "seqid") %>% 
  filter(!is.na(seqid)) %>% 
  mutate(similarity = ifelse(seqid == 0, '=', 'n'))

#fix the transcript ids
# seq_comparison_df_long$ref_orf <- gsub(' ', '_', seq_comparison_df_long$ref_orf)
# seq_comparison_df_long$ref_orf <- gsub('\\.', '_', seq_comparison_df_long$ref_orf)
# seq_comparison_df_long$qry_orf <- gsub('\\.', '_', seq_comparison_df_long$qry_orf)

seq_comparison_df_long$ref_orf_original <- 
  seq_comparison_df_long$ref_orf

seq_comparison_df_long$qry_orf_original <- 
  seq_comparison_df_long$qry_orf

#change transcript names to isoform names
transcript_names <- unique(seq_comparison_df_long$qry_orf)

for (i in 1:length(transcript_names)){
  
  isoform_name <- paste('isoform', i, sep='_')
  
  transcript_name <- transcript_names[i]
  
  same_isoform_list <- 
    seq_comparison_df_long[which(seq_comparison_df_long$qry_orf == transcript_name & 
                                   seq_comparison_df_long$similarity == '='), 1]
  
  same_isoform_list <- same_isoform_list$ref_orf
  
  seq_comparison_df_long$ref_orf <- ifelse(seq_comparison_df_long$ref_orf %in% same_isoform_list, 
                                           isoform_name, seq_comparison_df_long$ref_orf)
  
  seq_comparison_df_long$qry_orf  <- ifelse(seq_comparison_df_long$qry_orf  %in% same_isoform_list, 
                                            isoform_name, seq_comparison_df_long$qry_orf )
  
}

################################################################################

sequence_isoforms <- seq_comparison_df_long[, c('qry_orf', 
                                                'qry_orf_original')]

sequence_isoforms <- data.frame(distinct(sequence_isoforms))

colnames(sequence_isoforms) <- c('sequence_isoform', 'transcript_id')

structural_isoforms <- tmap[,c('qry_gene_id', 'qry_original', 'qry_species')]

structural_isoforms <- data.frame(distinct(structural_isoforms))

# structural_isoforms$qry_original <- gsub('\\.', '_', structural_isoforms$qry_original)

colnames(structural_isoforms) <- c('structural_isoform', 'transcript_id', 'species')

all_isoforms <-
  merge(structural_isoforms, sequence_isoforms)

################################################################################
#fix isofrom names

original_values = unique(all_isoforms$structural_isoform)

mapping <- setNames(paste0("isoform_", 1:length(original_values)), original_values)

all_isoforms$structural_isoform <- mapping[all_isoforms$structural_isoform]

original_values = unique(all_isoforms$sequence_isoform)

mapping <- setNames(paste0("isoform_", 1:length(original_values)), original_values)

all_isoforms$sequence_isoform <- mapping[all_isoforms$sequence_isoform]

write.csv(all_isoforms, 
          "C:/Users/gresh/Downloads/temp/RBMY/all_isoforms.csv", 
          row.names = FALSE)

################################################################################
#plot structural isoforms

isoform_summary <-
  all_isoforms %>% 
  dplyr::select(c(transcript_id, structural_isoform, species)) %>% 
  distinct() %>% 
  group_by(structural_isoform, species) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

isoform_summary_wide <-
  isoform_summary %>% 
  pivot_wider(names_from=species,
              values_from=n)

#manually adding species that are not present by defalut
# isoform_summary_wide$bonobo <- 0
# isoform_summary_wide$chimp <- 0
# isoform_summary_wide$gor <- 0
# isoform_summary_wide$human <- 0
# isoform_summary_wide$bor_orang <- 0
# isoform_summary_wide$sum_orang <- 0

isoform_summary_mod <-
  isoform_summary_wide %>% 
  mutate_if(is.numeric, ~coalesce(., 0)) %>% 
  pivot_longer(cols = -structural_isoform, names_to = "species", values_to = "n")

colnames(isoform_summary_mod)[1] <- 'isoform_id'

isoform_summary_mod$isoform_id <- gsub('isoform_', '', isoform_summary_mod$isoform_id)

isoform_summary_mod <- 
  isoform_summary_mod %>% 
  mutate(species = fct_recode(species,
                             "chimpanzee" = "chimp",
                             "gorilla" = "gor",
                             "Bornean orangutan" = "bor_orang",
                             "Sumatran orangutan" = "sum_orang"))

plot <-
  ggplot(isoform_summary_mod, aes(species, isoform_id, width = 1, height = 1)) +
  geom_tile(aes(fill = n), colour = "white") +
  geom_text(aes(label = ifelse(n == 0, NA, n)), 
            vjust = 0.5, hjust = 0.5, size = 3, alpha = 0.75)+
  scale_fill_gradient(low = "white", high = "red",
                      limits = c(0, 42))+
  scale_x_discrete(limits=c('bonobo', 'chimpanzee', 'gorilla', 
                            'Bornean orangutan', 'Sumatran orangutan', 'human'))+
  scale_y_discrete(limits=
                     paste0(1:length(unique(isoform_summary$structural_isoform))),
                   position = "right")+
  labs(title = expression(bolditalic("RBMY")),subtitle = "structural isoforms") +
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5,size = 16 ),
        plot.subtitle = element_text(hjust = 0.5, size = 12,face = "plain", color = "red"))

ggsave("C:/Users/gresh/Downloads/temp/RBMY/plot_structural_isoforms.png", 
       plot, width = 6, height = 30, units = "cm")

################################################################################
#plot sequence isoforms

isoform_summary <-
  all_isoforms %>% 
  dplyr::select(c(transcript_id, sequence_isoform, species)) %>% 
  distinct() %>% 
  group_by(sequence_isoform, species) %>% 
  summarise(n = n()) %>% 
  as.data.frame()

isoform_summary_wide <-
  isoform_summary %>% 
  pivot_wider(names_from=species,
              values_from=n)

#manually adding species that are not present by defalut
# isoform_summary_wide$bonobo <- 0
# isoform_summary_wide$chimp <- 0
# isoform_summary_wide$gor <- 0
# isoform_summary_wide$human <- 0
# isoform_summary_wide$bor_orang <- 0
# isoform_summary_wide$sum_orang <- 0

isoform_summary_mod <-
  isoform_summary_wide %>% 
  mutate_if(is.numeric, ~coalesce(., 0)) %>% 
  pivot_longer(cols = -sequence_isoform, names_to = "species", values_to = "n")

colnames(isoform_summary_mod)[1] <- 'isoform_id'

isoform_summary_mod$isoform_id <- gsub('isoform_', '', isoform_summary_mod$isoform_id)

isoform_summary_mod <- 
  isoform_summary_mod %>% 
  mutate(species = fct_recode(species,
                              "chimpanzee" = "chimp",
                              "gorilla" = "gor",
                              "Bornean orangutan" = "bor_orang",
                              "Sumatran orangutan" = "sum_orang"))

plot <-
  ggplot(isoform_summary_mod, aes(species, isoform_id,
                                  width = 1, height = 1)) +
  
  geom_tile(aes(fill = n), colour = "white") +
  geom_text(aes(label = ifelse(n == 0, NA, n)), 
            vjust = 0.5, hjust = 0.5, size = 3, alpha = 0.75)+
  scale_fill_gradient(low = "white", high = "blue",
                      limits = c(0, 25))+
  scale_x_discrete(limits=c('bonobo', 'chimpanzee', 'gorilla', 
                            'Bornean orangutan', 'Sumatran orangutan', 'human'))+
  scale_y_discrete(limits=
                     paste0(1:length(unique(isoform_summary$sequence_isoform))),
                   position = "right")+
  labs(title = expression(bolditalic("RBMY")),subtitle = "sequence isoforms") +
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5,size = 16 ),
        plot.subtitle = element_text(hjust = 0.5, size = 12,face = "plain", color = "blue"       
        ))

ggsave("C:/Users/gresh/Downloads/temp/RBMY/plot_sequence_isoforms.png", 
       plot, width = 6, height = 55, units = "cm")

################################################################################
#relationship between sequence isoforms

#make isoform names consitent

mapping <- setNames(all_isoforms$sequence_isoform, all_isoforms$transcript_id)


seq_comparison_df_long$ref_orf <- 
  seq_comparison_df_long$ref_orf_original

seq_comparison_df_long$ref_orf <- mapping[seq_comparison_df_long$ref_orf]


seq_comparison_df_long$qry_orf <- 
  seq_comparison_df_long$qry_orf_original

seq_comparison_df_long$qry_orf <- mapping[seq_comparison_df_long$qry_orf]

df_seqid <-
  seq_comparison_df_long %>% 
  dplyr::select(c(ref_orf, qry_orf, seqid)) %>% 
  distinct %>% 
  pivot_wider(names_from = qry_orf, values_from=seqid) %>% 
  data.frame()

rownames(df_seqid) <- df_seqid$ref_orf

df_seqid <- df_seqid[,-1]

png("C:/Users/gresh/Downloads/gffcompare/RBMY/plot_seqid.png", width = 800, height = 800, units = "px")

heatmap.2(as.matrix(df_seqid), 
          trace = "none",  # Disable row and column dendrograms
          dendrogram = 'row',
          symm=TRUE,
          colsep = c(1:19),
          rowsep = c(1:19),
          sepcolor="white",
          sepwidth=c(0.005,0.005),
          key = FALSE,
          col = colorRampPalette(c("blue", "white", "red"))(20)
)

dev.off()

write.csv(df_seqid, 
          "C:/Users/gresh/Downloads/gffcompare/RBMY/isoform_seqid.csv", 
          row.names = TRUE)

################################################################################
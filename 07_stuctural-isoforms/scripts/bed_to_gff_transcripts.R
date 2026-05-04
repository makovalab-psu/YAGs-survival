library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)

# GENE_FAM <- 'XKRY'
# bed_data <- read.table(paste("C:/Users/gresh/Downloads/temp/uLTRA/", GENE_FAM, "_ultra/reads.sorted.bed", sep=''), sep = "\t", header = FALSE)
# output_path <- paste("C:/Users/gresh/Downloads/temp/", GENE_FAM,"/VCY_transcripts.gff", sep='')
# output_directory <- paste("C:/Users/gresh/Downloads/temp/", GENE_FAM,"/", sep='')
# dir.create(output_directory)

GENE_FAM <- commandArgs(trailingOnly = TRUE)[1]
bed_data <- read.table(commandArgs(trailingOnly = TRUE)[2], sep = "\t", header = FALSE)
output_path <- commandArgs(trailingOnly = TRUE)[3]
output_directory <- commandArgs(trailingOnly = TRUE)[4]

colnames(bed_data) <- c('name', 'start', 'end', 'transcript_name', 'score', 'strand')
species <- c("gor", "chimp", "sum_orang", "bor_orang", "human", "bonobo")
pattern <- paste(species, collapse = "|")
bed_data$species <- str_extract(bed_data$transcript_name, pattern)

# Step 1: Create transcript features
transcripts <- bed_data %>%
  group_by(name, transcript_name, strand, species) %>%
  summarize(
    start = min(start) + 1,  # Convert to 1-based indexing
    end = max(end),
    type = "transcript",
    score = unique(score),
    .groups = "drop"
  )

# Step 2: Add attributes as separate columns
transcripts <- transcripts %>%
  mutate(
    ID = transcript_name,
    Name = transcript_name,
    Species = species
  )

exons <- bed_data %>%
  mutate(
    type = "exon",
    start = start + 1,  # Convert to 1-based indexing
    ID = NA,  # Exons do not have unique IDs
    Name = transcript_name,
    Species = species,
    Parent = transcript_name,
    score = score
  )

# Step 3: Combine transcripts and exons
combined <- bind_rows(
  transcripts %>%
    mutate(Parent = NA),  # Transcripts do not have a parent
  exons
) %>%
  arrange(transcript_name, start)

combined_df <-
  combined %>%
  data.frame() %>%
  dplyr::rename(
    seqid = name,  # Rename "a" to "col1"
    type = type,
    start = start,
    end = end,
    score = score,
    strand = strand
  ) %>%
  mutate(Phase = '.') %>%
  mutate(Source = 'bed_to_gff')

export(GenomicRanges::GRanges(combined_df,
                              type = combined_df$type,
                              strand = combined_df$strand,
                              Parent = combined_df$Parent,
                              ID = combined_df$ID,
                              Name = combined_df$Name),
       output_path, format = "gff3")

###########################################################################################################################
#transform bed to gff and extract information about each transcript into a separate file

split_gff <- combined_df %>%
  group_by(Name) %>%
  group_split()

for (i in 1:length(split_gff)) {
  
  transcript_data <- split_gff[[i]]
  
  output_filename <- paste0(output_directory, "transcript_", i, ".gff")
  
  export(GenomicRanges::GRanges(transcript_data,
                                type = transcript_data$type,
                                strand = transcript_data$strand,
                                Parent = transcript_data$Parent,
                                ID = transcript_data$ID,
                                Name = transcript_data$Name), output_filename, format = "gff3")
  
}


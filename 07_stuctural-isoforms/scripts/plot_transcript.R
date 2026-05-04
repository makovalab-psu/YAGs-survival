library(ggplot2)
library(dplyr)
library(ggtranscript)
library(ggpubr)

GENE_FAM <- commandArgs(trailingOnly = TRUE)[1]
bed_data <- read.table(commandArgs(trailingOnly = TRUE)[2], sep = "\t", header = FALSE)
output_path <- commandArgs(trailingOnly = TRUE)[3]

# GENE_FAM <- 'HSFY'
# bed_data <- read.table("C:/Users/gresh/Downloads/temp/uLTRA/HSFY_ultra/reads.sorted.bed", sep = "\t", header = FALSE)
# output_path <- commandArgs(trailingOnly = TRUE)[3]

colnames(bed_data) <- c('chromosome', 'start', 'end', 'name', 'score', 'strand')

bed_data$name <- gsub('bor_orang', 'bororang', bed_data$name)
bed_data$name <- gsub('sum_orang', 'sumorang', bed_data$name)

split_data <- strsplit(bed_data$name, "_")

all_params <- as.data.frame(do.call(rbind, split_data))

bed_data$species <- all_params$V3
bed_data$type <- 'exon'

bed_data$species <- gsub('bororang', 'bor_orang', bed_data$species)
bed_data$species <- gsub('sumorang', 'sum_orang', bed_data$species)

bed_data$species <- factor(bed_data$species, 
                           levels = c('human', 'chimp', 'bonobo', 'gor', 'bor_orang', 'sum_orang'))

bed_data$transcript_name <- paste(all_params$V3, all_params$V1, all_params$V2, sep='_')

x_axis_label <- paste("human", GENE_FAM, "gene coordinates", sep=' ') 

bed_plot <-
  bed_data %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = species)
  ) + labs(x = x_axis_label, y = "Transcripts") +
  geom_intron(
    data = to_intron(bed_data, "transcript_name"),
    aes(strand = strand) 
  )+
  theme_pubclean()+
  scale_fill_manual(values = c("human" = "#009E73",
                               "bonobo" = "#AA4499",
                               "chimp" = "#F0E442",
                               "gor" = "#0072B2",
                               "bor_orang" = "#D55E00",
                               "sum_orang" = "#888888"))

ggsave(output_path, bed_plot, scale=4)



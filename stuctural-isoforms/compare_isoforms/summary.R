library(vegan)
library(rstatix)
library(dplyr)
library(tidyr)

XKRY_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/XKRY/all_isoforms.csv")
XKRY_isoforms$gene_fam <- 'XKRY'

VCY_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/VCY/all_isoforms.csv")
VCY_isoforms$gene_fam <- 'VCY'

TSPY_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/TSPY/all_isoforms.csv")
TSPY_isoforms$gene_fam <- 'TSPY'

RBMY_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/RBMY/all_isoforms.csv")
RBMY_isoforms$gene_fam <- 'RBMY'

PRY_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/PRY/all_isoforms.csv")
PRY_isoforms$gene_fam <- 'PRY'

DAZ_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/DAZ/all_isoforms.csv")
DAZ_isoforms$gene_fam <- 'DAZ'

CDY_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/CDY/all_isoforms.csv")
CDY_isoforms$gene_fam <- 'CDY'

BPY2_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/BPY2/all_isoforms.csv")
BPY2_isoforms$gene_fam <- 'BPY2'

HSFY_isoforms <- 
  read.csv("C:/Users/gresh/Downloads/temp/HSFY/all_isoforms.csv")
HSFY_isoforms$gene_fam <- 'HSFY'

all_isoforms <-
  rbind(HSFY_isoforms, BPY2_isoforms, CDY_isoforms, DAZ_isoforms,
        PRY_isoforms, RBMY_isoforms, TSPY_isoforms, VCY_isoforms,
        XKRY_isoforms)

#############################################################################

calculate_diversity <- function(isoforms, gene_fam){
  
  data <-
    isoforms %>% 
    group_by(structural_isoform, species) %>% 
    summarise(n=n()) %>% 
    pivot_wider(names_from = structural_isoform,
                values_from = n) %>% 
    mutate_if(is.numeric, ~coalesce(., 0)) %>% 
    data.frame()
  
  rownames(data) <- data$species
  
  data <- data[,-1]
  
  structural_diversity <- data.frame(diversity(data, index = "invsimpson"))
  
  structural_diversity$species <- rownames(structural_diversity)
  
  structural_diversity$type <- 'structural'
  
  structural_diversity$gene_fam <- gene_fam
  
  ##########################################################
  
  data <-
    isoforms %>% 
    group_by(sequence_isoform, species) %>% 
    summarise(n=n()) %>% 
    pivot_wider(names_from = sequence_isoform,
                values_from = n) %>% 
    mutate_if(is.numeric, ~coalesce(., 0)) %>% 
    data.frame()
  
  rownames(data) <- data$species
  
  data <- data[,-1]
  
  sequence_diversity <- data.frame(diversity(data, index = "invsimpson"))
  
  sequence_diversity$species <- rownames(sequence_diversity)
  
  sequence_diversity$type <- 'sequence'
  
  sequence_diversity$gene_fam <- gene_fam
  
  diversity <- rbind(structural_diversity, 
                     sequence_diversity)
  
  return(diversity)
  
}

HSFY_diversity <- calculate_diversity(HSFY_isoforms, 'HSFY')
XKRY_diversity <- calculate_diversity(XKRY_isoforms, 'XKRY')
VCY_diversity <- calculate_diversity(VCY_isoforms, 'VCY')
TSPY_diversity <- calculate_diversity(TSPY_isoforms, 'TSPY')
RBMY_diversity <- calculate_diversity(RBMY_isoforms, 'RBMY')
PRY_diversity <- calculate_diversity(PRY_isoforms, 'PRY')
DAZ_diversity <- calculate_diversity(DAZ_isoforms, 'DAZ')
CDY_diversity <- calculate_diversity(CDY_isoforms, 'CDY')
BPY2_diversity <- calculate_diversity(BPY2_isoforms, 'BPY2')

PRY_diversity$species <- 'human'

all_gene_fam_diversity <-
  rbind(HSFY_diversity, XKRY_diversity, VCY_diversity,
        TSPY_diversity, RBMY_diversity, PRY_diversity,
        DAZ_diversity, CDY_diversity, BPY2_diversity)

rownames(all_gene_fam_diversity) <- 1:nrow(all_gene_fam_diversity)

colnames(all_gene_fam_diversity) <-
  c('shannon', 'species', 'type', 'gene_fam')

custom_colors = c("sequence" = "#1f77b4", "structural" =  "#ff7f0e")
custom_order <- c("bonobo", "chimpanzee", "gorilla",
                  "Bornean orangutan", "Sumatran orangutan", "human")

library('forcats')

all_gene_fam_diversity <- 
  all_gene_fam_diversity %>% 
  mutate(species = fct_recode(species,
                              'chimpanzee' = 'chimp',
                              'Bornean orangutan'  = 'bor_orang',
                              'Sumatran orangutan' = 'sum_orang',
                              'gorilla' = 'gor', 'human' = 'human',
                              'bonobo' = 'bonobo'))

all_gene_fam_diversity$species <-
  factor(all_gene_fam_diversity$species, custom_order)

all_gene_fam_diversity$shannon_round <- 
  round(all_gene_fam_diversity$shannon, 2)

diversity <-
  ggplot(all_gene_fam_diversity,
         aes(x=species, y=shannon_round, fill=type, label=shannon_round))+
  geom_col(position = position_dodge2(.9,  preserve = "single"),
           color='black', size=0.5)+
  geom_text(vjust=-0.5,
            position = position_dodge2(.9,  preserve = "single"), 
            size=2)+
  scale_fill_manual(values = custom_colors)+
  facet_wrap(vars(gene_fam))+
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.text = element_text(face = "bold.italic"))+
  ylab('Isoform diversity, inverse Simpson index')

ggsave("C:/Users/gresh/Downloads/temp/diversity.png", 
       diversity,
       width = 12, height = 8, dpi = 500)

#########################################################################

n_structural_isoform <-
  all_isoforms %>% 
  dplyr::select(c(gene_fam, species, structural_isoform)) %>% 
  distinct() %>% 
  group_by(gene_fam, species) %>% 
  summarise(n_structural = n()) %>% 
  as.data.frame()

n_sequence_isoform <-
  all_isoforms %>% 
  dplyr::select(c(gene_fam, species, sequence_isoform)) %>% 
  distinct() %>% 
  group_by(gene_fam, species) %>% 
  summarise(n_sequence = n()) %>% 
  as.data.frame()

isoform_count <-
  merge(n_structural_isoform, n_sequence_isoform)

########################################################################

isoform_count_long <-
  isoform_count %>% 
  pivot_longer(n_structural:n_sequence)

isoform_count_long <- 
  isoform_count_long %>% 
  mutate(species = fct_recode(species,
                              'chimpanzee' = 'chimp',
                              'Bornean orangutan'  = 'bor_orang',
                              'Sumatran orangutan' = 'sum_orang',
                              'gorilla' = 'gor', 'human' = 'human',
                              'bonobo' = 'bonobo'))

custom_colors = c("n_sequence" = "#1f77b4", "n_structural" =  "#ff7f0e")
custom_order <- c("bonobo", "chimpanzee", "gorilla",
                  "Bornean orangutan", "Sumatran orangutan", "human")

isoform_count_long$species <-
  factor(isoform_count_long$species, custom_order)

isoform_count_plot <-
  ggplot(isoform_count_long,
         aes(x=species, y=value, fill=name, label=value))+
  geom_col(position = position_dodge2(.9,  preserve = "single"),
           color='black', size=0.5)+
  geom_text(vjust=-0.5,
            position = position_dodge2(.9,  preserve = "single"), 
            size=2)+
  scale_fill_manual(values = custom_colors)+
  facet_wrap(vars(gene_fam))+
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(face = "bold.italic"))+
  scale_x_discrete(limits = custom_order)

ggsave("C:/Users/gresh/Downloads/temp/isoform_count.png", 
       isoform_count_plot,
       width = 12, height = 8, dpi = 500)

########################################################################

isoform_count_total <-
  isoform_count_long %>% 
  data.frame() %>% 
  pivot_wider(names_from=name, values_from=value) %>% 
  group_by(gene_fam) %>% 
  summarise(n_structural = sum(n_structural),
            n_sequence = sum(n_sequence)) %>% 
  pivot_longer(cols=n_structural:n_sequence)

plot_isoform_count <-
  isoform_count_total %>% 
  mutate(name=fct_recode(name,
                         'sequence isoforms' = 'n_sequence',
                         'structural isoforms' = 'n_structural')) %>% 
  ggplot(aes(x = gene_fam, y = value, group = name)) +
  geom_col(aes(fill = name),
           stat = "identity", 
           position = position_dodge(width = 1),
           color = 'black') +
  geom_text(aes(label = value), 
            position = position_dodge(width = 1), 
            vjust = -0.5, # Adjust vertical position above the bars
            size = 3) + # Adjust text size
  theme_pubclean()+
  scale_fill_manual(values=c("sequence isoforms" = "#1f77b4", 
                             "structural isoforms" =  "#ff7f0e"))+
  xlab('Ampliconic gene family')+
  ylab('Isoform count') +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(face = "bold.italic"))

ggsave("C:/Users/gresh/Downloads/temp/diversity.png", 
       diversity,
       width = 12, height = 8, dpi = 500)

plots_a_b <-
  ggarrange(plot_isoform_count, diversity,
          labels = c("A", "B"),  
          ncol = 1,              
          nrow = 2,
          heights = c(1, 3),
          align='hv')  

ggsave("C:/Users/gresh/Downloads/temp/count_diversity.png", 
       plots_a_b,
       width = 10, height = 13, dpi = 500)

########################################################################

diversity_ind <-
  all_gene_fam_diversity %>% 
  dplyr::select(-shannon_round) %>% 
  pivot_wider(names_from = type, values_from=shannon) %>% 
  as.data.frame()

isoform_count <- 
  isoform_count %>% 
  mutate(species = fct_recode(species,
                              'chimpanzee' = 'chimp',
                              'Bornean orangutan'  = 'bor_orang',
                              'Sumatran orangutan' = 'sum_orang',
                              'gorilla' = 'gor', 'human' = 'human',
                              'bonobo' = 'bonobo'))

isoform_count_diversity <-
  merge(isoform_count, diversity_ind)

expression <- read.csv("C:/Users/gresh/Downloads/expr.csv")

expression$exp <- as.numeric(expression$exp)

expression <- 
  expression %>% 
  mutate(species = fct_recode(species,
                              'chimpanzee' = 'chimp',
                              'Bornean orangutan'  = 'bor_orang',
                              'Sumatran orangutan' = 'sum_orang',
                              'gorilla' = 'gor', 'human' = 'human1',
                              'bonobo' = 'bonobo'))

mean_expression <-
  expression %>% 
  group_by(gene_fam, species) %>% 
  summarise(mean_exp = mean(exp, na.rm = TRUE)) %>% 
  as.data.frame()

expression_isoform_count <-
  merge(isoform_count_diversity, mean_expression)

expression_isoform_count$mean_exp_log <-
  log(expression_isoform_count$mean_exp)

gene_copies <- read.csv("C:/Users/gresh/Downloads/gene_copies_ncbi.csv")

gene_copies <- 
  gene_copies %>% 
  mutate(species = fct_recode(species,
                              'chimpanzee' = 'chimp',
                              'Bornean orangutan'  = 'bor_orang',
                              'Sumatran orangutan' = 'sum_orang',
                              'gorilla' = 'gor', 'human' = 'human1',
                              'bonobo' = 'bonobo'))



# #######################################################################
# 
# png("C:/Users/gresh/Downloads/temp/bonobo_cor.png", 
#     width = 500, height = 500, units = "px")
# 
# bonobo_cor <-
#   expression_isoform_count_copies %>% 
#   filter(species == 'bonobo') %>% 
#   select(c(structural, sequence, mean_exp_log, n_genes, n_coding)) %>% 
#   cor_mat() %>% 
#   cor_reorder() %>%
#   pull_lower_triangle() %>% 
#   cor_plot(method = "color")
# 
# dev.off()
# 
# #####################################################################
# 
# png("C:/Users/gresh/Downloads/temp/chimp_cor.png", 
#     width = 500, height = 500, units = "px")
# 
# expression_isoform_count_copies %>% 
#   filter(species == 'chimp') %>% 
#   select(c(structural, sequence, mean_exp_log, n_genes, n_coding)) %>% 
#   cor_mat() %>% 
#   cor_reorder() %>%
#   pull_lower_triangle() %>% 
#   cor_plot(method = "color")
# 
# dev.off()
# 
# #####################################################################
# 
# png("C:/Users/gresh/Downloads/temp/sum_orang_cor.png", 
#     width = 500, height = 500, units = "px")
# 
# expression_isoform_count_copies %>% 
#   filter(species == 'sum_orang') %>% 
#   select(c(structural, sequence, mean_exp_log, n_genes, n_coding)) %>% 
#   cor_mat() %>% 
#   cor_reorder() %>%
#   pull_lower_triangle() %>% 
#   cor_plot(method = "color")
# 
# dev.off()
# 
# #####################################################################
# 
# png("C:/Users/gresh/Downloads/temp/bor_orang_cor.png", 
#     width = 500, height = 500, units = "px")
# 
# expression_isoform_count_copies %>% 
#   filter(species == 'bor_orang') %>% 
#   select(c(structural, sequence, mean_exp_log, n_genes, n_coding)) %>% 
#   cor_mat() %>% 
#   cor_reorder() %>%
#   pull_lower_triangle() %>% 
#   cor_plot(method = "color")
# 
# dev.off()
# 
# #####################################################################
# 
# png("C:/Users/gresh/Downloads/temp/gor_cor.png", 
#     width = 500, height = 500, units = "px")
# 
# expression_isoform_count_copies %>% 
#   filter(species == 'gor') %>% 
#   select(c(structural, sequence, mean_exp_log, n_genes, n_coding)) %>% 
#   cor_mat() %>% 
#   cor_reorder() %>%
#   pull_lower_triangle() %>% 
#   cor_plot(method = "color")
# 
# dev.off()
# 
# #####################################################################
# 
# png("C:/Users/gresh/Downloads/temp/human1_cor.png", 
#     width = 500, height = 500, units = "px")
# 
# expression_isoform_count_copies %>% 
#   filter(species == 'human1') %>% 
#   select(c(structural, sequence, mean_exp_log, n_genes, n_coding)) %>% 
#   cor_mat() %>% 
#   cor_reorder() %>%
#   pull_lower_triangle() %>% 
#   cor_plot(method = "color")
# 
# dev.off()
# 
# #####################################################################
# 
# ########################################################################
#isoform diversity

expression_isoform_count_copies <-
  merge(expression_isoform_count, gene_copies)

expression_isoform_count_copies <-
  expression_isoform_count_copies %>% 
  filter(mean_exp_log > 0)

expression_isoform_count_copies <-
  expression_isoform_count_copies %>% 
  filter(!gene_fam %in% c('BPY2', 'VCY', 'PRY', 'HSFY'))

expression_isoform_count_copies <-
  expression_isoform_count_copies %>% 
  mutate(n_genes_log = log(n_genes),
         structural_log = log(structural),
         sequence_log = log(sequence),
         n_sequence_log = log(n_sequence),
         n_structural_log = log(n_structural))

library(ggpubr)

sequence_n_genes <- 
ggscatter(expression_isoform_count_copies, 
          x = "sequence_log", 
          y = "n_genes_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = 0.25, label.y = 4)+
  xlab('Sequence isoform diversity, log-transformed')+
  ylab('Gene copies count, log-transformed')+
  labs(color = "species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                               "bonobo" = "#AA4499",
                               "chimpanzee" = "#F0E442",
                               "gorilla" = "#0072B2",
                               "Bornean orangutan" = "#D55E00",
                               "Sumatran orangutan" = "#888888"))

structural_n_genes <- 
  ggscatter(expression_isoform_count_copies, 
          x = "structural_log", 
          y = "n_genes_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = 0.25, label.y = 4)+
  xlab('Structural isoform diversity, log-transformed')+
  ylab('Gene copies count, log-transformed')+
  labs(color = "Species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                                "bonobo" = "#AA4499",
                                "chimpanzee" = "#F0E442",
                                "gorilla" = "#0072B2",
                                "Bornean orangutan" = "#D55E00",
                                "Sumatran orangutan" = "#888888"))

sequence_mean_expression <- 
ggscatter(expression_isoform_count_copies, 
          x = "sequence_log", 
          y = "mean_exp_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = 0.25, label.y = 7.5)+
  xlab('Sequence isoform diversity, log-transformed')+
  ylab('Mean expression, log-transformed, log-transformed')+
  labs(color = "Species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                                "bonobo" = "#AA4499",
                                "chimpanzee" = "#F0E442",
                                "gorilla" = "#0072B2",
                                "Bornean orangutan" = "#D55E00",
                                "Sumatran orangutan" = "#888888"))

structural_mean_expression <-
  ggscatter(expression_isoform_count_copies, 
          x = "structural_log", 
          y = "mean_exp_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = .25, label.y = 8)+
  xlab('Structural isoform diversity, log-transformed')+
  ylab('Mean expression, log-transformed')+
  labs(color = "Species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                                "bonobo" = "#AA4499",
                                "chimpanzee" = "#F0E442",
                                "gorilla" = "#0072B2",
                                "Bornean orangutan" = "#D55E00",
                                "Sumatran orangutan" = "#888888"))

########################################################################
#isoform count

n_sequence_n_genes <-
ggscatter(expression_isoform_count_copies, 
          x = "n_sequence_log", 
          y = "n_genes_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = 0.25, label.y = 5)+
  xlab('Sequence isoform count, log-transformed')+
  ylab('Gene copies count, log-transformed')+
  labs(color = "Species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                                "bonobo" = "#AA4499",
                                "chimpanzee" = "#F0E442",
                                "gorilla" = "#0072B2",
                                "Bornean orangutan" = "#D55E00",
                                "Sumatran orangutan" = "#888888"))

n_structural_n_genes <-
ggscatter(expression_isoform_count_copies, 
          x = "n_structural_log", 
          y = "n_genes_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = .25, label.y = 4)+
  xlab('Structural isoform count, log-transformed')+
  ylab('Gene copies count, log-transformed')+
  labs(color = "Species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                                "bonobo" = "#AA4499",
                                "chimpanzee" = "#F0E442",
                                "gorilla" = "#0072B2",
                                "Bornean orangutan" = "#D55E00",
                                "Sumatran orangutan" = "#888888"))

n_sequence_mean_expression <-
  ggscatter(expression_isoform_count_copies, 
          x = "n_sequence_log", 
          y = "mean_exp_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = .25, label.y = 8)+
  xlab('Sequence isoform count, log-transformed')+
  ylab('Mean expression, log-transformed')+
  labs(color = "Species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                                "bonobo" = "#AA4499",
                                "chimpanzee" = "#F0E442",
                                "gorilla" = "#0072B2",
                                "Bornean orangutan" = "#D55E00",
                                "Sumatran orangutan" = "#888888"))

n_structural_mean_expression <- 
  ggscatter(expression_isoform_count_copies, 
          x = "n_structural_log", 
          y = "mean_exp_log",
          color = "species",
          add = "reg.line", facet.by = "gene_fam",                               
          conf.int = FALSE,                                 
          add.params = list(color = "black",
                            fill = "lightgray",
                            size=0.5))+
  stat_cor(method = "pearson",
           label.x = .25, label.y = 8)+
  xlab('Structural isoform count, log-transformed')+
  ylab('Mean expression, log-transformed')+
  labs(color = "Species")+
  theme(strip.text = element_text(face = "bold.italic"))+
  scale_color_manual(values = c("human" = "#009E73",
                                "bonobo" = "#AA4499",
                                "chimpanzee" = "#F0E442",
                                "gorilla" = "#0072B2",
                                "Bornean orangutan" = "#D55E00",
                                "Sumatran orangutan" = "#888888"))


ggsave("C:/Users/gresh/Downloads/temp/n_sequence.png",
       ggarrange(n_sequence_mean_expression + rremove('xlab') + rremove('legend'), n_sequence_n_genes + rremove('legend'), 
                 ncol=1, nrow=2, labels=c('A', 'B')), 
       width=8, height=10, dpi=500)

ggsave("C:/Users/gresh/Downloads/temp/sequence.png",
       ggarrange(sequence_mean_expression + rremove('xlab') + rremove('legend'), sequence_n_genes + rremove('legend'), 
                 ncol=1, nrow=2, labels=c('A', 'B')), 
       width=8, height=10, dpi=500)

ggsave("C:/Users/gresh/Downloads/temp/n_structural.png",
       ggarrange(n_structural_mean_expression + rremove('xlab') + rremove('legend'), n_structural_n_genes + rremove('legend'), 
                 ncol=1, nrow=2, labels=c('A', 'B')), 
       width=8, height=10, dpi=500)

ggsave("C:/Users/gresh/Downloads/temp/structural.png",
       ggarrange(structural_mean_expression + rremove('xlab') + rremove('legend'), structural_n_genes + rremove('legend'), 
                 ncol=1, nrow=2, labels=c('A', 'B')), 
       width=8, height=10, dpi=500)

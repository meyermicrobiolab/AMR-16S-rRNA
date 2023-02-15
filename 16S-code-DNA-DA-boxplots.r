library(corncob)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(vegan)
library(reshape2)
library(cowplot)
library(grid)
library(scales)
library(RColorBrewer)
library(randomcoloR)
library(microViz)

otu <- read.table("psmin5-ASV-table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("taxa-table-renamed.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadataDNA2.txt",sep="\t",header=TRUE,row.names=1)

OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))



----------------------------------------------------
#Collapse table by Phylum 


ps_phylum <- ps %>% tax_glom("Phylum")

ps_phylum_ra <- transform_sample_counts(ps_phylum, function(otu) otu/sum(otu)*100)

ps_phylum_ra_subset <- subset_taxa(ps_phylum_ra, rownames(tax_table(ps_phylum)) %in%  c("ASV1","ASV2","ASV5","ASV8","ASV15","ASV28","ASV30","ASV55","ASV101","ASV154"))
phyloseq::taxa_names(ps_phylum_ra_subset) <- phyloseq::tax_table(ps_phylum_ra_subset)[, "Phylum"]

colors = c("darkseagreen4","cadetblue3","coral2")


phyloseq::psmelt(ps_phylum_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy","antibiotic","diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 14)) +   
  theme(axis.text = element_text(size = 12)) 
  

ggsave("DNA-sign-DA-phyla.tiff", units="in", width=10, height=10, dpi=300, compression = 'lzw')

-----------------------------------------------------------------------------------------------------

ps_class <- ps %>% tax_glom("Class")

ps_class_ra <- transform_sample_counts(ps_class, function(otu) otu/sum(otu)*100)

ps_class_ra_subset <- subset_taxa(ps_class_ra, rownames(tax_table(ps_class)) %in% c("ASV1","ASV2","ASV4","ASV15","ASV28","ASV30","ASV47","ASV98","ASV111","ASV135","ASV154"))

phyloseq::taxa_names(ps_class_ra_subset) <- phyloseq::tax_table(ps_class_ra_subset)[, "Class"]

colors = c("darkseagreen4","cadetblue3","coral2")


phyloseq::psmelt(ps_class_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy","antibiotic","diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 12)) +   
  theme(axis.text = element_text(size = 10))
  
 
ggsave("DNA-sign-DA-classes.tiff", units="in", width=10, height=10, dpi=300, compression = 'lzw')


--------------------------------------------------------------------------------------------------------------

ps_order <- ps %>% tax_glom("Order")

ps_order_ra <- transform_sample_counts(ps_order, function(otu) otu/sum(otu)*100)

ps_order_ra_subset <- subset_taxa(ps_order_ra, rownames(tax_table(ps_order)) %in% c("ASV2","ASV4","ASV13","ASV15","ASV18","ASV21","ASV27","ASV30","ASV50","ASV86","ASV91","ASV98","ASV101","ASV110","ASV135","ASV154","ASV167"))


phyloseq::taxa_names(ps_order_ra_subset) <- phyloseq::tax_table(ps_order_ra_subset)[, "Order"]

colors = c("darkseagreen4","cadetblue3","coral2")

phyloseq::psmelt(ps_order_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy", "antibiotic", "diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 12)) +   
  theme(axis.text = element_text(size = 10)) 

ggsave("DNA-sign-DA-orders.tiff", units="in", width=10, height=12, dpi=300, compression = 'lzw')





-----------------------------------------------------------------------------------------------------------------

ps_family <- ps %>% tax_glom("Family")

ps_family_ra <- transform_sample_counts(ps_family, function(otu) otu/sum(otu)*100)

ps_family_ra_subset <- subset_taxa(ps_family_ra, rownames(tax_table(ps_family_ra)) %in% c("ASV2","ASV4","ASV6","ASV13","ASV15","ASV17","ASV21","ASV22","ASV25","ASV27","ASV30","ASV34","ASV39","ASV50","ASV52","ASV60","ASV62","ASV70","ASV86","ASV91","ASV98","ASV101","ASV110","ASV116","ASV135","ASV154","ASV167"))

phyloseq::taxa_names(ps_family_ra_subset) <- phyloseq::tax_table(ps_family_ra_subset)[, "Family"]

colors = c("darkseagreen4","cadetblue3","coral2")

phyloseq::psmelt(ps_family_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy","antibiotic","diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 12)) +   
  theme(axis.text = element_text(size = 10)) 

ggsave("DNA-sign-DA-families.tiff", units="in", width=10, height=12, dpi=300, compression = 'lzw')


phyloseq::tax_table(ps_family)[, "Family"]


------------------------------------------------------------------------

ps_genus <- ps %>% tax_glom("Genus")

ps_genus_ra <- transform_sample_counts(ps_genus, function(otu) otu/sum(otu)*100)


ps_genus_ra_subset <- subset_taxa(ps_genus_ra, rownames(tax_table(ps_genus_ra)) %in% c("ASV2","ASV4","ASV6","ASV15","ASV17","ASV21","ASV25","ASV27","ASV30","ASV34","ASV41","ASV50","ASV52","ASV57","ASV60","ASV62","ASV70","ASV82","ASV86","ASV90","ASV91","ASV98","ASV101","ASV110","ASV116","ASV135","ASV154","ASV167","ASV207"))

phyloseq::taxa_names(ps_genus_ra_subset) <- phyloseq::tax_table(ps_genus_ra_subset)[, "Genus"]

colors = c("darkseagreen4","cadetblue3","coral2")

phyloseq::psmelt(ps_genus_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy","antibiotic","diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 12)) +   
  theme(axis.text = element_text(size = 10)) 
  
  
ggsave("sign-DA-genera-by-treatment-RNA.tiff", units="in", width=10, height=16, dpi=300, compression = 'lzw')


-----------------------------------------------------------------------------------------------------------
#only enriched in diseased 

ps_family <- ps %>% tax_glom("Family")

ps_family_ra <- transform_sample_counts(ps_family, function(otu) otu/sum(otu)*100)

ps_family_ra_subset <- subset_taxa(ps_family_ra, rownames(tax_table(ps_family_ra)) %in% c("ASV2","ASV13","ASV22","ASV60","ASV110"))

phyloseq::taxa_names(ps_family_ra_subset) <- phyloseq::tax_table(ps_family_ra_subset)[, "Family"]

colors = c("darkseagreen4","cadetblue3","coral2")

phyloseq::psmelt(ps_family_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy","antibiotic","diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 12)) +   
  theme(axis.text = element_text(size = 10)) 

ggsave("DNA-sign-enriched-in-diseased.tiff", units="in", width=8, height=5, dpi=300, compression = 'lzw')

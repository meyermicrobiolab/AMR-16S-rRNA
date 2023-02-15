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
samples <-read.table("metadataRNA2.txt",sep="\t",header=TRUE,row.names=1)

OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))

----------------------------------------------------
#Collapse table by Phylum 

ps_phylum <- ps %>% tax_glom("Phylum")

ntaxa(ps_phylum)
#21
nsamples(ps_phylum)
#29

sample_names(ps_phylum)[1:10]

sample_variables(ps_phylum)

otu_table(ps_phylum)[1:10, 1:10]

ps_phylum_ra <- transform_sample_counts(ps_phylum, function(otu) otu/sum(otu)*100)

#[1] "ASV1"  "ASV8"  "ASV11" "ASV15" "ASV30"


ps_phylum_ra_subset <- subset_taxa(ps_phylum_ra, rownames(tax_table(ps_phylum)) %in%  c("ASV1","ASV8","ASV11","ASV15","ASV30"))


otu_table(ps_phylum_ra_subset)[1:5, 1:5]


phyloseq::taxa_names(ps_phylum_ra_subset) <- phyloseq::tax_table(ps_phylum_ra_subset)[, "Phylum"]


otu_table(ps_phylum_ra_subset)[1:5, 1:5]

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
  

ggsave("RNA-sign-DA-phyla.tiff", units="in", width=10, height=6, dpi=300, compression = 'lzw')


--------------------------------------------------------------------------------------------------------------

ps_class <- ps %>% tax_glom("Class")

otu_table(ps_class)[1:10, 1:10]

ps_class_ra <- transform_sample_counts(ps_class, function(otu) otu/sum(otu)*100)

otu_table(ps_class_ra)[1:10, 1:10]

ps_class_ra_subset <- subset_taxa(ps_class_ra, rownames(tax_table(ps_class)) %in% c("ASV1","ASV11","ASV15","ASV16","ASV30"))


phyloseq::taxa_names(ps_class_ra_subset) <- phyloseq::tax_table(ps_class_ra_subset)[, "Class"]

otu_table(ps_class_ra_subset)[1:5, 1:5]


colors = c("darkseagreen4","cadetblue3","coral2")

phyloseq::psmelt(ps_class_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy","antibiotic","diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 14)) +   
  theme(axis.text = element_text(size = 10)) 
  

ggsave("RNA-sign-DA-classes.tiff", units="in", width=10, height=6, dpi=300, compression = 'lzw')


--------------------------------------------------------------------------------------------------------------


ps_order <- ps %>% tax_glom("Order")

ps_order_ra <- transform_sample_counts(ps_order, function(otu) otu/sum(otu)*100)

ps_order_ra_subset <- subset_taxa(ps_order_ra, rownames(tax_table(ps_order)) %in% c("ASV1","ASV3","ASV6","ASV11","ASV15","ASV16","ASV24","ASV30","ASV41","ASV50","ASV86","ASV191"))

phyloseq::taxa_names(ps_order_ra_subset) <- phyloseq::tax_table(ps_order_ra_subset)[, "Order"]

colors = c("darkseagreen4","cadetblue3","coral2")

phyloseq::psmelt(ps_order_ra_subset) %>%
mutate(treatment = factor(treatment, levels=c("healthy","antibiotic","diseased"))) %>%
ggplot(data = ., aes(x = treatment, y = Abundance)) +
  geom_boxplot(aes(color = treatment),outlier.shape  = NA ) +
  geom_jitter(aes(color = treatment), height = 0, width = .2) +
  labs(x = "", y = "Relative Abundance (%)") +
  facet_wrap(~ OTU, ncol=3, scales = "free") +
  scale_color_manual(values=colors) + 
  theme(text = element_text(size = 12)) +   
  theme(axis.text = element_text(size = 10)) 
  

ggsave("RNA-sign-DA-orders3by4-.tiff", units="in", width=9, height=8, dpi=300, compression = 'lzw')


-----------------------------------------------------------------------------------------------------------------
#boxplot on diff abund families 

ps_family <- ps %>% tax_glom("Family")

ps_family_ra <- transform_sample_counts(ps_family, function(otu) otu/sum(otu)*100)

ps_family_ra_subset <- subset_taxa(ps_family_ra, rownames(tax_table(ps_family_ra)) %in% c("ASV1","ASV3","ASV6","ASV11","ASV15","ASV16","ASV30","ASV33","ASV41","ASV50","ASV52","ASV62","ASV70","ASV86","ASV191"))

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
  
  
ggsave("RNA-sign-DA-families.tiff", units="in", width=9, height=10, dpi=300, compression = 'lzw')


phyloseq::tax_table(ps_family)[, "Family"]


------------------------------------------------------------------------
ps_genus <- ps %>% tax_glom("Genus")

ps_genus_ra <- transform_sample_counts(ps_genus, function(otu) otu/sum(otu)*100)


ps_genus_ra_subset <- subset_taxa(ps_genus_ra, rownames(tax_table(ps_genus_ra)) %in% c("ASV1","ASV3","ASV6","ASV11","ASV15","ASV16","ASV30","ASV33","ASV50","ASV52","ASV62","ASV69","ASV70","ASV86","ASV99","ASV188","ASV191"))

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
  
ggsave("RNA-sign-DA-genera.tiff", units="in", width=10, height=12, dpi=300, compression = 'lzw')



-----------------------------------------------------------------------------------------------------
#only enriched in diseased 

ps_family <- ps %>% tax_glom("Family")

ps_family_ra <- transform_sample_counts(ps_family, function(otu) otu/sum(otu)*100)

ps_family_ra_subset <- subset_taxa(ps_family_ra, rownames(tax_table(ps_family_ra)) %in% c("ASV1","ASV3","ASV41"))

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

ggsave("RNA-sign-enriched-in-diseased.tiff", units="in", width=8, height=3, dpi=300, compression = 'lzw')


----------------------------------------------
#only enriched after antibiotic 

ps_family <- ps %>% tax_glom("Family")

ps_family_ra <- transform_sample_counts(ps_family, function(otu) otu/sum(otu)*100)

ps_family_ra_subset <- subset_taxa(ps_family_ra, rownames(tax_table(ps_family_ra)) %in% c("ASV15","ASV52"))

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

ggsave("RNA-sign-enriched-antibiotic.tiff", units="in", width=8, height=3, dpi=300, compression = 'lzw')

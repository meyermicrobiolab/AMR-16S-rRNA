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


https://www.yanh.org/2021/01/01/microbiome-r/#abundance-bar-plot

otu <- read.table("psmin5-ASV-table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("taxa-table-renamed.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadataRNA2.txt",sep="\t",header=TRUE,row.names=1)

OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
               
sample_data(ps)$treatment<-factor(sample_data(ps)$treatment,levels=c("healthy","antibiotic","diseased"))


ps_phylum <- ps %>% tax_glom("Phylum")

ps_phylum_ra <- transform_sample_counts(ps_phylum, function(otu) otu/sum(otu)*100)


#change the data format for ggplot plotting 

ps.melt <- psmelt(ps_phylum_ra)

# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%
  group_by(treatment, Phylum) %>%
  mutate(median=median(Abundance))
  
  
# select group median > 1
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"

#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,treatment,Phylum) %>%
  summarise(Abundance=sum(Abundance))

colors = c("antiquewhite3","coral1","darkolivegreen2","cornflowerblue","darkgoldenrod1","aquamarine3","darkslategray","darkslateblue")



ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="Samples", y= "Relative Abundance %") +
  facet_wrap(~treatment, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(),axis.text.x.bottom = element_text(angle = -90))+
  scale_fill_manual(values=colors)+
  theme(text = element_text(size = 14)) +   
  theme(axis.text = element_text(size = 12)) 
  
ggsave("top-phyla-bar-graph-by-treatment-RNA-v3.tiff", units="in", width=8, height=4, dpi=300, compression = 'lzw')

-------------------------------------------------------------------------------------------------------------------

ps_class <- ps %>% tax_glom("Class")

ps_class_ra <- transform_sample_counts(ps_class, function(otu) otu/sum(otu)*100)

ps.melt <- psmelt(ps_class_ra)

# change to character for easy-adjusted level
ps.melt$Class <- as.character(ps.melt$Class)

ps.melt <- ps.melt %>%
  group_by(treatment, Class) %>%
  mutate(median=median(Abundance))
  
  
# select group median > 1
keep <- unique(ps.melt$Class[ps.melt$median > 1])
ps.melt$Class[!(ps.melt$Class %in% keep)] <- "< 1%"

#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,treatment,Class) %>%
  summarise(Abundance=sum(Abundance))

        
colors = c("antiquewhite3","coral1","darkolivegreen2","cornflowerblue","darkgoldenrod1","darkseagreen1",
        "cyan","aquamarine3","darkmagenta","darkslategray","darkslateblue")


ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill=Class)) + 
  labs(x="", y="Relative Abundance (%)") +
  facet_wrap(~treatment, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(),axis.text.x.bottom = element_text(angle = -90))+
  scale_fill_manual(values=colors) +    
  theme(text = element_text(size = 14)) +   
  theme(axis.text = element_text(size = 12)) 
  
  
ggsave("top-classes-bar-graph-by-treatment-RNA-v3.tiff", units="in", width=8, height=4, dpi=300, compression = 'lzw')

------------------------------------------------------------------------------------------------------------------------

ps_order <- ps %>% tax_glom("Order")

ps_order_ra <- transform_sample_counts(ps_order, function(otu) otu/sum(otu)*100)

ps.melt <- psmelt(ps_order_ra)

# change to character for easy-adjusted level
ps.melt$Order <- as.character(ps.melt$Order)

ps.melt <- ps.melt %>%
  group_by(treatment, Order) %>%
  mutate(median=median(Abundance))
  
  
# select group median > 1
keep <- unique(ps.melt$Order[ps.melt$median > 2])
ps.melt$Order[!(ps.melt$Order %in% keep)] <- "< 2%"

#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,treatment,Order) %>%
  summarise(Abundance=sum(Abundance))


#trying to match DNA-graph colors (<2%, Actinomarinales, SAR11, Blastocatellales, Cyanobacteriales, Cytophygales, Enterbacerales, 
# Flavobacteriales, Proteobacteria unclass, Pseudomonadales, Punicceispirillales, Rhodobacterilaes, Rhodospirillase, Sphinobacteriales, Synechoco) 

colors = c("antiquewhite3","coral1","darkolivegreen2","cornflowerblue","darkgoldenrod1","darkseagreen1","lightpink1",
    "cyan","brown1","bisque","aquamarine3","darkmagenta","darkslategray","lightgoldenrod1","darkslateblue")

#there are 17 orders in RNA (<2%, Actinomarinales, SAR11, Blastocatellales, Burkholderiales, Campylobacteriles, Cyanobacteriales, Cytophygales, Enterbacerales, 
# Flavobacteriales, "Lentisphareia", Proteobacteria unclass, Pseudomonadales, Punicceispirillales, Rhodobacterillales, Sphinobacteriales, Synechoco) 



colors = c("antiquewhite3","coral1","darkolivegreen2","cornflowerblue","cadetblue4","darkorchid","darkgoldenrod1","darkseagreen1","lightpink1",
    "cyan","lightgray","brown1","bisque","aquamarine3","darkmagenta","limegreen","darkslateblue")


ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Order)) + 
  geom_bar(stat = "identity", aes(fill=Order)) + 
  labs(x="Sample", y="Relative Abundance (%)") +
  facet_wrap(~treatment, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90))+
  scale_fill_manual(values=colors) +    
  theme(text = element_text(size = 14)) +   
  theme(axis.text = element_text(size = 12)) 

ggsave("top-orders-bar-graph-by-mean-abundance-RNA-v3.tiff", units="in", width=9, height=5, dpi=300, compression = 'lzw')


----------------------------------------------------------------------------------------------------------------------------
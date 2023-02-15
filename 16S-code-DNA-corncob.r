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

otu <- read.table("final-ASV-DNA-table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("final-taxa-table.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadataDNA.txt",sep="\t",header=TRUE,row.names=1)


OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
                tax_table(taxon))
ntaxa(ps) 
length(get_taxa_unique(ps, "Phylum"))
length(get_taxa_unique(ps, "Class"))
length(get_taxa_unique(ps, "Order"))
length(get_taxa_unique(ps, "Family"))
length(get_taxa_unique(ps, "Genus"))


#unflitered ASVs=6209, phyla=50, classes=118, orders=282, families=482, genera=892


#filter low abudant ASVs 
    
psmin5 <-filter_taxa(ps, function(x) mean(x) >5, TRUE)

write.table(as(otu_table(psmin5), "matrix"),"psmin5-ASV-table.txt",sep="\t",col.names=NA)

ntaxa(psmin5) 
length(get_taxa_unique(psmin5, "Phylum"))
length(get_taxa_unique(psmin5, "Class"))
length(get_taxa_unique(psmin5, "Order"))
length(get_taxa_unique(psmin5, "Family"))
length(get_taxa_unique(psmin5, "Genus"))


#flitered ASVs=350, phyla=20, classes=32, orders=67, families=106, genera=171


#alternativly min 10 reads on average 

psmin10 <- filter_taxa(ps, function(x) mean(x) >10, TRUE)
ntaxa(psmin10) 
length(get_taxa_unique(psmin10, "Phylum")) 
length(get_taxa_unique(psmin10, "Class")) 
length(get_taxa_unique(psmin10, "Order")) 
length(get_taxa_unique(psmin10, "Family")) 
length(get_taxa_unique(psmin10, "Genus")) 
#flitered10 ASVs=212, phyla=14, classes=21, orders=49, families=80, genera=118


#convert to relative abudance and save for excel if needed 

ps_asv_ra <- transform_sample_counts(psmin5, function(otu) otu/sum(otu)*100)
write.table(as(otu_table(ps_asv_ra), "matrix"),"ps_asv_ra.txt",sep="\t",col.names=NA)


# corncob run on filtered by min mean reads 5 

#diff abundance test per treatment (note: one of the treatments has to be a control - 
#but the test picks alphabetically -- this is way the label is ahealthy )

set.seed(1)
da_analysis <- differentialTest(formula = ~ treatment,
	phi.formula = ~ treatment,
	formula_null = ~ 1,
	phi.formula_null = ~ treatment,
	test = "Wald", boot = FALSE,
	data = psmin5,
	fdr_cutoff = 0.05)

da_analysis$significant_taxa 

sign.da.asv=otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = psmin5)
write.table(sign.da.asv,"sign.da.asv-treatment.txt",sep="\t",col.names=NA)
 
 #[1] "ASV4"   "ASV6"   "ASV7"   "ASV15"  "ASV17"  "ASV21"  "ASV25"  "ASV27"  "ASV30"  "ASV34" 
#[11] "ASV38"  "ASV41"  "ASV50"  "ASV52"  "ASV57"  "ASV62"  "ASV68"  "ASV70"  "ASV71"  "ASV74" 
#[21] "ASV84"  "ASV86"  "ASV90"  "ASV91"  "ASV98"  "ASV105" "ASV115" "ASV118" "ASV148" "ASV152"
#[31] "ASV154" "ASV167" "ASV191" "ASV203" "ASV208" "ASV261" "ASV269" "ASV395"


plot(da_analysis)
ggsave("da_ASVs_analysis.tiff", units="in", width=18, height=8, dpi=300, compression = 'lzw')	


sign1.out <- da_analysis[["significant_models"]][[1]][["coefficients"]]
sign2.out <- da_analysis[["significant_models"]][[2]][["coefficients"]]
sign3.out <- da_analysis[["significant_models"]][[3]][["coefficients"]]
sign4.out <- da_analysis[["significant_models"]][[4]][["coefficients"]]
sign5.out <- da_analysis[["significant_models"]][[5]][["coefficients"]]
sign6.out <- da_analysis[["significant_models"]][[6]][["coefficients"]]
sign7.out <- da_analysis[["significant_models"]][[7]][["coefficients"]]
sign8.out <- da_analysis[["significant_models"]][[8]][["coefficients"]]
sign9.out <- da_analysis[["significant_models"]][[9]][["coefficients"]]
sign10.out <- da_analysis[["significant_models"]][[10]][["coefficients"]]
sign11.out <- da_analysis[["significant_models"]][[11]][["coefficients"]]
sign12.out <- da_analysis[["significant_models"]][[12]][["coefficients"]]
sign13.out <- da_analysis[["significant_models"]][[13]][["coefficients"]]
sign14.out <- da_analysis[["significant_models"]][[14]][["coefficients"]]
sign15.out <- da_analysis[["significant_models"]][[15]][["coefficients"]]
sign16.out <- da_analysis[["significant_models"]][[16]][["coefficients"]]
sign17.out <- da_analysis[["significant_models"]][[17]][["coefficients"]]
sign18.out <- da_analysis[["significant_models"]][[18]][["coefficients"]]
sign19.out <- da_analysis[["significant_models"]][[19]][["coefficients"]]
sign20.out <- da_analysis[["significant_models"]][[20]][["coefficients"]]
sign21.out <- da_analysis[["significant_models"]][[21]][["coefficients"]]
sign22.out <- da_analysis[["significant_models"]][[22]][["coefficients"]]
sign23.out <- da_analysis[["significant_models"]][[23]][["coefficients"]]
sign24.out <- da_analysis[["significant_models"]][[24]][["coefficients"]]
sign25.out <- da_analysis[["significant_models"]][[25]][["coefficients"]]
sign26.out <- da_analysis[["significant_models"]][[26]][["coefficients"]]
sign27.out <- da_analysis[["significant_models"]][[27]][["coefficients"]]
sign28.out <- da_analysis[["significant_models"]][[28]][["coefficients"]]
sign29.out <- da_analysis[["significant_models"]][[29]][["coefficients"]]
sign30.out <- da_analysis[["significant_models"]][[30]][["coefficients"]]
sign31.out <- da_analysis[["significant_models"]][[31]][["coefficients"]]
sign32.out <- da_analysis[["significant_models"]][[32]][["coefficients"]]
sign33.out <- da_analysis[["significant_models"]][[33]][["coefficients"]]
sign34.out <- da_analysis[["significant_models"]][[34]][["coefficients"]]
sign35.out <- da_analysis[["significant_models"]][[35]][["coefficients"]]
sign36.out <- da_analysis[["significant_models"]][[36]][["coefficients"]]
sign37.out <- da_analysis[["significant_models"]][[37]][["coefficients"]]
sign38.out <- da_analysis[["significant_models"]][[38]][["coefficients"]]


sign.ASVs.all.coeff <- rbind(sign1.out,sign2.out,sign3.out,sign4.out,sign5.out,sign6.out,sign7.out,sign8.out,sign9.out,sign10.out,sign11.out,sign12.out,sign13.out,sign14.out,sign15.out,sign16.out,sign17.out,sign18.out,sign19.out,sign20.out,sign21.out,sign22.out,sign23.out,sign24.out,sign25.out,sign26.out,sign27.out,sign28.out,sign29.out,sign30.out,sign31.out,sign32.out,sign33.out,sign34.out,sign35.out,sign36.out,sign37.out,sign38.out)

write.table(sign.ASVs.all.coeff,"sign.ASVs.coeff.treatment.txt",sep="\t",col.names=NA)

  


--------------------------------------
#collapse data by phylum 
set.seed(1)
ps_phylum <- psmin5 %>% tax_glom("Phylum")


ps_phylum_ra <- transform_sample_counts(ps_phylum, function(otu) otu/sum(otu)*100)

write.table(as(otu_table(ps_phylum_ra), "matrix"),"ps_phylum-ra-table.txt",sep="\t",col.names=NA)

da_phylum <- differentialTest(formula = ~ treatment,
	phi.formula = ~ treatment,
	formula_null = ~ 1,
	phi.formula_null = ~ treatment,
	test = "Wald", boot = FALSE,
	data = ps_phylum,
	fdr_cutoff = 0.05)
    

sign.phyla.treatment <- otu_to_taxonomy(OTU = da_phylum$significant_taxa, data = ps_phylum)
write.table(sign.phyla.treatment,"sign.phyla.treatment.txt",sep="\t",col.names=NA)

ASV1.out <- da_phylum[["significant_models"]][[1]][["coefficients"]]

da_phylum$significant_taxa 

#[1] "ASV1"   "ASV2"   "ASV5"   "ASV8"   "ASV15"  "ASV28"  "ASV30"  "ASV55"  "ASV101" "ASV154"


ps.phyla.sign <- subset_taxa(ps_phylum, rownames(tax_table(ps_phylum)) %in% c("ASV1","ASV2","ASV5","ASV8","ASV15","ASV28","ASV30","ASV55","ASV101","ASV154"))
write.table(as(otu_table(ps.phyla.sign), "matrix"),"sign-phyla.txt",sep="\t",col.names=NA)

plot(da_phylum)
ggsave("da_phylum.tiff", units="in", width=12, height=8, dpi=300, compression = 'lzw')

#the complete corncob outout is complicate but here is one way for getting it out 

sign1.out <- da_phylum[["significant_models"]][[1]][["coefficients"]]
sign2.out <- da_phylum[["significant_models"]][[2]][["coefficients"]]
sign3.out <- da_phylum[["significant_models"]][[3]][["coefficients"]]
sign4.out <- da_phylum[["significant_models"]][[4]][["coefficients"]]
sign5.out <- da_phylum[["significant_models"]][[5]][["coefficients"]]
sign6.out <- da_phylum[["significant_models"]][[6]][["coefficients"]]
sign7.out <- da_phylum[["significant_models"]][[7]][["coefficients"]]
sign8.out <- da_phylum[["significant_models"]][[8]][["coefficients"]]
sign9.out <- da_phylum[["significant_models"]][[9]][["coefficients"]]
sign10.out <- da_phylum[["significant_models"]][[10]][["coefficients"]]


sign.phyla.all.coeff <- rbind(sign1.out,sign2.out,sign3.out,sign4.out,sign5.out,sign6.out,sign7.out,sign8.out,sign9.out,sign10.out)

write.table(sign.phyla.all.coeff,"sign.phyla.coeff.treatment.txt",sep="\t",col.names=NA)



-------------------------------------
#to create class-level-table

ps_class <- psmin5 %>% tax_glom("Class")
ps_class_ra <- transform_sample_counts(ps_class, function(otu) otu/sum(otu)*100)
write.table(as(otu_table(ps_class_ra), "matrix"),"ps_class-level.txt",sep="\t",col.names=NA)

da_class <- differentialTest(formula = ~ treatment,
	phi.formula = ~ treatment,
	formula_null = ~ 1,
	phi.formula_null = ~ treatment,
	test = "Wald", boot = FALSE,
	data = ps_class,
	fdr_cutoff = 0.05)
    
da_class$significant_taxa 
#[1] "ASV1"   "ASV2"   "ASV4"   "ASV15"  "ASV28"  "ASV30"  "ASV47"  "ASV98"  "ASV101" "ASV135"
#[11] "ASV154"

plot(da_class)
ggsave("da_class.tiff", units="in", width=12, height=8, dpi=300, compression = 'lzw')


sign.classes.treatment <- otu_to_taxonomy(OTU = da_class$significant_taxa, data = ps_class)
write.table(sign.classes.treatment,"sign.classes.treatment.txt",sep="\t",col.names=NA)

#the complete corncob outout is complicate but here is one way for getting it out 

sign1.out <- da_class[["significant_models"]][[1]][["coefficients"]]
sign2.out <- da_class[["significant_models"]][[2]][["coefficients"]]
sign3.out <- da_class[["significant_models"]][[3]][["coefficients"]]
sign4.out <- da_class[["significant_models"]][[4]][["coefficients"]]
sign5.out <- da_class[["significant_models"]][[5]][["coefficients"]]
sign6.out <- da_class[["significant_models"]][[6]][["coefficients"]]
sign7.out <- da_class[["significant_models"]][[7]][["coefficients"]]
sign8.out <- da_class[["significant_models"]][[8]][["coefficients"]]
sign9.out <- da_class[["significant_models"]][[9]][["coefficients"]]
sign10.out <- da_class[["significant_models"]][[10]][["coefficients"]]
sign11.out <- da_class[["significant_models"]][[11]][["coefficients"]]

sign.classes.all.coeff <- rbind(sign1.out,sign2.out,sign3.out,sign4.out,sign5.out,sign6.out,sign7.out,sign8.out,sign9.out,sign10.out,sign11.out)

write.table(sign.classes.all.coeff,"sign.classes.coeff.treatment.txt",sep="\t",col.names=NA)

-----------------------------------------------------------------------------------------
#to create order-level-table

ps_order <- psmin5 %>% tax_glom("Order")
ps_order_ra <- transform_sample_counts(ps_order, function(otu) otu/sum(otu)*100)
write.table(as(otu_table(ps_order_ra), "matrix"),"ps_order-level.txt",sep="\t",col.names=NA)

da_order <- differentialTest(formula = ~ treatment,
	phi.formula = ~ treatment,
	formula_null = ~ 1,
	phi.formula_null = ~ treatment,
	test = "Wald", boot = FALSE,
	data = ps_order,
	fdr_cutoff = 0.05)
	
da_order$significant_taxa   
#[1] "ASV2"   "ASV4"   "ASV13"  "ASV15"  "ASV18"  "ASV21"  "ASV27"  "ASV30"  "ASV50"  "ASV86"  "ASV91"  "ASV98" 
#[13] "ASV101" "ASV110" "ASV135" "ASV154" "ASV167"

plot(da_order)
ggsave("da_order.tiff", units="in", width=15, height=8, dpi=300, compression = 'lzw')



sign.orders.treatment <- otu_to_taxonomy(OTU = da_order$significant_taxa, data = ps_order)
write.table(sign.orders.treatment,"sign.orders.treatment.txt",sep="\t",col.names=NA)
  
  
sign1.out <- da_order[["significant_models"]][[1]][["coefficients"]]
sign2.out <- da_order[["significant_models"]][[2]][["coefficients"]]
sign3.out <- da_order[["significant_models"]][[3]][["coefficients"]]
sign4.out <- da_order[["significant_models"]][[4]][["coefficients"]]
sign5.out <- da_order[["significant_models"]][[5]][["coefficients"]]
sign6.out <- da_order[["significant_models"]][[6]][["coefficients"]]
sign7.out <- da_order[["significant_models"]][[7]][["coefficients"]]
sign8.out <- da_order[["significant_models"]][[8]][["coefficients"]]
sign9.out <- da_order[["significant_models"]][[9]][["coefficients"]]
sign10.out <- da_order[["significant_models"]][[10]][["coefficients"]]
sign11.out <- da_order[["significant_models"]][[11]][["coefficients"]]
sign12.out <- da_order[["significant_models"]][[12]][["coefficients"]]
sign13.out <- da_order[["significant_models"]][[13]][["coefficients"]]
sign14.out <- da_order[["significant_models"]][[14]][["coefficients"]]
sign15.out <- da_order[["significant_models"]][[15]][["coefficients"]]
sign16.out <- da_order[["significant_models"]][[16]][["coefficients"]]
sign17.out <- da_order[["significant_models"]][[17]][["coefficients"]]


sign.orders.all.coeff <- rbind(sign1.out,sign2.out,sign3.out,sign4.out,sign5.out,sign6.out,sign7.out,sign8.out,sign9.out,sign10.out,sign11.out,sign12.out,sign13.out,sign14.out,sign15.out,sign16.out,sign17.out)

write.table(sign.orders.all.coeff,"sign.orders.coeff.treatment.txt",sep="\t",col.names=NA)

  
  
  
-------------------------------------------------------------------------------------------- 

#to create family-level-table 

ps_family <- psmin5 %>% tax_glom("Family")
ps_family_ra <- transform_sample_counts(ps_family, function(otu) otu/sum(otu)*100)
write.table(as(otu_table(ps_family_ra), "matrix"),"ps_family-table.txt",sep="\t",col.names=NA)

da_family <- differentialTest(formula = ~ treatment,
	phi.formula = ~ treatment,
	formula_null = ~ 1,
	phi.formula_null = ~ treatment,
	test = "Wald", boot = FALSE,
	data = ps_family,
	fdr_cutoff = 0.05)

da_family$significant_taxa   
#[1] "ASV2"   "ASV4"   "ASV6"   "ASV13"  "ASV15"  "ASV17"  "ASV21"  "ASV22"  "ASV25"  "ASV27"  "ASV30" 
#[12] "ASV34"  "ASV39"  "ASV50"  "ASV52"  "ASV60"  "ASV62"  "ASV70"  "ASV86"  "ASV91"  "ASV98"  "ASV101"
#[23] "ASV110" "ASV116" "ASV135" "ASV154" "ASV167"

plot(da_family)
ggsave("da_family.tiff", units="in", width=15, height=8, dpi=300, compression = 'lzw')	

sign.families.treatment <- otu_to_taxonomy(OTU = da_family$significant_taxa, data = ps_family)
write.table(sign.families.treatment,"sign.families.treatment.txt",sep="\t",col.names=NA)


sign1.out <- da_family[["significant_models"]][[1]][["coefficients"]]
sign2.out <- da_family[["significant_models"]][[2]][["coefficients"]]
sign3.out <- da_family[["significant_models"]][[3]][["coefficients"]]
sign4.out <- da_family[["significant_models"]][[4]][["coefficients"]]
sign5.out <- da_family[["significant_models"]][[5]][["coefficients"]]
sign6.out <- da_family[["significant_models"]][[6]][["coefficients"]]
sign7.out <- da_family[["significant_models"]][[7]][["coefficients"]]
sign8.out <- da_family[["significant_models"]][[8]][["coefficients"]]
sign9.out <- da_family[["significant_models"]][[9]][["coefficients"]]
sign10.out <- da_family[["significant_models"]][[10]][["coefficients"]]
sign11.out <- da_family[["significant_models"]][[11]][["coefficients"]]
sign12.out <- da_family[["significant_models"]][[12]][["coefficients"]]
sign13.out <- da_family[["significant_models"]][[13]][["coefficients"]]
sign14.out <- da_family[["significant_models"]][[14]][["coefficients"]]
sign15.out <- da_family[["significant_models"]][[15]][["coefficients"]]
sign16.out <- da_family[["significant_models"]][[16]][["coefficients"]]
sign17.out <- da_family[["significant_models"]][[17]][["coefficients"]]
sign18.out <- da_family[["significant_models"]][[18]][["coefficients"]]
sign19.out <- da_family[["significant_models"]][[19]][["coefficients"]]
sign20.out <- da_family[["significant_models"]][[20]][["coefficients"]]
sign21.out <- da_family[["significant_models"]][[21]][["coefficients"]]
sign22.out <- da_family[["significant_models"]][[22]][["coefficients"]]
sign23.out <- da_family[["significant_models"]][[23]][["coefficients"]]
sign24.out <- da_family[["significant_models"]][[24]][["coefficients"]]
sign25.out <- da_family[["significant_models"]][[25]][["coefficients"]]
sign26.out <- da_family[["significant_models"]][[26]][["coefficients"]]
sign27.out <- da_family[["significant_models"]][[27]][["coefficients"]]

sign.families.all.coeff <- rbind(sign1.out,sign2.out,sign3.out,sign4.out,sign5.out,sign6.out,sign7.out,sign8.out,sign9.out,sign10.out,sign11.out,sign12.out,sign13.out,sign14.out,sign15.out,sign16.out,sign17.out,sign18.out,sign19.out,sign20.out,sign21.out,sign22.out,sign23.out,sign24.out,sign25.out,sign26.out,sign27.out)

write.table(sign.families.all.coeff,"sign.families.coeff.treatment.txt",sep="\t",col.names=NA)

  
------------------------------------------------------------------------------

ps_genus <- psmin5 %>% tax_glom("Genus")
ps_genus_ra <- transform_sample_counts(ps_genus, function(otu) otu/sum(otu)*100)
write.table(as(otu_table(ps_genus_ra), "matrix"),"ps_genus-level.txt",sep="\t",col.names=NA)

da_genus <- differentialTest(formula = ~ treatment,
	phi.formula = ~ treatment,
	formula_null = ~ 1,
	phi.formula_null = ~ treatment,
	test = "Wald", boot = FALSE,
	data = ps_genus,
	fdr_cutoff = 0.05)
	
da_genus$significant_taxa 
#[1] "ASV2"   "ASV4"   "ASV6"   "ASV15"  "ASV17"  "ASV21"  "ASV25"  "ASV27"  "ASV30"  "ASV34"  "ASV41"  "ASV50"  "ASV52" 
#[14] "ASV57"  "ASV60"  "ASV62"  "ASV70"  "ASV82"  "ASV86"  "ASV90"  "ASV91"  "ASV98"  "ASV101" "ASV110" "ASV116" "ASV135"
#[27] "ASV154" "ASV167" "ASV207"

plot(da_genus)
ggsave("da_genus.tiff", units="in", width=15, height=8, dpi=300, compression = 'lzw')	


sign1.out <- da_genus[["significant_models"]][[1]][["coefficients"]]
sign2.out <- da_genus[["significant_models"]][[2]][["coefficients"]]
sign3.out <- da_genus[["significant_models"]][[3]][["coefficients"]]
sign4.out <- da_genus[["significant_models"]][[4]][["coefficients"]]
sign5.out <- da_genus[["significant_models"]][[5]][["coefficients"]]
sign6.out <- da_genus[["significant_models"]][[6]][["coefficients"]]
sign7.out <- da_genus[["significant_models"]][[7]][["coefficients"]]
sign8.out <- da_genus[["significant_models"]][[8]][["coefficients"]]
sign9.out <- da_genus[["significant_models"]][[9]][["coefficients"]]
sign10.out <- da_genus[["significant_models"]][[10]][["coefficients"]]
sign11.out <- da_genus[["significant_models"]][[11]][["coefficients"]]
sign12.out <- da_genus[["significant_models"]][[12]][["coefficients"]]
sign13.out <- da_genus[["significant_models"]][[13]][["coefficients"]]
sign14.out <- da_genus[["significant_models"]][[14]][["coefficients"]]
sign15.out <- da_genus[["significant_models"]][[15]][["coefficients"]]
sign16.out <- da_genus[["significant_models"]][[16]][["coefficients"]]
sign17.out <- da_genus[["significant_models"]][[17]][["coefficients"]]
sign18.out <- da_genus[["significant_models"]][[18]][["coefficients"]]
sign19.out <- da_genus[["significant_models"]][[19]][["coefficients"]]
sign20.out <- da_genus[["significant_models"]][[20]][["coefficients"]]
sign21.out <- da_genus[["significant_models"]][[21]][["coefficients"]]
sign22.out <- da_genus[["significant_models"]][[22]][["coefficients"]]
sign23.out <- da_genus[["significant_models"]][[23]][["coefficients"]]
sign24.out <- da_genus[["significant_models"]][[24]][["coefficients"]]
sign25.out <- da_genus[["significant_models"]][[25]][["coefficients"]]
sign26.out <- da_genus[["significant_models"]][[26]][["coefficients"]]
sign27.out <- da_genus[["significant_models"]][[27]][["coefficients"]]
sign28.out <- da_genus[["significant_models"]][[28]][["coefficients"]]
sign29.out <- da_genus[["significant_models"]][[29]][["coefficients"]]



sign.genera.all.coeff <- rbind(sign1.out,sign2.out,sign3.out,sign4.out,sign5.out,sign6.out,sign7.out,sign8.out,sign9.out,sign10.out,sign11.out,sign12.out,sign13.out,sign14.out,sign15.out,sign16.out,sign17.out,sign18.out,sign19.out,sign20.out,sign21.out,sign22.out,sign23.out,sign24.out,sign25.out,sign26.out,sign27.out,sign28.out,sign29.out)

write.table(sign.genera.all.coeff,"sign.genera.coeff.treatment.txt",sep="\t",col.names=NA)

 

sign.genera.treatment <- otu_to_taxonomy(OTU = da_genus$significant_taxa, data = ps_genus)
 
write.table(sign.genera.treatment,"sign.genera.treatment.txt",sep="\t",col.names=NA)



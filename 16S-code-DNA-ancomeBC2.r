
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ANCOMBC")
https://github.com/FrederickHuangLin/ANCOMBC


library(ANCOMBC)
library(microbiome)
library(phyloseq)
library(DT)

#note the otu/asv-table must be in format of samples names in first column and ASVs in first row (no the other way ! will not work) 

otu <- read.table("final-ASV-DNA-table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("final-taxa-table.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadataDNA.txt",sep="\t",header=TRUE,row.names=1)

OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))

ps <-filter_taxa(ps, function(x) mean(x) >5, TRUE)

               
#Visualize data

sample_names(ps)
sample_variables(ps)

#run AncomeBC2 


out = ancombc2(data = ps, assay_name = "counts", tax_level = "Phylum",
               fix_formula = "treatment",
               rand_formula = NULL,
               p_adj_method = "BH")
               
               
res = out$res

#the below will display resluts (beta, se, w, p_val, q_val, diff_abn) -- you can copy and paste or html tables =code below 
out$res 
#Result from the ANCOM-BC log-linear model to determine taxa that are differentially abundant 
#according to the covariate of interest. It contains: 1) coefficients; 2) standard errors; 
#3) test statistics; 4) p-values; 5) adjusted p-values; 
#6) indicators whether the taxon is differentially abundant (TRUE) or not (FALSE).

write.table(res,file="ancomBC-results-BH-phyla-by-treatment.txt")

--------------------------------------------------------------------------------------
ancome_class_out = ancombc2(data = ps, assay_name = "counts", tax_level = "Class",
               fix_formula = "treatment",
               rand_formula = NULL,
               p_adj_method = "holm")
               
res_class = ancome_class_out$res     

write.table(res_class,file="ancomBC-results-class_holm-by-treatment.txt")   



ancome_class_out = ancombc2(data = ps, assay_name = "counts", tax_level = "Class",
               fix_formula = "treatment",
               rand_formula = NULL,
               p_adj_method = "BH")
               
res_class = ancome_class_out$res     

write.table(res_class,file="ancomBC-results-class-BH-by-treatment.txt")   


----------------------------------------------------------------------------------------

ancome_order_out = ancombc2(data = ps, assay_name = "counts", tax_level = "Order",
               fix_formula = "treatment",
               rand_formula = NULL,
               p_adj_method = "BH")
               
res_order = ancome_order_out$res     

write.table(res_order,file="ancomBC-results-Order-BH-by-treatment.txt")   

-------------------------------------------------------------------------------------------
ancome_family_out = ancombc2(data = ps, assay_name = "counts", tax_level = "Family",
               fix_formula = "treatment",
               rand_formula = NULL,
               p_adj_method = "BH")
               
res_order = ancome_order_out$res     

write.table(res_order,file="ancomBC-results-Order-BH-by-treatment.txt")   




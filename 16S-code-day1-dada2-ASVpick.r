#installation of R packages 

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2",force=TRUE)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ALDEx2",force=TRUE)

install.packages("exactRankTests")
install.packages("ShortRead")
install.packages("zCompositions")
install.packages("car")
install.packages("propr")
install.packages("randomcoloR")


----------------------------------------------
library(dada2)
library(ShortRead)
library(ggplot2)
library(phyloseq)
library(vegan)
library(knitr)
library(ALDEx2)
library(zCompositions)
library(igraph)
library(car)
library(grDevices)
library(propr)
library(cowplot)
library(randomcoloR)
library(dplyr)
library(reshape2)
library(tibble)
library(exactRankTests)
library(nlme)
library(data.table)
library(scales)



###################### Quality-filter reads and create Amplicon Sequence Variant tables 

###adjust path name as appropriate
getwd()
 
path <- "D:/AMR_RNAseq_16S"

list.files(path)
#you have to see your fasta.gz files -- if you get character() it means dada2 cannot see your sequencing data 

# Samplename is everything before the first underscore
fnFs <- sort(list.files(path, pattern="_R1_cut.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_*.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Perform filtering and trimming
# Assign the filenames for the filtered fastq.gz files.
# Make directory and filenames for the filtered fastqs
# Place filtered files in filtered/ subdirectory
filt_path <- file.path(path, "filtered") 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter the forward and reverse reads
# WINDOWS USERS: set multithread=FALSE
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)



head(out)
dim(out)
write.table(out,file="filter-information")

# Learn the Error Rates, it TAKES TIME! do first forward and then reverse
# Forward reads
errF <- learnErrors(filtFs, multithread=FALSE)
# Reverse reads
errR <- learnErrors(filtRs, multithread=FALSE)

# visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


#save plot error files 
tiff(file="ploerrF.tiff",
width=6, height=4, units="in", res=300)
plotErrors(errF, nominalQ=TRUE)
dev.off()

tiff(file="ploterrR.tiff",
width=6, height=4, units="in", res=300)
plotErrors(errR, nominalQ=TRUE)
dev.off()



# Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

save(dadaFs, file = "dadaFs.RData")
save(dadaRs, file = "dadaRs.RData")

head(dadaFs)
dim(dadaFs)
head(dadaRs)
dim(dadaRs)


# Inspecting the dada-class object returned by dada:
dadaFs[[1]]

# Merge the denoised forward and reverse reads:

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
save(mergers, file = "mergers.RData")

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
write.table(seqtab,file="seqtab")

# Inspect distribution of sequence lengths
AVSseqlengthtab <- table(nchar(getSequences(seqtab)))
write.table(AVSseqlengthtab,file="AVSseqlengthtab")


#Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
write.table(seqtab.nochim,file="seqtab.nochim")
save(seqtab.nochim, file = "seqtab.nochim.RData")



# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, "dada_read_stats.txt",sep="\t",col.names=NA)

# SAVE THIS FILE SO YOU DON'T HAVE TO REPEAT ALL OF THE ABOVE STEPS, adjust name

saveRDS(seqtab.nochim, file="seqtab.nochim.rds")




# RELOAD THE SAVED INFO FROM HERE (if you have closed the project):

seqtab.nochim <- readRDS("seqtab.nochim.rds")

###################### ASSIGNING THE TAXONOMY

# Make sure the appropriate database is available in the DADA2 directory

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)

taxon <- as.data.frame(taxa,stringsAsFactors=FALSE)

# FIX the NAs in the taxa table

taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]


write.table(taxon,"silva_taxa_table.txt",sep="\t",col.names=NA)
write.table(seqtab.nochim, "silva_ASVs_table.txt",sep="\t",col.names=NA)



##################### REMOVING MITOCHONDIRAL AND CHLOROPLAST READS

# Create phyloseq object from otu and taxonomy tables from dada2, along with the sample metadata

otu <- read.table("silva_ASVs_table.txt",sep="\t",header=TRUE, row.names=1)

taxon <- read.table("silva_taxa_table.txt",sep="\t",header=TRUE,row.names=1)

samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)


OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))


get_taxa_unique(ps, "Family") 

get_taxa_unique(ps, "Order") 
get_taxa_unique(ps, "Kingdom") 
ps <- subset_taxa(ps, Family !="Mitochondria")
ps <- subset_taxa(ps, Order !="Chloroplast")
ps <- subset_taxa(ps, Kingdom !="Eukaryota")
ps <- subset_taxa(ps, Kingdom !="NA")
get_taxa_unique(ps, "Family") 
get_taxa_unique(ps, "Order") 
get_taxa_unique(ps, "Kingdom") 
ps 

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 6209 taxa and 52 samples ]
sample_data() Sample Data:       [ 52 samples by 3 sample variables ]
tax_table()   Taxonomy Table:    [ 6209 taxa by 6 taxonomic ranks ]



# filtered taxa with phyloseq, now export cleaned otu and taxa tables from phyloseq
otu = as(otu_table(ps), "matrix")
taxon = as(tax_table(ps), "matrix")
metadata = as(sample_data(ps), "matrix")
write.table(otu,"silva_nochloronomito_otu_table.txt",sep="\t",col.names=NA)
write.table(taxon,"silva_nochloronomito_taxa_table.txt",sep="\t",col.names=NA)


#you are done for the day


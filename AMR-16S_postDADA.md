---
title: "AMR-16S_postDADA"
author: "JM"
date: "2026-01-04"
output:
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
    number_sections: true
    theme: united
---



Load libraries


``` r
library(dada2)
```

```
## Loading required package: Rcpp
```

``` r
library(CoDaSeq)
```

```
## Loading required package: ALDEx2
```

```
## Loading required package: zCompositions
```

```
## Loading required package: MASS
```

```
## Loading required package: truncnorm
```

```
## Loading required package: survival
```

```
## Loading required package: lattice
```

```
## Loading required package: latticeExtra
```

```
## Loading required package: car
```

```
## Loading required package: carData
```

``` r
library(ggplot2)
```

```
## 
## Attaching package: 'ggplot2'
```

```
## The following object is masked from 'package:latticeExtra':
## 
##     layer
```

``` r
library(phyloseq)
library(vegan)
```

```
## Loading required package: permute
```

``` r
library(ANCOMBC)
library(corncob)
```

```
## 
## Attaching package: 'corncob'
```

```
## The following object is masked from 'package:car':
## 
##     logit
```

``` r
library(knitr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:car':
## 
##     recode
```

```
## The following object is masked from 'package:MASS':
## 
##     select
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

``` r
library(nloptr)
library(tibble)
library(randomcoloR)
library(cowplot)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

# Read in ASV table, taxonomy table, and metadata to create phyloseq object.


``` r
# load in data and create phyloseq object
otu <- read.table("silva_nochloronomito_otu_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps # 52 samples
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 6137 taxa and 52 samples ]
## sample_data() Sample Data:       [ 52 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 6137 taxa by 6 taxonomic ranks ]
```

``` r
psnb = subset_samples(ps, colony != "blank") #remove sample blanks
psnb # 50 samples
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 6137 taxa and 50 samples ]
## sample_data() Sample Data:       [ 50 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 6137 taxa by 6 taxonomic ranks ]
```

``` r
# remove disease samples to compare before and after amox treatment
ps1 = subset_samples(psnb, type != "disease") 
ps1 # 35 samples
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 6137 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 6137 taxa by 6 taxonomic ranks ]
```

``` r
otu = as(otu_table(ps1), "matrix")
taxon = as(tax_table(ps1), "matrix")
metadata = as(sample_data(ps1), "matrix")
write.table(otu,"silva_nochloronomito_otu_table_ps1.txt",sep="\t",col.names=NA)
write.table(taxon,"silva_nochloronomito_taxa_table_ps1.txt",sep="\t",col.names=NA)
write.table(metadata,"metadata_ps1.txt",sep="\t",col.names=NA)

ps5<-filter_taxa(ps1, function(x) mean(x) >5, TRUE)
ntaxa(ps5)
```

```
## [1] 112
```

``` r
ps10<-filter_taxa(ps1, function(x) mean(x) >10, TRUE)
ntaxa(ps10)
```

```
## [1] 78
```

``` r
get_taxa_unique(ps1, "Genus") 
```

```
##   [1] "Marinifilum"                                       
##   [2] "Vibrio"                                            
##   [3] "Clade III"                                         
##   [4] "Fusibacter"                                        
##   [5] "Clade Ia"                                          
##   [6] "P.palmC41"                                         
##   [7] "Algicola"                                          
##   [8] "Halarcobacter"                                     
##   [9] "Proteobacteria"                                    
##  [10] "Rhodobacteraceae"                                  
##  [11] "Thioflexothrix"                                    
##  [12] "Candidatus Actinomarina"                           
##  [13] "Sphingopyxis"                                      
##  [14] "Tenacibaculum"                                     
##  [15] "Amphritea"                                         
##  [16] "Thalassolituus"                                    
##  [17] "Ruegeria"                                          
##  [18] "Alteromonadaceae"                                  
##  [19] "Cryomorphaceae"                                    
##  [20] "Synechococcus CC9902"                              
##  [21] "Symphothece PCC-7002"                              
##  [22] "Litoricola"                                        
##  [23] "Hormoscilla SI04-45"                               
##  [24] "AEGEAN-169 marine group"                           
##  [25] "Delftia"                                           
##  [26] "Spirulina DRTO-55.2"                               
##  [27] "Blastocatellaceae"                                 
##  [28] "Pseudomonas"                                       
##  [29] "Clade II"                                          
##  [30] "Kordia"                                            
##  [31] "Lentisphaera"                                      
##  [32] "Cyanobacteriia"                                    
##  [33] "Fluviicola"                                        
##  [34] "HIMB11"                                            
##  [35] "Porticoccus"                                       
##  [36] "Phormidium ETS-05"                                 
##  [37] "Cribrihabitans"                                    
##  [38] "Roseibacillus"                                     
##  [39] "Rubritalea"                                        
##  [40] "SAR116 clade"                                      
##  [41] "Bacteroidetes VC2.1 Bac22"                         
##  [42] "SAR86 clade"                                       
##  [43] "Halodesulfovibrio"                                 
##  [44] "Microscilla"                                       
##  [45] "NS5 marine group"                                  
##  [46] "Aurantivirga"                                      
##  [47] "Cohaesibacter"                                     
##  [48] "Nautella"                                          
##  [49] "Campylobacterales"                                 
##  [50] "Thalassotalea"                                     
##  [51] "Spirulina P7"                                      
##  [52] "OM60(NOR5) clade"                                  
##  [53] "Enterobacterales"                                  
##  [54] "Pseudoteredinibacter"                              
##  [55] "Gammaproteobacteria"                               
##  [56] "Flavobacteriaceae"                                 
##  [57] "Oleiphilus"                                        
##  [58] "Oleibacter"                                        
##  [59] "S25-593"                                           
##  [60] "NS4 marine group"                                  
##  [61] "NS11-12 marine group"                              
##  [62] "Aquibacter"                                        
##  [63] "Aestuariibacter"                                   
##  [64] "Trichodesmium IMS101"                              
##  [65] "Urania-1B-19 marine sediment group"                
##  [66] "Sandaracinaceae"                                   
##  [67] "Candidatus Puniceispirillum"                       
##  [68] "MBIC10086"                                         
##  [69] "NS9 marine group"                                  
##  [70] "Bacteria"                                          
##  [71] "Cyclobacteriaceae"                                 
##  [72] "Bdellovibrio"                                      
##  [73] "Pseudomonadales"                                   
##  [74] "Hahella"                                           
##  [75] "MBMPE27"                                           
##  [76] "Candidatus Paraholospora"                          
##  [77] "Fulvivirga"                                        
##  [78] "Fabibacter"                                        
##  [79] "Agaribacter"                                       
##  [80] "Propionigenium"                                    
##  [81] "Blastocatella"                                     
##  [82] "Fokiniaceae"                                       
##  [83] "Pseudoruegeria"                                    
##  [84] "Balneola"                                          
##  [85] "KI89A clade"                                       
##  [86] "Synechococcus PCC-7336"                            
##  [87] "Profundimonas"                                     
##  [88] "Hyphomonas"                                        
##  [89] "Desulfocella"                                      
##  [90] "Pseudoalteromonas"                                 
##  [91] "Ferrimonas"                                        
##  [92] "Beggiatoaceae"                                     
##  [93] "Marinimicrobia (SAR406 clade)"                     
##  [94] "Actibacter"                                        
##  [95] "Cellvibrionaceae"                                  
##  [96] "Unknown Family"                                    
##  [97] "Alphaproteobacteria"                               
##  [98] "Rhodopirellula"                                    
##  [99] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
## [100] "Francisellaceae"                                   
## [101] "Woeseia"                                           
## [102] "Roseobacter clade CHAB-I-5 lineage"                
## [103] "PS1 clade"                                         
## [104] "Crocinitomicaceae"                                 
## [105] "Kordiimonas"                                       
## [106] "Puniceicoccaceae"                                  
## [107] "Bradymonadales"                                    
## [108] "Terasakiellaceae"                                  
## [109] "Staphylococcus"                                    
## [110] "Azotobacter"                                       
## [111] "Vibrionaceae"                                      
## [112] "Pedosphaeraceae"                                   
## [113] "Photobacterium"                                    
## [114] "Alteromonas"                                       
## [115] "Candidatus Amoebophilus"                           
## [116] "Paramoritella"                                     
## [117] "Pseudohaliea"                                      
## [118] "Coraliomargarita"                                  
## [119] "Saprospiraceae"                                    
## [120] "Serratia"                                          
## [121] "Desulfobacter"                                     
## [122] "Ekhidna"                                           
## [123] "Candidatus Cryptoprodotis"                         
## [124] "Rickettsiales"                                     
## [125] "Planktotalea"                                      
## [126] "Rubripirellula"                                    
## [127] "Thalassobaculales"                                 
## [128] "HTCC5015"                                          
## [129] "Comamonadaceae"                                    
## [130] "Aureitalea"                                        
## [131] "Shimia"                                            
## [132] "Acinetobacter"                                     
## [133] "JGI 0000069-P22"                                   
## [134] "Zeaxanthinibacter"                                 
## [135] "Reichenbachiella"                                  
## [136] "Xenococcus PCC-7305"                               
## [137] "Halioxenophilus"                                   
## [138] "Erythrobacter"                                     
## [139] "Desulfovibrio"                                     
## [140] "Arenicellaceae"                                    
## [141] "Allofrancisella"                                   
## [142] "Amoebophilaceae"                                   
## [143] "Limibaculum"                                       
## [144] "DEV008"                                            
## [145] "Cognatishimia"                                     
## [146] "Halomicronema TFEP1"                               
## [147] "Kiloniella"                                        
## [148] "Desulfocapsaceae"                                  
## [149] "Nitrincolaceae"                                    
## [150] "Denitrovibrio"                                     
## [151] "Marine Methylotrophic Group 3"                     
## [152] "Geminicoccaceae"                                   
## [153] "Arcobacteraceae"                                   
## [154] "Marine Group II"                                   
## [155] "Desulfosarcinaceae"                                
## [156] "Actinomarinales"                                   
## [157] "Neptuniibacter"                                    
## [158] "Dyella"                                            
## [159] "Phaeodactylibacter"                                
## [160] "Winogradskyella"                                   
## [161] "BD2-7"                                             
## [162] "Vicingus"                                          
## [163] "Saccharospirillaceae"                              
## [164] "HOC36"                                             
## [165] "IS-44"                                             
## [166] "R76-B128"                                          
## [167] "Labrenzia"                                         
## [168] "SCGC AAA164-E04"                                   
## [169] "Acrophormium PCC-7375"                             
## [170] "A4b"                                               
## [171] "Cyanobium PCC-6307"                                
## [172] "Portibacter"                                       
## [173] "NB1-j"                                             
## [174] "Chitinophagales"                                   
## [175] "Rhodospirillales"                                  
## [176] "Marinoscillum"                                     
## [177] "Gracilibacteria"                                   
## [178] "Seonamhaeicola"                                    
## [179] "Filomicrobium"                                     
## [180] "Desulfobotulus"                                    
## [181] "Phormidesmiaceae"                                  
## [182] "Lewinella"                                         
## [183] "Crocinitomix"                                      
## [184] "Kiloniellaceae"                                    
## [185] "Muricauda"                                         
## [186] "vadinHA49"                                         
## [187] "Flexithrix"                                        
## [188] "Hyphomonadaceae"                                   
## [189] "Nioella"                                           
## [190] "Lawsonella"                                        
## [191] "Aureicoccus"                                       
## [192] "Marinicella"                                       
## [193] "DEV007"                                            
## [194] "Halioglobus"                                       
## [195] "Aliikangiella"                                     
## [196] "Flammeovirga"                                      
## [197] "Corynebacterium"                                   
## [198] "Kapabacteriales"                                   
## [199] "Alcanivorax"                                       
## [200] "Endozoicomonas"                                    
## [201] "Chromatiales"                                      
## [202] "Micavibrionaceae"                                  
## [203] "[Caedibacter] taeniospiralis group"                
## [204] "Nannocystaceae"                                    
## [205] "NS10 marine group"                                 
## [206] "Thiohalorhabdaceae"                                
## [207] "Cryomorpha"                                        
## [208] "Halieaceae"                                        
## [209] "Arenicella"                                        
## [210] "OM190"                                             
## [211] "CI75cm.2.12"                                       
## [212] "Blastopirellula"                                   
## [213] "Phycisphaeraceae"                                  
## [214] "Tistlia"                                           
## [215] "Bacteroidia"                                       
## [216] "Limnothrix"                                        
## [217] "Phaeocystidibacter"                                
## [218] "SM2D12"                                            
## [219] "Gilvibacter"                                       
## [220] "OM182 clade"                                       
## [221] "SM1A02"                                            
## [222] "Leptospira"                                        
## [223] "Phycisphaera"                                      
## [224] "Persicobacter"                                     
## [225] "BD7-11"                                            
## [226] "Tropicimonas"                                      
## [227] "Desulfobulbus"                                     
## [228] "Rubidimonas"                                       
## [229] "Oceanospirillum"                                   
## [230] "OM27 clade"                                        
## [231] "Bacteroidales"                                     
## [232] "Stappiaceae"                                       
## [233] "Flavobacteriales"                                  
## [234] "Phormidium MBIC10003"                              
## [235] "Candidatus Megaira"                                
## [236] "Rothia"                                            
## [237] "Izemoplasmatales"                                  
## [238] "BD1-7 clade"                                       
## [239] "Parahaliea"                                        
## [240] "Absconditabacteriales (SR1)"                       
## [241] "Pirellula"                                         
## [242] "Gimesiaceae"                                       
## [243] "Flavirhabdus"                                      
## [244] "BD7-8"                                             
## [245] "Sumerlaea"                                         
## [246] "Leptospiraceae"                                    
## [247] "Oligoflexaceae"                                    
## [248] "Pla3 lineage"                                      
## [249] "Legionellaceae"                                    
## [250] "Sandaracinus"                                      
## [251] "Aquicella"                                         
## [252] "Algimonas"                                         
## [253] "Oscillatoriaceae"                                  
## [254] "Macellibacteroides"                                
## [255] "Pla4 lineage"                                      
## [256] "Micavibrionales"                                   
## [257] "Sva0996 marine group"                              
## [258] "BD2-11 terrestrial group"                          
## [259] "Halomonas"                                         
## [260] "Candidatus Navis"                                  
## [261] "Aureisphaera"                                      
## [262] "Enhydrobacter"                                     
## [263] "Flammeovirgaceae"                                  
## [264] "BIrii41"                                           
## [265] "Ardenticatenales"                                  
## [266] "PRD18C08"                                          
## [267] "Pseudenhygromyxa"                                  
## [268] "SU2 symbiont group"                                
## [269] "P3OB-42"                                           
## [270] "Sva0081 sediment group"                            
## [271] "Defluviicoccus"                                    
## [272] "Microcystaceae"                                    
## [273] "Acidovorax"                                        
## [274] "Carboxylicivirga"                                  
## [275] "Planktothricoides SR001"                           
## [276] "Cyanobacteriales"                                  
## [277] "Spongiibacteraceae"                                
## [278] "Leptolyngbyaceae"                                  
## [279] "possible genus 03"                                 
## [280] "Methyloceanibacter"                                
## [281] "Saprospira"                                        
## [282] "Chlamydiaceae"                                     
## [283] "Pirellulaceae"                                     
## [284] "Mycoplasma"                                        
## [285] "Roseibacterium"                                    
## [286] "Endothiovibrio"                                    
## [287] "Maritimibacter"                                    
## [288] "Haemophilus"                                       
## [289] "Cytophagales"                                      
## [290] "Ilumatobacter"                                     
## [291] "Anaerolineaceae UCG-001"                           
## [292] "Nannocystis"                                       
## [293] "Pir4 lineage"                                      
## [294] "Litorimicrobium"                                   
## [295] "Aquabacterium"                                     
## [296] "Myxococcota"                                       
## [297] "AT-s2-59"                                          
## [298] "Peredibacter"                                      
## [299] "Thiothrix"                                         
## [300] "Haloferula"                                        
## [301] "Hoeflea"                                           
## [302] "Lentimonas"                                        
## [303] "Reinekea"                                          
## [304] "Obscuribacteraceae"                                
## [305] "Verruc-01"                                         
## [306] "Planctomycetales"                                  
## [307] "Phormidesmiales"                                   
## [308] "Halovulum"                                         
## [309] "Congregibacter"                                    
## [310] "Cm1-21"                                            
## [311] "PB19"                                              
## [312] "Haliangium"                                        
## [313] "Pleurocapsa PCC-7319"                              
## [314] "Thermosynechococcales"                             
## [315] "Luminiphilus"                                      
## [316] "Oceanobacterium"                                   
## [317] "Subgroup 22"                                       
## [318] "UBA10353 marine group"                             
## [319] "Gven-F17"                                          
## [320] "Leptolyngbya PCC-6406"                             
## [321] "B2706-C7"                                          
## [322] "Terasakiella"                                      
## [323] "Chroococcidiopsis PCC-6712"                        
## [324] "Fuerstia"                                          
## [325] "Agarilytica"                                       
## [326] "Pseudohongiella"                                   
## [327] "Methylobacterium-Methylorubrum"                    
## [328] "Enhygromyxa"                                       
## [329] "Ascidiimonas"                                      
## [330] "Verrucomicrobiota"                                 
## [331] "Pseudobacteriovorax"                               
## [332] "Draconibacterium"                                  
## [333] "Spirochaeta 2"                                     
## [334] "Psychrosphaera"                                    
## [335] "Bacteriovoracaceae"                                
## [336] "FS142-36B-02"                                      
## [337] "SBR1031"                                           
## [338] "Clostridiaceae"                                    
## [339] "Persicirhabdus"                                    
## [340] "Nitrosococcaceae"                                  
## [341] "Candidatus Nitrosopumilus"                         
## [342] "SAR324 clade(Marine group B)"                      
## [343] "Sericytochromatia"                                 
## [344] "Romboutsia"                                        
## [345] "WPS-2"                                             
## [346] "Kordiimonadales"                                   
## [347] "Rhizobiales"                                       
## [348] "Agaribacterium"                                    
## [349] "Agarivorans"                                       
## [350] "Aureispira"                                        
## [351] "Cellulosilyticum"                                  
## [352] "Oceanicaulis"                                      
## [353] "Holosporaceae"                                     
## [354] "Thaumasiovibrio"                                   
## [355] "Edaphobaculum"                                     
## [356] "SH-PL14"                                           
## [357] "Algitalea"                                         
## [358] "Streptococcus"                                     
## [359] "Rubinisphaeraceae"                                 
## [360] "Latescibacterota"                                  
## [361] "Bacteroidetes BD2-2"                               
## [362] "C1-B045"                                           
## [363] "Rubritaleaceae"                                    
## [364] "Nodosilineaceae"                                   
## [365] "Anderseniella"                                     
## [366] "Geminobacterium"                                   
## [367] "Planctomycetota"                                   
## [368] "Cerasicoccus"                                      
## [369] "Chlamydiales"                                      
## [370] "UASB-TL25"                                         
## [371] "[Desulfobacterium] catecholicum group"             
## [372] "Pelagibius"                                        
## [373] "Cloacibacterium"                                   
## [374] "Parvularculaceae"                                  
## [375] "Salinirepens"                                      
## [376] "Shewanella"                                        
## [377] "NS7 marine group"                                  
## [378] "Phormidesmis ANT.LACV5.1"                          
## [379] "Caldilineaceae"                                    
## [380] "Microbacteriaceae"                                 
## [381] "Salinimonas"                                       
## [382] "Bermanella"                                        
## [383] "Bythopirellula"                                    
## [384] "Lentisphaeria"                                     
## [385] "Idiomarina"                                        
## [386] "Neorickettsia"                                     
## [387] "Neisseria"                                         
## [388] "Candidatus Tenderia"                               
## [389] "Granulosicoccales"                                 
## [390] "Sporolactobacillaceae"                             
## [391] "IheB3-7"                                           
## [392] "Cetobacterium"                                     
## [393] "Sphingobium"                                       
## [394] "Rhodothermaceae"                                   
## [395] "Magnetospira"                                      
## [396] "Blfdi19"                                           
## [397] "Turneriella"                                       
## [398] "Eel-36e1D6"                                        
## [399] "Gemmataceae"                                       
## [400] "C86"                                               
## [401] "Nitrospira"                                        
## [402] "Imperialibacter"                                   
## [403] "Niveispirillum"                                    
## [404] "Nitrosopumilaceae"                                 
## [405] "MBAE14"                                            
## [406] "Polyangiales"                                      
## [407] "Thiotrichaceae"                                    
## [408] "WCHB1-41"                                          
## [409] "Coxiella"                                          
## [410] "Desulfospira"                                      
## [411] "Candidatus Thiobios"                               
## [412] "Massilia"                                          
## [413] "Ga0077536"                                         
## [414] "Subgroup 9"                                        
## [415] "Thalassobaculum"                                   
## [416] "Nisaea"                                            
## [417] "Taibaiella"                                        
## [418] "Owenweeksia"                                       
## [419] "Defluviicoccales"                                  
## [420] "Plesiocystis"                                      
## [421] "Azospira"                                          
## [422] "Marinomonas"                                       
## [423] "Pleionea"                                          
## [424] "Colwelliaceae"                                     
## [425] "SAR202 clade"                                      
## [426] "Ascidiaceihabitans"                                
## [427] "Desulfobacterota"                                  
## [428] "Synechococcales"                                   
## [429] "Roseimarinus"                                      
## [430] "Flavobacterium"                                    
## [431] "Desulfosarcina"                                    
## [432] "Rubinisphaera"                                     
## [433] "Rickettsiaceae"                                    
## [434] "B2M28"                                             
## [435] "Haliea"                                            
## [436] "Ponticaulis"                                       
## [437] "Paraclostridium"                                   
## [438] "Desulfobacteraceae"                                
## [439] "Caenarcaniphilales"                                
## [440] "Acholeplasma"                                      
## [441] "Clostridium sensu stricto 1"                       
## [442] "Rhizobiaceae"                                      
## [443] "Sedimenticolaceae"                                 
## [444] "Gramella"                                          
## [445] "Weissella"                                         
## [446] "Bacillaceae"                                       
## [447] "Parvularcula"                                      
## [448] "Verrucomicrobiae"                                  
## [449] "SGST604"                                           
## [450] "Puniceispirillales"                                
## [451] "BBMC-4"                                            
## [452] "Roseibium"                                         
## [453] "CL500-3"                                           
## [454] "Candidatus Hepatoplasma"                           
## [455] "Subgroup 10"                                       
## [456] "Cyanobiaceae"                                      
## [457] "Altererythrobacter"                                
## [458] "Burkholderiales"                                   
## [459] "SCGC AAA286-E23"                                   
## [460] "Xenococcaceae"                                     
## [461] "Lachnospiraceae"                                   
## [462] "SS1-B-02-17"                                       
## [463] "Perexilibacter"                                    
## [464] "MBNT15"                                            
## [465] "Spirochaetaceae"                                   
## [466] "Prochlorothrix PCC-9006"                           
## [467] "Algisphaera"                                       
## [468] "Cyanothece PCC-8801"                               
## [469] "Aliagarivorans"                                    
## [470] "Helicobacteraceae"                                 
## [471] "Lactococcus"                                       
## [472] "Hydrogenedensaceae"                                
## [473] "0319-6G20"                                         
## [474] "SB-5"                                              
## [475] "Lactobacillus"                                     
## [476] "Spirochaetota"                                     
## [477] "Mangrovimonas"                                     
## [478] "Aliicoccus"                                        
## [479] "Cloacimonadales"                                   
## [480] "Formosa"                                           
## [481] "Schizothrix LEGE 07164"                            
## [482] "Hydrogenophaga"                                    
## [483] "Sphingomonas"                                      
## [484] "Sphingobacterium"                                  
## [485] "Izemoplasmataceae"                                 
## [486] "Latilactobacillus"                                 
## [487] "Lyngbya PCC-7419"                                  
## [488] "Oceanicella"                                       
## [489] "AB1"                                               
## [490] "Cyanobacteriaceae"                                 
## [491] "Marinicaulis"                                      
## [492] "Sphingobacteriales"                                
## [493] "Chryseobacterium"                                  
## [494] "Psychroflexus"                                     
## [495] "Methylacidiphilaceae"                              
## [496] "Sediminitomix"                                     
## [497] "Lineage IIb"                                       
## [498] "Marinobacter"                                      
## [499] "IheB2-31"                                          
## [500] "Sunxiuqinia"                                       
## [501] "Fusobacterium"                                     
## [502] "Caulobacter"                                       
## [503] "Clostridia"                                        
## [504] "Rheinheimera"                                      
## [505] "Fibrobacterales"                                   
## [506] "Run-SP154"                                         
## [507] "MD3-55"                                            
## [508] "Roseofilum AO1-A"                                  
## [509] "Maricaulis"                                        
## [510] "Candidatus Kaiserbacteria"                         
## [511] "Desulfofrigus"                                     
## [512] "Beggiatoa"                                         
## [513] "OM43 clade"                                        
## [514] "Alloprevotella"                                    
## [515] "Verrucomicrobiales"                                
## [516] "LCP-80"                                            
## [517] "Prevotella_7"                                      
## [518] "Desulfotalea"                                      
## [519] "Desulfatibacillum"                                 
## [520] "Zunongwangia"                                      
## [521] "EPR3968-O8a-Bc78"                                  
## [522] "Marinimicrobium"                                   
## [523] "Varunaivibrio"                                     
## [524] "Dietzia"                                           
## [525] "Algivirga"                                         
## [526] "Chalicogloea CCALA 975"                            
## [527] "Kocuria"                                           
## [528] "Micrococcus"                                       
## [529] "EV818SWSAP88"                                      
## [530] "Desertifilaceae"                                   
## [531] "Natranaerovirga"                                   
## [532] "Luteolibacter"                                     
## [533] "Acanthopleuribacter"                               
## [534] "Polyangia"                                         
## [535] "Bernardetiaceae"                                   
## [536] "Magnetospiraceae"                                  
## [537] "Prolixibacteraceae"                                
## [538] "Lutimaribacter"                                    
## [539] "Stenotrophomonas"                                  
## [540] "Woesearchaeales"                                   
## [541] "Schleiferia"                                       
## [542] "Sarcina"                                           
## [543] "Skermanella"                                       
## [544] "Subgroup 26"                                       
## [545] "Anaerolineae"                                      
## [546] "PAUC26f"                                           
## [547] "Sungkyunkwania"                                    
## [548] "MSBL3"                                             
## [549] "Lishizhenia"                                       
## [550] "Auricoccus-Abyssicoccus"                           
## [551] "SUP05 cluster"                                     
## [552] "Marinifilaceae"                                    
## [553] "endosymbionts"                                     
## [554] "Desulforhopalus"                                   
## [555] "Candidatus Hepatincola"                            
## [556] "Roseitalea"                                        
## [557] "Legionella"                                        
## [558] "Actibacterium"                                     
## [559] "Pseudovibrio"                                      
## [560] "NS2b marine group"                                 
## [561] "IMCC26256"                                         
## [562] "Simkaniaceae"                                      
## [563] "028H05-P-BN-P5"                                    
## [564] "Bradyrhizobium"                                    
## [565] "Caldithrix"                                        
## [566] "Cystobacter"                                       
## [567] "EC3"                                               
## [568] "mle1-8"                                            
## [569] "oc32"                                              
## [570] "Bacillus"                                          
## [571] "SM1A07"                                            
## [572] "Bosea"                                             
## [573] "Prochlorococcus MIT9313"                           
## [574] "Hirschia"                                          
## [575] "Vallitalea"                                        
## [576] "SCGC AAA011-D5"                                    
## [577] "Latescibacteraceae"                                
## [578] "Thalassomonas"                                     
## [579] "MB11C04 marine group"                              
## [580] "Paeniclostridium"                                  
## [581] "Blattabacteriaceae"                                
## [582] "Desulfopila"                                       
## [583] "Alkanindiges"                                      
## [584] "Xanthomonadales"                                   
## [585] "3PJM14"                                            
## [586] "Fibrobacteraceae"                                  
## [587] "Geopsychrobacteraceae"                             
## [588] "Tepidibacter"                                      
## [589] "Thermicanus"                                       
## [590] "Calditrichaceae"                                   
## [591] "Magnetococcus"                                     
## [592] "Blastococcus"                                      
## [593] "Temperatibacter"                                   
## [594] "Algibacter"                                        
## [595] "KD4-96"                                            
## [596] "Neiella"                                           
## [597] "1174-901-12"                                       
## [598] "Mastigocoleus BC008"                               
## [599] "Chromatiaceae"                                     
## [600] "Bacteroides"                                       
## [601] "Planomicrobium"                                    
## [602] "Symploca PCC-8002"                                 
## [603] "Aerococcaceae"                                     
## [604] "Spongiimonas"                                      
## [605] "Gemmobacter"                                       
## [606] "Labilibacter"                                      
## [607] "Zoogloea"                                          
## [608] "Rubrobacter"                                       
## [609] "Pelomonas"                                         
## [610] "Modestobacter"                                     
## [611] "Dstr-E11"                                          
## [612] "Silvanigrellaceae"                                 
## [613] "Pelagicoccus"                                      
## [614] "Sulfurimonas"                                      
## [615] "Thermoplasmata"                                    
## [616] "4572-13"                                           
## [617] "Rapidithrix"                                       
## [618] "Aureimarina"                                       
## [619] "Moorea"                                            
## [620] "Sediminispirochaeta"                               
## [621] "Robertkochia"                                      
## [622] "Subgroup 17"                                       
## [623] "Candidatus Peribacteria"                           
## [624] "Desulfatiglans"                                    
## [625] "Virgibacillus"                                     
## [626] "Actinomyces"                                       
## [627] "Balneolaceae"                                      
## [628] "Paracaedibacteraceae"                              
## [629] "Gaiellales"                                        
## [630] "Saccharopolyspora"                                 
## [631] "Sulfurospirillum"                                  
## [632] "Anaerococcus"                                      
## [633] "Candidatus Electrothrix"                           
## [634] "Lachnospiraceae UCG-007"                           
## [635] "Candidatus Campbellbacteria"                       
## [636] "Aestuariicella"                                    
## [637] "SZB85"                                             
## [638] "Exiguobacterium"                                   
## [639] "Thermomonas"                                       
## [640] "XY-R5"                                             
## [641] "GW2011_GWC1_47_15"                                 
## [642] "Geitlerinema PCC-9228"                             
## [643] "Candidatus Omnitrophus"                            
## [644] "Zixibacteria"                                      
## [645] "EF100-94H03"                                       
## [646] "Desulfobacterales"                                 
## [647] "Lentimicrobiaceae"                                 
## [648] "Deferribacteraceae"                                
## [649] "SG8-4"                                             
## [650] "Truepera"                                          
## [651] "Novosphingobium"                                   
## [652] "FW22"                                              
## [653] "Hellea"                                            
## [654] "AncK6"                                             
## [655] "Candidatus Entotheonella"                          
## [656] "Diplosphaera"                                      
## [657] "Pseudomonadaceae"                                  
## [658] "Porphyromonas"                                     
## [659] "Clostridiales"                                     
## [660] "Deinococcaceae"                                    
## [661] "Amphiplicatus"                                     
## [662] "Opitutaceae"                                       
## [663] "Desulfococcaceae"                                  
## [664] "Kytococcus"                                        
## [665] "Gemella"                                           
## [666] "Magnetovibrio"                                     
## [667] "Candidatus Nitrosopelagicus"                       
## [668] "Snodgrassella"                                     
## [669] "Chroococcidiopsaceae"                              
## [670] "Hyphomicrobium"                                    
## [671] "Candidatus Endoecteinascidia"                      
## [672] "koll11"                                            
## [673] "Lineage IV"                                        
## [674] "Salinispira"                                       
## [675] "Kiloniellales"                                     
## [676] "Thiovulum"                                         
## [677] "Dermacoccaceae"                                    
## [678] "Planctomicrobium"                                  
## [679] "Aquimarina"                                        
## [680] "Rivularia PCC-7116"                                
## [681] "Nostocaceae"                                       
## [682] "Francisella"                                       
## [683] "Saccharibacillus"                                  
## [684] "Mycoplasmataceae"                                  
## [685] "Methyloligellaceae"                                
## [686] "Hyphomicrobiaceae"                                 
## [687] "PAUC34f"                                           
## [688] "Anaerolineaceae"                                   
## [689] "Aestuariispira"                                    
## [690] "CCM11a"                                            
## [691] "Kiritimatiellaceae"                                
## [692] "Akkermansia"                                       
## [693] "Bauldia"                                           
## [694] "AKAU3564 sediment group"                           
## [695] "[Synechococcus] spongiarum group"                  
## [696] "Lenti-02"                                          
## [697] "Clostridia vadinBB60 group"                        
## [698] "Firmicutes"                                        
## [699] "LS-NOB"                                            
## [700] "Candidatus Gigarickettsia"                         
## [701] "Arctic97B-4 marine group"                          
## [702] "Oceanicoccus"                                      
## [703] "Gemmatimonadaceae"                                 
## [704] "TG3"                                               
## [705] "Desulfobulbaceae"                                  
## [706] "Mycobacterium"                                     
## [707] "Geobacillus"                                       
## [708] "TSBb06"                                            
## [709] "S-70"                                              
## [710] "Vampirovibrionales"                                
## [711] "Brevundimonas"                                     
## [712] "Celerinatantimonas"                                
## [713] "Abyssivirga"                                       
## [714] "Candidatus Peregrinibacteria"                      
## [715] "Syntrophobacterales"                               
## [716] "Xenococcus CRM"                                    
## [717] "Steroidobacteraceae"                               
## [718] "Solimonadaceae"                                    
## [719] "Candidatus Photodesmus"                            
## [720] "Malaciobacter"                                     
## [721] "Psychromonas"                                      
## [722] "Hypnocyclicus"                                     
## [723] "Cylindrospermopsis CRJ1"                           
## [724] "MVP-88"                                            
## [725] "Candidatus Captivus"                               
## [726] "Entotheonellaceae"                                 
## [727] "Candidatus Paenicardinium"                         
## [728] "ZOR0006"                                           
## [729] "Margulisbacteria"                                  
## [730] "Muriicola"                                         
## [731] "Aureibacter"                                       
## [732] "Archaea"                                           
## [733] "Halanaerobium"                                     
## [734] "Candidatus Magasanikbacteria"                      
## [735] "Parcubacteria"                                     
## [736] "AqS1"                                              
## [737] "Pseudarthrobacter"                                 
## [738] "Subgroup 21"                                       
## [739] "Nguyenibacter"                                     
## [740] "Maritalea"                                         
## [741] "Salinimicrobium"                                   
## [742] "Paramaledivibacter"                                
## [743] "Nitrospirillum"                                    
## [744] "Saccharimonadales"                                 
## [745] "Victivallales"                                     
## [746] "Alsobacter"                                        
## [747] "Brevibacillus"                                     
## [748] "Luteivirga"                                        
## [749] "Gastranaerophilales"                               
## [750] "MidBa8"                                            
## [751] "Cellulomonas"                                      
## [752] "Saccharicrinis"                                    
## [753] "Thermodesulfovibrionia"                            
## [754] "Spongiispira"                                      
## [755] "Nitrosomonas"                                      
## [756] "I3A"                                               
## [757] "LCP-89"                                            
## [758] "Diplorickettsiaceae"                               
## [759] "Oceaniserpentilla"                                 
## [760] "Microtrichaceae"                                   
## [761] "Sphingosinicella"                                  
## [762] "Spirochaeta"                                       
## [763] "Subgroup 23"                                       
## [764] "Omnitrophales"                                     
## [765] "Phenylobacterium"                                  
## [766] "MSB-3C8"                                           
## [767] "Chitinivibrio"                                     
## [768] "PAUC43f marine benthic group"                      
## [769] "Izimaplasma"                                       
## [770] "Kangiellaceae"                                     
## [771] "Marinilabiliaceae"                                 
## [772] "Cyanobacteria"                                     
## [773] "Catenovulum"                                       
## [774] "Chloroflexi"                                       
## [775] "Nevskia"                                           
## [776] "Omnitrophia"                                       
## [777] "mle1-7"                                            
## [778] "Endomicrobium"                                     
## [779] "Candidatus Falkowbacteria"                         
## [780] "Oceanibaculum"                                     
## [781] "Pseudobowmanella"                                  
## [782] "Lokiarchaeia"                                      
## [783] "Actinoplanes"                                      
## [784] "Blautia"                                           
## [785] "D90"                                               
## [786] "Simkania"                                          
## [787] "pItb-vmat-80"                                      
## [788] "Prevotella"                                        
## [789] "Calothrix PCC-6303"                                
## [790] "Micrococcales"                                     
## [791] "Ectothiorhodospiraceae"                            
## [792] "Finegoldia"                                        
## [793] "Cnuella"                                           
## [794] "Bacteroidota"                                      
## [795] "MVP-21"                                            
## [796] "GWE2-31-10"                                        
## [797] "Rikenella"                                         
## [798] "Candidatus Altiarchaeum"                           
## [799] "Acetobacteraceae"                                  
## [800] "BD72BR169"                                         
## [801] "Thioploca"                                         
## [802] "Devosia"                                           
## [803] "Milano-WF1B-44"                                    
## [804] "Catenococcus"                                      
## [805] "Marine Benthic Group D and DHVEG-1"                
## [806] "Cenarchaeum"                                       
## [807] "Chthoniobacterales"                                
## [808] "Flavisolibacter"                                   
## [809] "Clostridia UCG-014"                                
## [810] "EC94"                                              
## [811] "Bryobacter"                                        
## [812] "Cytophaga"                                         
## [813] "Cocleimonas"                                       
## [814] "Moduliflexaceae"                                   
## [815] "SAR11 clade"                                       
## [816] "Pseudanabaena PCC-7367"                            
## [817] "Defluviitaleaceae UCG-011"                         
## [818] "RF39"                                              
## [819] "Aliifodinibius"                                    
## [820] "Clade I"                                           
## [821] "Victivallaceae"                                    
## [822] "Sva0485"                                           
## [823] "Brachybacterium"                                   
## [824] "Litorilituus"                                      
## [825] "Campylobacter"                                     
## [826] "Pseudoxanthomonas"                                 
## [827] "MSBL8"                                             
## [828] "Hydrogenispora"                                    
## [829] "Patescibacteria"                                   
## [830] "SPG12-343-353-B69"                                 
## [831] "MVP-15"                                            
## [832] "Thalassospira"                                     
## [833] "Ellin6055"                                         
## [834] "Sumerlaeia"                                        
## [835] "Clade Ib"                                          
## [836] "Pelobacter"                                        
## [837] "Peptoniphilus"                                     
## [838] "Marispirillum"                                     
## [839] "cvE6"                                              
## [840] "Beijerinckiaceae"                                  
## [841] "Alterococcus"                                      
## [842] "Solirubrobacter"                                   
## [843] "Methylotenera"                                     
## [844] "Ardenticatenaceae"                                 
## [845] "Agrilactobacillus"                                 
## [846] "Pedobacter"                                        
## [847] "Tahibacter"                                        
## [848] "Arcticibacter"                                     
## [849] "Lentimicrobium"                                    
## [850] "Deep Sea Euryarchaeotic Group(DSEG)"               
## [851] "Aestuariimonas"                                    
## [852] "Incertae Sedis"                                    
## [853] "Chthoniobacteraceae"                               
## [854] "Oligosphaeraceae"                                  
## [855] "SBYC"                                              
## [856] "BD2-3"                                             
## [857] "LD1-PA32"                                          
## [858] "Kiritimatiellae"                                   
## [859] "Defluviitaleaceae"                                 
## [860] "Rs-M59 termite group"                              
## [861] "Hymenobacter"                                      
## [862] "Conexibacter"                                      
## [863] "Candidatus Fritschea"                              
## [864] "Brevinema"                                         
## [865] "SM23-30"                                           
## [866] "Gemmata"                                           
## [867] "Ga0074140"                                         
## [868] "Salinihabitans"                                    
## [869] "AKIW781"                                           
## [870] "Erysipelotrichales"                                
## [871] "OPB41"                                             
## [872] "Leuconostoc"                                       
## [873] "Acidobacteriae"                                    
## [874] "Lactiplantibacillus"                               
## [875] "Candidatus Terrybacteria"                          
## [876] "Sulfurovum"                                        
## [877] "PeM15"                                             
## [878] "Buchnera"                                          
## [879] "Bathyarchaeia"                                     
## [880] "Isosphaeraceae"
```

``` r
get_taxa_unique(ps5, "Genus") 
```

```
##  [1] "Vibrio"                                            
##  [2] "Clade III"                                         
##  [3] "Clade Ia"                                          
##  [4] "Proteobacteria"                                    
##  [5] "Candidatus Actinomarina"                           
##  [6] "Sphingopyxis"                                      
##  [7] "Synechococcus CC9902"                              
##  [8] "Litoricola"                                        
##  [9] "AEGEAN-169 marine group"                           
## [10] "Delftia"                                           
## [11] "Blastocatellaceae"                                 
## [12] "Pseudomonas"                                       
## [13] "Clade II"                                          
## [14] "HIMB11"                                            
## [15] "Rubritalea"                                        
## [16] "SAR116 clade"                                      
## [17] "SAR86 clade"                                       
## [18] "NS5 marine group"                                  
## [19] "Campylobacterales"                                 
## [20] "OM60(NOR5) clade"                                  
## [21] "S25-593"                                           
## [22] "Cryomorphaceae"                                    
## [23] "NS4 marine group"                                  
## [24] "NS11-12 marine group"                              
## [25] "Urania-1B-19 marine sediment group"                
## [26] "Flavobacteriaceae"                                 
## [27] "Candidatus Puniceispirillum"                       
## [28] "NS9 marine group"                                  
## [29] "Bacteria"                                          
## [30] "Pseudomonadales"                                   
## [31] "MBMPE27"                                           
## [32] "Blastocatella"                                     
## [33] "Fokiniaceae"                                       
## [34] "Balneola"                                          
## [35] "Pseudoalteromonas"                                 
## [36] "Marinimicrobia (SAR406 clade)"                     
## [37] "Alphaproteobacteria"                               
## [38] "Rhodopirellula"                                    
## [39] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
## [40] "Francisellaceae"                                   
## [41] "PS1 clade"                                         
## [42] "Bradymonadales"                                    
## [43] "Staphylococcus"                                    
## [44] "Azotobacter"                                       
## [45] "Pedosphaeraceae"                                   
## [46] "Alteromonas"                                       
## [47] "Candidatus Amoebophilus"                           
## [48] "Coraliomargarita"                                  
## [49] "Cohaesibacter"                                     
## [50] "Candidatus Cryptoprodotis"                         
## [51] "Rickettsiales"                                     
## [52] "Comamonadaceae"                                    
## [53] "Ekhidna"                                           
## [54] "Acinetobacter"                                     
## [55] "Marine Methylotrophic Group 3"                     
## [56] "Marine Group II"                                   
## [57] "Dyella"                                            
## [58] "R76-B128"                                          
## [59] "SCGC AAA164-E04"                                   
## [60] "Gammaproteobacteria"                               
## [61] "Propionigenium"
```

``` r
get_taxa_unique(ps10, "Genus") 
```

```
##  [1] "Clade III"                                         
##  [2] "Clade Ia"                                          
##  [3] "Proteobacteria"                                    
##  [4] "Candidatus Actinomarina"                           
##  [5] "Sphingopyxis"                                      
##  [6] "Synechococcus CC9902"                              
##  [7] "Litoricola"                                        
##  [8] "AEGEAN-169 marine group"                           
##  [9] "Delftia"                                           
## [10] "Blastocatellaceae"                                 
## [11] "Pseudomonas"                                       
## [12] "Clade II"                                          
## [13] "HIMB11"                                            
## [14] "SAR116 clade"                                      
## [15] "SAR86 clade"                                       
## [16] "Vibrio"                                            
## [17] "NS5 marine group"                                  
## [18] "Campylobacterales"                                 
## [19] "OM60(NOR5) clade"                                  
## [20] "S25-593"                                           
## [21] "Cryomorphaceae"                                    
## [22] "NS4 marine group"                                  
## [23] "NS11-12 marine group"                              
## [24] "Urania-1B-19 marine sediment group"                
## [25] "Flavobacteriaceae"                                 
## [26] "Candidatus Puniceispirillum"                       
## [27] "NS9 marine group"                                  
## [28] "Bacteria"                                          
## [29] "Pseudomonadales"                                   
## [30] "MBMPE27"                                           
## [31] "Blastocatella"                                     
## [32] "Fokiniaceae"                                       
## [33] "Balneola"                                          
## [34] "Pseudoalteromonas"                                 
## [35] "Marinimicrobia (SAR406 clade)"                     
## [36] "Alphaproteobacteria"                               
## [37] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
## [38] "Francisellaceae"                                   
## [39] "PS1 clade"                                         
## [40] "Staphylococcus"                                    
## [41] "Azotobacter"                                       
## [42] "Pedosphaeraceae"                                   
## [43] "Cohaesibacter"                                     
## [44] "Candidatus Cryptoprodotis"                         
## [45] "Rickettsiales"                                     
## [46] "Comamonadaceae"
```

``` r
# Now export filtered otu and taxa tables from phyloseq for future reference
otu_ps5 = as(otu_table(ps5), "matrix")
taxon_ps5 = as(tax_table(ps5), "matrix")
metadata = as(sample_data(ps5), "matrix")
write.table(otu_ps5,"silva_nochloronomito_otu_table_ps5.txt",sep="\t",col.names=NA)
write.table(taxon_ps5,"silva_nochloronomito_taxa_table_ps5.txt",sep="\t",col.names=NA)
write.table(metadata,"metadata_ps5.txt",sep="\t",col.names=NA) 

#relative abundance
ps5_ra<-transform_sample_counts(ps5, function(OTU) OTU/sum(OTU))
otu_ps5_ra = as(otu_table(ps5_ra), "matrix")
write.table(otu_ps5_ra,"silva_nochloronomito_otu_table_ps5_RA.txt",sep="\t",col.names=NA)
```

# Alpha Diversity on unfiltered dataset

``` r
# load in data for 35 samples, no low abundance taxa removed
otu <- read.table("silva_nochloronomito_otu_table_ps1.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps1.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps1.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps1 <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 6137 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 6137 taxa by 6 taxonomic ranks ]
```

``` r
# calculate alpha diversity measures with phyloseq wrapper
alpha.diversity <- estimate_richness(ps1, measures = c("Observed","Shannon","Simpson"))

# combine alpha diversity measures with sample metadata
alpha.diversity <- tibble::rownames_to_column(alpha.diversity, "sample")
# fix sample names
alpha.diversity$sample<- gsub("X", "", alpha.diversity$sample)
alpha.diversity$sample<- gsub(".", "-", alpha.diversity$sample,fixed=TRUE)
samples <- tibble::rownames_to_column(samples, "sample")
alpha_meta <- merge(alpha.diversity, samples, by="sample", all=TRUE)

# set aesthetics
theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
alpha_meta$type<-factor(alpha_meta$type,levels=c("before","after"))
cols<-c("before"="#6CD5D9","after"="#107A86")


# plot box and whiskers with points on top
obs <- ggplot(alpha_meta,aes(x=type,y=Observed))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(fill=type,shape=fraction))+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank())+
  theme(text=element_text(size=14,face="bold"))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=cols)

shan <- ggplot(alpha_meta,aes(x=type,y=Shannon))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(fill=type,shape=fraction))+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank())+
  theme(text=element_text(size=14,face="bold"))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=cols)

simp <- ggplot(alpha_meta,aes(x=type,y=Simpson))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(aes(fill=type,shape=fraction))+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank())+
  theme(text=element_text(size=14,face="bold"))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=cols)

#pdf("AlphaDiversity.pdf", bg ="white", width=8,height=4)
plot_grid(obs,shan,simp, labels=c("A","B","C"), ncol=3)
```

<img src="../figures/16S-fig-unnamed-chunk-3-1.png" width="672" />

``` r
#dev.off()
```


# Test for statistical differences in alpha diversity.


``` r
# First test for normality using Shapiro test. If the value of p is equal to or less than 0.05, then the hypothesis of normality will be rejected. 
shapiro.test(alpha_meta$Observed) # normality rejected
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  alpha_meta$Observed
## W = 0.87757, p-value = 0.001041
```

``` r
# 2-way ANOVA on ranked data for non-normally distributed diversity measures
anova1 <- aov(rank(Observed) ~ type, alpha_meta) 
summary(anova1)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)  
## type         1    405   404.9   4.223 0.0479 *
## Residuals   33   3165    95.9                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# significant at p=0.01 level

anova2 <- aov(rank(Observed) ~ fraction, alpha_meta) 
summary(anova2)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)
## fraction     1    148   148.2    1.43   0.24
## Residuals   33   3421   103.7
```

``` r
# n.s.
```




``` r
# First test for normality using Shapiro test. If the value of p is equal to or less than 0.05, then the hypothesis of normality will be rejected. 
shapiro.test(alpha_meta$Shannon)
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  alpha_meta$Shannon
## W = 0.94362, p-value = 0.07218
```

``` r
# 2-way ANOVA test for normally distributed diversity measure - Shannon diversity
anova3 <- aov(Shannon ~ type, alpha_meta) 
summary(anova3)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)
## type         1  0.110  0.1105    0.48  0.493
## Residuals   33  7.596  0.2302
```

``` r
# not sign. diff.

anova4 <- aov(Shannon ~ fraction, alpha_meta) 
summary(anova4)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)  
## fraction     1  0.860  0.8602   4.146 0.0498 *
## Residuals   33  6.846  0.2075                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# significant at p=0.01 level
```



``` r
# First test for normality using Shapiro test. If the value of p is equal to or less than 0.05, then the hypothesis of normality will be rejected. 
shapiro.test(alpha_meta$Simpson) # normality rejected
```

```
## 
## 	Shapiro-Wilk normality test
## 
## data:  alpha_meta$Simpson
## W = 0.78493, p-value = 1.048e-05
```

``` r
# 2-way ANOVA on ranked data for non-normally distributed diversity measures
anova5 <- aov(rank(Simpson) ~ type, alpha_meta) 
summary(anova5)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)
## type         1      0    0.46   0.004  0.949
## Residuals   33   3570  108.17
```

``` r
# not sign. diff.

anova6 <- aov(rank(Simpson) ~ fraction, alpha_meta) 
summary(anova6)
```

```
##             Df Sum Sq Mean Sq F value Pr(>F)  
## fraction     1  439.7   439.7   4.635 0.0387 *
## Residuals   33 3130.3    94.9                 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
# significant at p=0.01 level
```

# Beta Diversity (Perform center-log-ratio transformation on ASVs and calculate Aitchison Distance and principal components)


``` r
# load in data and create phyloseq object
otu <- read.table("silva_nochloronomito_otu_table_ps5.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 112 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 112 taxa by 6 taxonomic ranks ]
```

``` r
# First, replace 0 values with an estimate (because normalization is taking log, can't have 0)
# Also transposing here, need samples as rows
d.czm <- cmultRepl(t(otu), method="CZM", label=0, z.warning=1)
```

```
## No. adjusted imputations:  312
```

``` r
# Perform the center-log-ratio (CLR) transformation 
d.clr <- codaSeq.clr(d.czm)
# transpose matrix of CLR transformed data for ordination and dendrogram
E.clr <- t(d.clr)
# plot compositional PCA biplot (perform a singular value decomposition)
d.pcx <- prcomp(E.clr)
# calculate percent variance explained for the axis labels
pc1 <- round(d.pcx$sdev[1]^2/sum(d.pcx$sdev^2),2)
pc2 <- round(d.pcx$sdev[2]^2/sum(d.pcx$sdev^2),2)
xlab <- paste("PC1: ", pc1, sep="")
ylab <- paste("PC2: ", pc2, sep="")
#biplot(d.pcx, cex=c(0.6,0.4), var.axes=F,scale=1, xlab=xlab, ylab=ylab)
summary(d.pcx)
```

```
## Importance of components:
##                            PC1     PC2     PC3    PC4    PC5     PC6     PC7
## Standard deviation     14.9061 6.81851 5.64438 5.2543 4.7606 4.30324 3.66796
## Proportion of Variance  0.4137 0.08656 0.05932 0.0514 0.0422 0.03448 0.02505
## Cumulative Proportion   0.4137 0.50026 0.55958 0.6110 0.6532 0.68765 0.71270
##                            PC8     PC9    PC10   PC11    PC12    PC13    PC14
## Standard deviation     3.53436 3.48748 3.36524 3.1946 3.14683 3.08590 2.89233
## Proportion of Variance 0.02326 0.02265 0.02109 0.0190 0.01844 0.01773 0.01558
## Cumulative Proportion  0.73596 0.75861 0.77969 0.7987 0.81713 0.83486 0.85044
##                           PC15    PC16    PC17    PC18    PC19    PC20    PC21
## Standard deviation     2.80305 2.74131 2.61766 2.49548 2.43424 2.32117 2.26858
## Proportion of Variance 0.01463 0.01399 0.01276 0.01159 0.01103 0.01003 0.00958
## Cumulative Proportion  0.86507 0.87906 0.89182 0.90341 0.91445 0.92448 0.93406
##                           PC22    PC23    PC24    PC25    PC26    PC27    PC28
## Standard deviation     2.15156 2.07123 1.93513 1.86187 1.75575 1.67392 1.65349
## Proportion of Variance 0.00862 0.00799 0.00697 0.00645 0.00574 0.00522 0.00509
## Cumulative Proportion  0.94268 0.95067 0.95764 0.96409 0.96983 0.97505 0.98014
##                          PC29    PC30    PC31    PC32    PC33   PC34     PC35
## Standard deviation     1.5719 1.52521 1.32589 1.27627 1.16377 1.0621 6.62e-15
## Proportion of Variance 0.0046 0.00433 0.00327 0.00303 0.00252 0.0021 0.00e+00
## Cumulative Proportion  0.9847 0.98907 0.99235 0.99538 0.99790 1.0000 1.00e+00
```

``` r
str(d.pcx)
```

```
## List of 5
##  $ sdev    : num [1:35] 14.91 6.82 5.64 5.25 4.76 ...
##  $ rotation: num [1:112, 1:35] -0.0355 -0.1387 -0.148 -0.1505 -0.0679 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:112] "TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAATTGCATTTGAAACT"| __truncated__ "TACGAAGGGACCTAGCGTAATTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTAAGTAAGTTAATTGTGAAAGCCCGAAGCTCAACTTCGGAATTGCAATTAAAACT"| __truncated__ "TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACT"| __truncated__ "TACGAAGGGACCTAGCGTAATTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTAAGTAAGTTAATTGTGAAAGCCCAAAGCTCAACTTTGGAATTGCAATTAAAACT"| __truncated__ ...
##   .. ..$ : chr [1:35] "PC1" "PC2" "PC3" "PC4" ...
##  $ center  : Named num [1:112] -3.81e-17 3.93e-16 1.14e-16 -6.34e-17 3.30e-16 ...
##   ..- attr(*, "names")= chr [1:112] "TACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCATGCAGGTGGTTTGTTAAGTCAGATGTGAAAGCCCGGGGCTCAACCTCGGAATTGCATTTGAAACT"| __truncated__ "TACGAAGGGACCTAGCGTAATTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTAAGTAAGTTAATTGTGAAAGCCCGAAGCTCAACTTCGGAATTGCAATTAAAACT"| __truncated__ "TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACT"| __truncated__ "TACGAAGGGACCTAGCGTAATTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGTTAAGTAAGTTAATTGTGAAAGCCCAAAGCTCAACTTTGGAATTGCAATTAAAACT"| __truncated__ ...
##  $ scale   : logi FALSE
##  $ x       : num [1:35, 1:35] -5.7 10.04 -20.29 14.01 3.94 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:35] "10A-DNA" "10A-rna" "10H-DNA" "10H-rna" ...
##   .. ..$ : chr [1:35] "PC1" "PC2" "PC3" "PC4" ...
##  - attr(*, "class")= chr "prcomp"
```

``` r
screeplot(d.pcx)
```

<img src="../figures/16S-fig-unnamed-chunk-7-1.png" width="672" />

``` r
# Make a pretty PCA plot with ggplot
df_out <- as.data.frame(d.pcx$x)
theme_set(theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))
cols<-c("before"="#6CD5D9","after"="#107A86")
samples$type<-factor(samples$type, levels=c("before","after"))
#pdf("PCA.pdf",bg ="white",width=8, height=4)
p<-ggplot(df_out,aes(x=PC1,y=PC2,fill=samples$type,shape=samples$fraction))
p<-p+geom_point(size=3)+
  theme(axis.title = element_text(size=14))+
  theme(axis.text=element_text(size=12))+
  theme(legend.title = element_text(size=14))+
  theme(legend.text = element_text(size=12))+
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=cols)+
  guides(fill = guide_legend(override.aes=list(shape=21)))+
  facet_grid(~samples$fraction)+
  ylim(-35,35)+
  xlim(-35,35)+
  coord_fixed(ratio=1)
p + labs(x=xlab, y=ylab, fill="Type",shape="Fraction") 
```

<img src="../figures/16S-fig-unnamed-chunk-7-2.png" width="672" />

``` r
#dev.off()

####### Use phyloseq/vegan to perform PERMANOVA
# set metadata as factors
type<-as.character(samples$type)
frac<-as.character(samples$fraction)
# permanova between groups using Aitchison distance
dist.clr <- dist(E.clr)
perm<-adonis2(dist.clr~type*frac,as(sample_data(ps),"data.frame"),by="terms")
print(perm)
```

```
## Permutation test for adonis under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## adonis2(formula = dist.clr ~ type * frac, data = as(sample_data(ps), "data.frame"), by = "terms")
##           Df SumOfSqs      R2      F Pr(>F)    
## type       1   1335.6 0.07314 3.0383  0.012 *  
## frac       1   2680.6 0.14679 6.0981  0.001 ***
## type:frac  1    617.8 0.03383 1.4054  0.170    
## Residual  31  13627.0 0.74624                  
## Total     34  18260.9 1.00000                  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```



Stacked bar charts

``` r
# load in data and create phyloseq object
otu <- read.table("silva_nochloronomito_otu_table_ps5.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 112 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 112 taxa by 6 taxonomic ranks ]
```

``` r
ps_ra<-transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
get_taxa_unique(ps_ra, "Class")
```

```
##  [1] "Gammaproteobacteria"           "Alphaproteobacteria"          
##  [3] "Proteobacteria"                "Acidimicrobiia"               
##  [5] "Cyanobacteriia"                "Blastocatellia"               
##  [7] "Verrucomicrobiae"              "Bacteroidia"                  
##  [9] "Campylobacteria"               "Phycisphaerae"                
## [11] "Bacteria"                      "Rhodothermia"                 
## [13] "Marinimicrobia (SAR406 clade)" "Planctomycetes"               
## [15] "Desulfuromonadia"              "Bacilli"                      
## [17] "Thermoplasmata"                "Kiritimatiellae"              
## [19] "Fusobacteriia"
```

``` r
get_taxa_unique(ps_ra, "Order")
```

```
##  [1] "Enterobacterales"              "SAR11 clade"                  
##  [3] "Proteobacteria"                "Actinomarinales"              
##  [5] "Sphingomonadales"              "Synechococcales"              
##  [7] "Pseudomonadales"               "Rhodospirillales"             
##  [9] "Burkholderiales"               "Blastocatellales"             
## [11] "Rhodobacterales"               "Verrucomicrobiales"           
## [13] "Puniceispirillales"            "Flavobacteriales"             
## [15] "Campylobacterales"             "Rickettsiales"                
## [17] "Sphingobacteriales"            "Phycisphaerales"              
## [19] "Bacteria"                      "MBMPE27"                      
## [21] "Balneolales"                   "Marinimicrobia (SAR406 clade)"
## [23] "Alphaproteobacteria"           "Pirellulales"                 
## [25] "Rhizobiales"                   "Francisellales"               
## [27] "Parvibaculales"                "Bradymonadales"               
## [29] "Staphylococcales"              "Pedosphaerales"               
## [31] "Cytophagales"                  "Opitutales"                   
## [33] "Nitrosococcales"               "Marine Group II"              
## [35] "Xanthomonadales"               "Kiritimatiellales"            
## [37] "Gammaproteobacteria"           "Fusobacteriales"
```

``` r
get_taxa_unique(ps_ra, "Family")
```

```
##  [1] "Vibrionaceae"                  "Clade III"                    
##  [3] "Clade I"                       "Proteobacteria"               
##  [5] "Actinomarinaceae"              "Sphingomonadaceae"            
##  [7] "Cyanobiaceae"                  "Litoricolaceae"               
##  [9] "AEGEAN-169 marine group"       "Comamonadaceae"               
## [11] "Blastocatellaceae"             "Pseudomonadaceae"             
## [13] "Clade II"                      "Rhodobacteraceae"             
## [15] "Rubritaleaceae"                "SAR116 clade"                 
## [17] "SAR86 clade"                   "Flavobacteriaceae"            
## [19] "Campylobacterales"             "Halieaceae"                   
## [21] "S25-593"                       "Cryomorphaceae"               
## [23] "NS11-12 marine group"          "Phycisphaeraceae"             
## [25] "NS9 marine group"              "Bacteria"                     
## [27] "Pseudomonadales"               "MBMPE27"                      
## [29] "Fokiniaceae"                   "Balneolaceae"                 
## [31] "Pseudoalteromonadaceae"        "Marinimicrobia (SAR406 clade)"
## [33] "Alphaproteobacteria"           "Pirellulaceae"                
## [35] "Rhizobiaceae"                  "Francisellaceae"              
## [37] "PS1 clade"                     "Bradymonadales"               
## [39] "Staphylococcaceae"             "Pedosphaeraceae"              
## [41] "Alteromonadaceae"              "Amoebophilaceae"              
## [43] "Puniceicoccaceae"              "Rickettsiaceae"               
## [45] "Rickettsiales"                 "Cyclobacteriaceae"            
## [47] "Moraxellaceae"                 "Methylophagaceae"             
## [49] "Marine Group II"               "Rhodanobacteraceae"           
## [51] "Kiritimatiellaceae"            "Gammaproteobacteria"          
## [53] "Fusobacteriaceae"
```

``` r
get_taxa_unique(ps_ra, "Genus")
```

```
##  [1] "Vibrio"                                            
##  [2] "Clade III"                                         
##  [3] "Clade Ia"                                          
##  [4] "Proteobacteria"                                    
##  [5] "Candidatus Actinomarina"                           
##  [6] "Sphingopyxis"                                      
##  [7] "Synechococcus CC9902"                              
##  [8] "Litoricola"                                        
##  [9] "AEGEAN-169 marine group"                           
## [10] "Delftia"                                           
## [11] "Blastocatellaceae"                                 
## [12] "Pseudomonas"                                       
## [13] "Clade II"                                          
## [14] "HIMB11"                                            
## [15] "Rubritalea"                                        
## [16] "SAR116 clade"                                      
## [17] "SAR86 clade"                                       
## [18] "NS5 marine group"                                  
## [19] "Campylobacterales"                                 
## [20] "OM60(NOR5) clade"                                  
## [21] "S25-593"                                           
## [22] "Cryomorphaceae"                                    
## [23] "NS4 marine group"                                  
## [24] "NS11-12 marine group"                              
## [25] "Urania-1B-19 marine sediment group"                
## [26] "Flavobacteriaceae"                                 
## [27] "Candidatus Puniceispirillum"                       
## [28] "NS9 marine group"                                  
## [29] "Bacteria"                                          
## [30] "Pseudomonadales"                                   
## [31] "MBMPE27"                                           
## [32] "Blastocatella"                                     
## [33] "Fokiniaceae"                                       
## [34] "Balneola"                                          
## [35] "Pseudoalteromonas"                                 
## [36] "Marinimicrobia (SAR406 clade)"                     
## [37] "Alphaproteobacteria"                               
## [38] "Rhodopirellula"                                    
## [39] "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
## [40] "Francisellaceae"                                   
## [41] "PS1 clade"                                         
## [42] "Bradymonadales"                                    
## [43] "Staphylococcus"                                    
## [44] "Azotobacter"                                       
## [45] "Pedosphaeraceae"                                   
## [46] "Alteromonas"                                       
## [47] "Candidatus Amoebophilus"                           
## [48] "Coraliomargarita"                                  
## [49] "Cohaesibacter"                                     
## [50] "Candidatus Cryptoprodotis"                         
## [51] "Rickettsiales"                                     
## [52] "Comamonadaceae"                                    
## [53] "Ekhidna"                                           
## [54] "Acinetobacter"                                     
## [55] "Marine Methylotrophic Group 3"                     
## [56] "Marine Group II"                                   
## [57] "Dyella"                                            
## [58] "R76-B128"                                          
## [59] "SCGC AAA164-E04"                                   
## [60] "Gammaproteobacteria"                               
## [61] "Propionigenium"
```

``` r
sample_data(ps_ra)$type<-factor(sample_data(ps_ra)$type, levels=c("before","after"))

n <- 19
# after plotting, you can re-run the next line to create a different selection of colors
palette <- distinctColorPalette(n)

#pdf("barchart_Class.pdf",width=11)
p1=plot_bar(ps_ra, "colony" ,fill="Class")+
  facet_grid(fraction~type,scales="free",space="free")+
  geom_bar(aes(fill=Class), stat="identity",position="stack")+
  theme_bw()+
  theme(strip.text=element_text(face="bold", size=12))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  scale_fill_manual(values=palette)+
  #theme(axis.title.x = element_blank())+
  theme(legend.position = "bottom")
p1
```

<img src="../figures/16S-fig-unnamed-chunk-8-1.png" width="672" />

``` r
#dev.off()



n <- 38
palette <- distinctColorPalette(n)

#pdf("barchart_Order.pdf",width=11)
p1=plot_bar(ps_ra, "colony" ,fill="Order")+
  facet_grid(fraction~type,scales="free",space="free")+
  geom_bar(aes(fill=Order), stat="identity",position="stack")+
  theme_bw()+
  theme(strip.text=element_text(face="bold", size=12))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  scale_fill_manual(values=palette)+
  #theme(axis.title.x = element_blank())+
  theme(legend.position = "bottom")
p1
```

<img src="../figures/16S-fig-unnamed-chunk-8-2.png" width="672" />

``` r
#dev.off()


n <- 53
palette <- distinctColorPalette(n)

#pdf("barchart_Family.pdf",width=11)
p1=plot_bar(ps_ra, "colony" ,fill="Family")+
  facet_grid(fraction~type,scales="free",space="free")+
  geom_bar(aes(fill=Family), stat="identity",position="stack")+
  theme_bw()+
  theme(strip.text=element_text(face="bold", size=12))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  scale_fill_manual(values=palette)+
  #theme(axis.title.x = element_blank())+
  theme(legend.position = "bottom")
p1
```

<img src="../figures/16S-fig-unnamed-chunk-8-3.png" width="672" />

``` r
#dev.off()

n <- 61
palette <- distinctColorPalette(n)

#pdf("barchart_Genus.pdf",width=11)
p1=plot_bar(ps_ra, "colony" ,fill="Genus")+
  facet_grid(fraction~type,scales="free",space="free")+
  geom_bar(aes(fill=Genus), stat="identity",position="stack")+
  theme_bw()+
  theme(strip.text=element_text(face="bold", size=12))+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  scale_fill_manual(values=palette)+
  #theme(axis.title.x = element_blank())+
  theme(legend.position = "bottom")
p1
```

<img src="../figures/16S-fig-unnamed-chunk-8-4.png" width="672" />

``` r
#dev.off()
```


Differential abundance testing using Corncob



``` r
# use an ASV table that has been filtered to remove low-abundance ASVs, no blanks

otu <- read.table("silva_nochloronomito_otu_table_ps5.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps5.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps5.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 112 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 112 taxa by 6 taxonomic ranks ]
```

``` r
# DA determined relative to first reference level. The default in R is alphabetical order, but I want to compare taxa "after" relative to "before", so I need to manually state the order of input.
sample_data(ps)$type<-factor(sample_data(ps)$type,levels=c("before","after"))
# You can verify the change by checking:
levels(sample_data(ps)$type)
```

```
## [1] "before" "after"
```

``` r
# I will also control for nucleic acid fraction since they are so different in composition

set.seed(1)
treatment.da <- differentialTest(formula = ~ type + fraction, 
                                 phi.formula = ~ type,
                                 formula_null = ~ fraction,
                                 phi.formula_null = ~ 1,
                                 test = "Wald", boot = FALSE,
                                 data = ps,
                                 fdr_cutoff = 0.05)
```

```
## Registered S3 methods overwritten by 'registry':
##   method               from 
##   print.registry_field proxy
##   print.registry_entry proxy
```

``` r
summary(treatment.da)
```

```
##                      Length Class    Mode     
## p                    112    -none-   numeric  
## p_fdr                112    -none-   numeric  
## significant_taxa       1    -none-   character
## significant_models     1    -none-   list     
## all_models           112    -none-   list     
## restrictions_DA        1    -none-   character
## restrictions_DV        1    -none-   character
## discriminant_taxa_DA  23    -none-   character
## discriminant_taxa_DV  17    -none-   character
## data                   1    phyloseq S4       
## full_output            0    -none-   NULL
```

``` r
treatment.da
```

```
## Object of class differentialTest 
## 
## $p: p-values 
## $p_fdr: FDR-adjusted p-values 
## $significant_taxa: taxa names of the statistically significant taxa 
## $significant_models: model summaries of the statistically significant taxa 
## $all_models: all model summaries 
## $restrictions_DA: covariates tested for differential abundance 
## $restrictions_DV: covariates tested for differential variability 
## $discriminant_taxa_DA: taxa for which at least one covariate associated with the abundance was perfectly discriminant 
## $discriminant_taxa_DV: taxa for which at least one covariate associated with the dispersion was perfectly discriminant 
## 
## plot( ) to see a plot of tested coefficients from significant taxa
```

``` r
treatment.da$significant_taxa
```

```
## [1] "TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATCCAAAACTGGCAAGCTAGAGTACGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGG"
```

``` r
plot(treatment.da, c("Family", "Genus"))
```

<img src="../figures/16S-fig-unnamed-chunk-9-1.png" width="672" />

``` r
# results: 1 significant ASV classified as Pseudomonas - that was lower abundance after treatment
```



Differential abundance testing using ANCOM-BC
https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html



``` r
# Start by creating phyloseq object using full dataset (low abundance ASVs are not removed).

otu <- read.table("silva_nochloronomito_otu_table_ps1.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps1.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps1.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps1 <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 6137 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 6137 taxa by 6 taxonomic ranks ]
```

``` r
# prune empty rows
ps_qf <- prune_taxa(taxa_sums(ps1) > 1, ps1) 
ps_qf
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1855 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 1855 taxa by 6 taxonomic ranks ]
```

``` r
# ANCOM-BC2 results will be expressed as a relative relationship - it will test variables in the order that they are input. The default in R is alphabetical order, but I want to compare taxa "after" relative to "before", so I need to manually state the order of input.
sample_data(ps_qf)$type<-factor(sample_data(ps_qf)$type,levels=c("before","after"))
# You can verify the change by checking:
levels(sample_data(ps_qf)$type)
```

```
## [1] "before" "after"
```

``` r
set.seed(123)
output = ancombc2(data = ps_qf, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "fraction + type", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "type", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, mdfdr_control = NULL, 
                  trend_control = NULL)
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## Sensitivity analysis for pseudo-count addition to 0s: ...
```

```
## ANCOM-BC2 primary results ...
```

```
## Warning in pt(abs(W), df = dof, lower.tail = FALSE): NaNs produced
```

``` r
res_prim = output$res
print(res_prim)
```

```
##                                                  taxon lfc_(Intercept)
## 1                                          Marinifilum    -0.911401175
## 2                                               Vibrio    -0.316358710
## 3                                            Clade III     1.167734948
## 4                                             Clade Ia     0.917826089
## 5                                       Proteobacteria    -0.057114454
## 6                                     Rhodobacteraceae    -0.826949816
## 7                              Candidatus Actinomarina     0.604779901
## 8                                         Sphingopyxis              NA
## 9                                        Tenacibaculum              NA
## 10                                      Cryomorphaceae     0.181741611
## 11                                Synechococcus CC9902     0.165289448
## 12                                          Litoricola     0.277426422
## 13                             AEGEAN-169 marine group     0.591986054
## 14                                             Delftia              NA
## 15                                   Blastocatellaceae     0.483404102
## 16                                         Pseudomonas    -2.017416238
## 17                                            Clade II     0.784631455
## 18                                              HIMB11     0.432035655
## 19                                       Roseibacillus    -0.878341754
## 20                                          Rubritalea    -0.071460338
## 21                                        SAR116 clade     0.790350702
## 22                                         SAR86 clade     0.467197726
## 23                                   Halodesulfovibrio    -0.872775340
## 24                                    NS5 marine group     0.939629783
## 25                                   Campylobacterales     0.172319973
## 26                                       Thalassotalea    -0.107854364
## 27                                    OM60(NOR5) clade     0.332015723
## 28                                    Enterobacterales    -0.563556044
## 29                                Pseudoteredinibacter    -0.285285489
## 30                                   Flavobacteriaceae     0.009830924
## 31                                             S25-593     0.605082100
## 32                                    NS4 marine group     0.315987100
## 33                                NS11-12 marine group     0.408995947
## 34                                     Aestuariibacter    -0.629740200
## 35                  Urania-1B-19 marine sediment group     0.407534749
## 36                                     Sandaracinaceae              NA
## 37                         Candidatus Puniceispirillum    -0.536105198
## 38                                    NS9 marine group     0.184342102
## 39                                            Bacteria    -0.103445031
## 40                                   Cyclobacteriaceae     0.075985889
## 41                                     Pseudomonadales     0.362258851
## 42                                             MBMPE27     0.637314753
## 43                                      Propionigenium     0.123467875
## 44                                         Fokiniaceae     0.299140919
## 45                                            Balneola     0.458231849
## 46                                         KI89A clade              NA
## 47                              Synechococcus PCC-7336    -0.406113397
## 48                                   Pseudoalteromonas    -0.908964435
## 49                       Marinimicrobia (SAR406 clade)     0.454348706
## 50                                 Alphaproteobacteria    -0.053958811
## 51                                      Rhodopirellula     0.040055080
## 52  Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium              NA
## 53                                     Francisellaceae     0.943668176
## 54                                             Woeseia    -0.092091194
## 55                                           PS1 clade     0.349745736
## 56                                      Bradymonadales    -0.210045595
## 57                                    Terasakiellaceae    -0.569845781
## 58                                      Staphylococcus    -0.494063935
## 59                                         Azotobacter    -0.358615935
## 60                                     Pedosphaeraceae     0.082362689
## 61                                      Photobacterium    -0.811109962
## 62                                         Alteromonas     0.345041717
## 63                             Candidatus Amoebophilus    -0.617406694
## 64                                    Coraliomargarita     0.052166613
## 65                                      Saprospiraceae     0.236921833
## 66                                       Cohaesibacter    -0.355566956
## 67                           Candidatus Cryptoprodotis     0.453305864
## 68                                       Rickettsiales     1.043292333
## 69                                      Comamonadaceae    -0.332490453
## 70                                        Lentisphaera    -0.043264116
## 71                                             Ekhidna    -0.917538802
## 72                                       Acinetobacter    -0.463577395
## 73                                 Xenococcus PCC-7305    -1.427985844
## 74                                      Cyanobacteriia              NA
## 75                                    Desulfocapsaceae     0.132559444
## 76                       Marine Methylotrophic Group 3     0.532484999
## 77                                     Marine Group II     0.304176986
## 78                                     Actinomarinales              NA
## 79                                            R76-B128     0.244656382
## 80                                     SCGC AAA164-E04    -0.329647613
## 81                                 Gammaproteobacteria    -0.139613538
## 82                                     Chitinophagales    -0.112629722
## 83                                               NB1-j    -1.165431086
## 84                                          Lawsonella    -1.159434855
## 85                                     Hyphomonadaceae    -0.297094713
## 86                                     Corynebacterium     2.562967360
## 87                                         Alcanivorax    -0.650083558
## 88                                      Endozoicomonas    -0.171803620
## 89                                        Chromatiales    -0.805688520
## 90                  [Caedibacter] taeniospiralis group    -0.029510756
## 91                                     JGI 0000069-P22     0.026049545
## 92                                              BD7-11    -0.186385306
## 93                                          OM27 clade     0.144620254
## 94                                       Bacteroidales    -1.075195739
## 95                                         Stappiaceae    -0.883197715
## 96                                  Candidatus Megaira     1.015396921
## 97                                              Rothia              NA
## 98                                    Phormidesmiaceae    -0.801818229
## 99                                      Oligoflexaceae     0.068464239
## 100                                       Pla4 lineage    -0.716583431
## 101                                      Enhydrobacter    -0.857592466
## 102                                     Microcystaceae              NA
## 103                                       Crocinitomix    -0.451828217
## 104                                      Pirellulaceae    -0.633797097
## 105                                        Bacteroidia    -0.643120339
## 106                                         Halieaceae     0.174900150
## 107                                      Halarcobacter              NA
## 108                                 Cyanobium PCC-6307              NA
## 109                                         Lentimonas     0.027329697
## 110                                              OM190    -0.425812733
## 111                                   Unknown Family_1    -0.731667042
## 112                                               PB19    -0.754742996
## 113                                   Cyanobacteriales     0.324657350
## 114                                          Lewinella    -0.072976521
## 115                     Methylobacterium-Methylorubrum    -0.979723245
## 116                                      Marinoscillum              NA
## 117                                     Psychrosphaera              NA
## 118                                 Bacteriovoracaceae    -0.125267756
## 119                                       FS142-36B-02    -0.550501953
## 120                       SAR324 clade(Marine group B)    -0.683673097
## 121                                         Aureispira    -0.942243976
## 122                           BD2-11 terrestrial group     1.911853965
## 123                                        OM182 clade              NA
## 124                                        Portibacter              NA
## 125                                    Cloacibacterium     0.371695330
## 126                                      Lentisphaeria    -0.693911377
## 127                                   Latescibacterota    -0.443144416
## 128                                        Sphingobium              NA
## 129                                        Turneriella    -0.435094595
## 130                                           PRD18C08    -0.891056919
## 131                                    Rhodothermaceae              NA
## 132                                     Niveispirillum              NA
## 133                                           Coxiella    -0.102225853
## 134                                  Rubinisphaeraceae              NA
## 135                                     Flavobacterium     0.534632469
## 136                                     Desulfosarcina    -0.429977949
## 137                                            BIrii41              NA
## 138                                   Cellvibrionaceae    -0.608548580
## 139                                       Acholeplasma              NA
## 140                        Clostridium sensu stricto 1    -0.127085957
## 141                                       Cytophagales    -0.862163462
## 142                                       Chlamydiales    -0.686990308
## 143                                          vadinHA49    -0.512011064
## 144                                      Spirochaeta 2              NA
## 145                                        BD1-7 clade    -0.791749792
## 146                                      Lactobacillus    -2.418653949
## 147                             Sva0081 sediment group    -0.071670542
## 148                          Candidatus Nitrosopumilus              NA
## 149                                       Pla3 lineage    -0.607380882
## 150                                   Flammeovirgaceae    -0.487363384
## 151                                   Chryseobacterium              NA
## 152                                    Lachnospiraceae     0.091272128
## 153                                                A4b    -0.660734278
## 154                                   Flavobacteriales    -0.235126386
## 155                                           WCHB1-41    -0.351027300
## 156                                   Ardenticatenales     0.039720496
## 157                                      SUP05 cluster    -1.487734541
## 158                                             DEV007    -0.330320731
## 159                                   Micavibrionaceae              NA
## 160                                     SCGC AAA011-D5              NA
## 161                                    Blastopirellula    -0.091246218
## 162                                    Amoebophilaceae    -0.245820101
## 163                                 Verrucomicrobiales              NA
## 164                                 Desulfobacteraceae    -0.812618479
## 165                                    Woesearchaeales     0.031720453
## 166                                    Spirochaetaceae    -0.190714265
## 167                                             SM2D12              NA
## 168                                Sediminispirochaeta    -0.555728888
## 169                                       Peredibacter     0.291497150
## 170                                        SS1-B-02-17    -0.473579163
## 171                                   Rhodospirillales    -0.659876365
## 172                                            P3OB-42              NA
## 173                                Bacteroidetes BD2-2    -0.432215911
## 174                        Absconditabacteriales (SR1)     0.103944910
## 175                                     Rickettsiaceae              NA
## 176                                   Margulisbacteria              NA
## 177                                 Caenarcaniphilales              NA
## 178                                   Defluviicoccales              NA
## 179                             Candidatus Omnitrophus              NA
## 180                                       Bdellovibrio    -0.712736320
## 181                                Gastranaerophilales    -0.170708162
##     lfc_fractionRNA lfc_typeafter se_(Intercept) se_fractionRNA se_typeafter
## 1        0.64064282   0.900731300     0.09700300      0.1934736    0.2057049
## 2        0.55961275   0.478963445     0.14694567      0.2382844    0.2403677
## 3       -1.77802050  -0.291574062     0.16259804      0.2700302    0.2778836
## 4       -1.09458225  -0.412717374     0.15967115      0.2796867    0.2763201
## 5        1.25635331  -0.878549394     0.34639017      0.4689974    0.4483086
## 6        0.60880438   2.075004174     0.09908828      0.1991959    0.2092184
## 7       -1.70434187   0.813630029     0.23705089      0.3695771    0.3577894
## 8                NA            NA             NA             NA           NA
## 9                NA            NA             NA             NA           NA
## 10      -0.17470620  -0.039115962     0.29154409      0.3748509    0.3433016
## 11       0.20187022  -0.037490452     0.20324805      0.3571738    0.3631015
## 12       0.07379387  -0.174599485     0.14736958      0.2977012    0.2898291
## 13      -1.14048428   0.452914710     0.15522109      0.3369019    0.3123730
## 14               NA            NA             NA             NA           NA
## 15      -1.26448785   0.764551670     0.40923196      0.4830141    0.4680264
## 16       3.57650709   0.199683206     0.14989420      0.2966020    0.3073263
## 17      -1.13367924  -0.354308267     0.16739990      0.2777617    0.2725745
## 18      -0.21924622  -0.171875592     0.13053528      0.2187623    0.2134869
## 19       2.03373830   1.187148781     0.20850968      0.2331639    0.2404466
## 20      -0.19263560   0.411182886     0.22684537      0.3800690    0.4385086
## 21      -0.48764535  -0.717178310     0.23594037      0.3420143    0.3219115
## 22      -0.82305244   0.334483583     0.14295375      0.2405450    0.2400629
## 23       1.13264448   0.578751747     0.10770083      0.2454972    0.2603051
## 24      -0.94207502  -0.619758052     0.19343968      0.3356714    0.3448674
## 25       0.40004265  -0.523986275     0.22933772      0.3841853    0.3665964
## 26       0.22511294   0.651779502     0.11722367      0.2330234    0.2560130
## 27      -0.14106240  -0.213602777     0.17647007      0.3395688    0.3180682
## 28      -2.43402582   2.627345627     0.12553130      0.1984684    0.2114401
## 29       1.04990221   0.254928554     0.11298087      0.2518459    0.2873095
## 30      -0.25366581   0.702066618     0.18571902      0.4147493    0.4434452
## 31      -0.07567918  -0.996722585     0.16820338      0.2761095    0.2873112
## 32      -0.29726076  -0.102350304     0.19073992      0.2601608    0.2611963
## 33      -1.00129424   0.299125002     0.12769753      0.2514591    0.2443614
## 34       1.14692746  -1.200958122     0.09749224      0.1912364    0.2054954
## 35      -0.88796005   0.088418574     0.12478004      0.2174133    0.2332603
## 36               NA            NA             NA             NA           NA
## 37       0.63800904  -0.322485040     0.08923967      0.1983592    0.2092790
## 38       0.02273643  -0.078392331     0.24879102      0.3332245    0.3229255
## 39       1.29181407  -0.774059091     0.29708271      0.4338981    0.4025304
## 40       0.35608640  -0.289743787     0.20237557      0.2698750    0.2810175
## 41       0.31696620  -1.299306149     0.20312854      0.2823624    0.2750905
## 42      -1.07804369  -0.164789959     0.26789868      0.3310333    0.3173378
## 43      -0.20351833  -0.408034351     0.28786943      0.3395379    0.3185462
## 44      -0.59340379  -0.115234252     0.35069356      0.3942046    0.3743292
## 45      -0.44306498  -0.058041035     0.16228793      0.2588491    0.2644318
## 46               NA            NA             NA             NA           NA
## 47      -0.38831695   0.742701733     0.08634368      0.1906050    0.2061415
## 48       0.18724621   1.319652708     0.14787120      0.2727195    0.2755876
## 49      -0.57303311  -0.127795300     0.14045025      0.2597387    0.2655595
## 50       0.22484351  -0.033632419     0.26504766      0.3114981    0.3121002
## 51      -0.83397927   0.530284382     0.19370258      0.2334961    0.2421001
## 52               NA            NA             NA             NA           NA
## 53      -1.27687988  -0.667258596     0.22495728      0.2608407    0.2607947
## 54      -0.06319596   0.753296876     0.12116302      0.2050653    0.2191235
## 55      -0.14013153  -0.576659153     0.17964661      0.2847281    0.2859925
## 56       0.07285945   0.780407065     0.19190976      0.2427098    0.2532505
## 57       0.45848527  -0.004429774     0.08721549      0.1899100    0.2050654
## 58       0.74896796   0.710838092     0.23531118      0.3175656    0.3309135
## 59       1.04247985  -0.094363776     0.09747473      0.2333145    0.2341650
## 60       0.49278745  -0.311443361     0.13051141      0.2390495    0.2547576
## 61      -1.25153028   2.206036706     0.09503647      0.2016778    0.2179781
## 62      -0.30630838   0.167365944     0.09520067      0.2215238    0.2295761
## 63       0.22066580   1.191573885     0.16942152      0.3788315    0.3942214
## 64      -0.54192597   0.519955630     0.17427673      0.2469054    0.2627314
## 65       0.03698000  -0.266943828     0.13850523      0.2284960    0.2405491
## 66       0.88286124   0.240505971     0.11097281      0.2364666    0.2400570
## 67      -0.13222993  -1.018196195     0.25049579      0.2568565    0.2593851
## 68       0.05727522  -2.023240244     0.14831529      0.3370131    0.3185072
## 69       1.47075356  -1.745992894     0.15315974      0.3563971    0.3414937
## 70       0.20915845   0.593574108     0.13241594      0.2157355    0.2329483
## 71      -0.21327378   4.658957616     0.08598511      0.1879910    0.2022562
## 72       0.98514699   0.052530217     0.22314716      0.4157165    0.3796086
## 73       1.86117502   1.888244163     0.09180373      0.2126424    0.2198816
## 74               NA            NA             NA             NA           NA
## 75      -0.06556510  -0.124246863     0.11324919      0.2096735    0.2274095
## 76       0.47972444  -0.396480335     0.20114656      0.2963711    0.2888332
## 77      -0.70660429  -0.136174457     0.24696349      0.3193958    0.3266289
## 78               NA            NA             NA             NA           NA
## 79      -0.07452302  -0.118052526     0.16060033      0.2530788    0.2739652
## 80       0.20206075  -0.218374649     0.10629642      0.1976208    0.2146277
## 81       0.39702545   0.168386814     0.14519024      0.3134338    0.2998790
## 82       0.23240193   0.231740586     0.22390909      0.3134382    0.2991253
## 83       1.61772369   2.449176040     0.11938647      0.2105112    0.2182733
## 84       1.83985853   0.843803617     0.09128574      0.2019600    0.2185774
## 85       0.83539958   0.105429965     0.14461103      0.2088110    0.2211005
## 86      -1.07144615  -1.861710339     0.08594885      0.1896983    0.2051487
## 87       0.59405728   0.265971010     0.10472592      0.1971562    0.2082730
## 88       0.49637862   0.320306482     0.22348032      0.2676469    0.2697039
## 89       0.79909105   0.122864353     0.11804121      0.1955841    0.2077610
## 90       0.52837000  -0.982976469     0.19469236      0.2678505    0.2607408
## 91      -0.06858957  -0.184792642     0.29295857      0.3705419    0.3411796
## 92       0.32042491  -0.053674585     0.23550371      0.2520809    0.2531442
## 93      -0.25515655  -0.504993641     0.18076180      0.2308146    0.2432010
## 94       0.79096439   1.241926842     0.09158016      0.2026368    0.2205695
## 95       1.03337784   0.965851923     0.08401151      0.1839223    0.1980934
## 96      -0.69060413  -0.888423383     0.11478106      0.2139991    0.2268004
## 97               NA            NA             NA             NA           NA
## 98       3.63362854  -2.985730194     0.27062583      0.2854786    0.2925674
## 99       0.08232153   0.108585346     0.15156078      0.2465250    0.2653614
## 100      0.85227201   0.927098158     0.08632830      0.1895678    0.2049251
## 101      0.90865520   0.972482907     0.09699905      0.2143602    0.2188893
## 102              NA            NA             NA             NA           NA
## 103      0.14920739   0.447305053     0.11472910      0.1922473    0.2054721
## 104     -0.79974940   0.821366083     0.12236299      0.2078235    0.2258358
## 105      0.38991690   1.492475368     0.12617752      0.2888644    0.2901119
## 106      0.25097353  -0.902990393     0.08649397      0.1891625    0.2045052
## 107              NA            NA             NA             NA           NA
## 108              NA            NA             NA             NA           NA
## 109     -1.21267237  -0.157163615     0.09008058      0.1896485    0.2034549
## 110      0.48694426   0.532349470     0.14105426      0.2447339    0.2478587
## 111      0.79338177   1.338412089     0.09727189      0.2112217    0.2184401
## 112      1.00175461   1.179840754     0.11118691      0.2237703    0.2455758
## 113      0.18026072            NA     0.08679269      0.1916362           NA
## 114     -0.12699114   0.603656017     0.12535829      0.2701518    0.2599528
## 115      2.16326130   0.270799865     0.08634368      0.1906050    0.2061415
## 116              NA            NA             NA             NA           NA
## 117              NA            NA             NA             NA           NA
## 118     -0.02361473   0.203550221     0.15473172      0.2119190    0.2256220
## 119      0.36795353   0.575877832     0.09223218      0.2117428    0.2189476
## 120      0.83121515   0.570036517     0.16172208      0.2328914    0.2462884
## 121      1.28707388   1.218052427     0.10621696      0.2152744    0.2401603
## 122      1.38427941  -2.382412807     0.08564453      0.1889995    0.2043898
## 123              NA            NA             NA             NA           NA
## 124              NA            NA             NA             NA           NA
## 125      0.67526122  -1.124326313     0.09068982      0.2094995    0.2170771
## 126     -0.21840804   0.887292837     0.09943441      0.1907886    0.2052940
## 127      0.19365931   0.753165940     0.12909168      0.2320214    0.2462248
## 128              NA            NA             NA             NA           NA
## 129      0.64085300  -0.231729498     0.09352177      0.1897731    0.2044109
## 130      0.85618333   1.085535732     0.08634368      0.1906050    0.2061415
## 131              NA            NA             NA             NA           NA
## 132              NA            NA             NA             NA           NA
## 133     -0.09515887   0.043994277     0.10490482      0.1970927    0.2092298
## 134              NA            NA             NA             NA           NA
## 135      0.25221410            NA     0.09806539      0.2175473           NA
## 136      0.86557336   0.189065714     0.08634368      0.1906050    0.2061415
## 137              NA            NA             NA             NA           NA
## 138      0.90346283   0.257535075     0.14590456      0.2047141    0.2181753
## 139              NA            NA             NA             NA           NA
## 140      0.68958802   1.193966902     0.09274632      0.1905893    0.2056332
## 141      1.55724002            NA     0.17644398      0.2201492           NA
## 142      1.45063492   0.340884200     0.08594885      0.1896983    0.2051487
## 143     -0.68506798   1.209287936     0.11807819      0.2034854    0.2144084
## 144              NA            NA             NA             NA           NA
## 145      1.95496263  -0.753562447     0.09047344      0.1940268    0.2096924
## 146      2.56410446   1.023917920     0.08672163      0.1932500    0.2070920
## 147      0.46144346            NA     0.08885063      0.1963638           NA
## 148              NA            NA             NA             NA           NA
## 149      1.23251125   0.918690117     0.15625494      0.2087378    0.2193016
## 150      0.94450218   0.891562482     0.08589877      0.1900203    0.2050228
## 151              NA            NA             NA             NA           NA
## 152     -0.20384518  -0.551195604     0.13542401      0.2031575    0.2156987
## 153      2.20366077   0.245305802     0.09476602      0.1942587    0.2113685
## 154     -0.22126336   1.912373264     0.08842380      0.2007260    0.2112364
## 155      0.76226786   0.144074788     0.12070558      0.2062544    0.2197078
## 156      0.32198592   0.170208513     0.11126406      0.1972260    0.2115145
## 157      0.82713272   1.691144227     0.08847482      0.1955004    0.2131813
## 158     -0.46860637   0.677647690     0.12820984      0.2047409    0.2163959
## 159              NA            NA             NA             NA           NA
## 160              NA            NA             NA             NA           NA
## 161     -0.76673217  -0.248817963     0.08559595      0.1874748    0.2026693
## 162      0.44053739   0.195487772     0.13612664      0.2158774    0.2271983
## 163              NA            NA             NA             NA           NA
## 164      0.20786330   2.043785774     0.08848717      0.1896040    0.2048538
## 165     -0.31065912   0.377687564     0.16665146      0.2278167    0.2354665
## 166      0.26060754   0.395213821     0.12335905      0.2001533    0.2124049
## 167              NA            NA             NA             NA           NA
## 168      0.55380911   0.863754549     0.08862776      0.1958518    0.2136541
## 169      0.44039595  -0.016818838     0.09042793      0.1931485    0.2083639
## 170      0.90376109   0.227481589     0.11112236      0.1954400    0.2094918
## 171     -0.14219059   1.284013832     0.15442288      0.2139272    0.2267530
## 172              NA            NA             NA             NA           NA
## 173      0.98458638   0.331017099     0.12112618      0.2028343    0.2150900
## 174      0.58315032   0.616362444     0.10119913      0.1927050    0.2073118
## 175              NA            NA             NA             NA           NA
## 176              NA            NA             NA             NA           NA
## 177              NA            NA             NA             NA           NA
## 178              NA            NA             NA             NA           NA
## 179              NA            NA             NA             NA           NA
## 180      1.22218074            NA     0.13832321      0.2041662           NA
## 181      0.92972391   0.253098993     0.08583520      0.1890195    0.2043923
##     W_(Intercept) W_fractionRNA  W_typeafter p_(Intercept) p_fractionRNA
## 1     -9.39559803    3.31126772   4.37875438  7.150682e-04  2.961929e-02
## 2     -2.15289581    2.34850761   1.99262848  4.371309e-02  2.924190e-02
## 3      7.18172851   -6.58452486  -1.04926696  2.016720e-07  8.244512e-07
## 4      5.74822754   -3.91360172  -1.49362031  6.347841e-06  6.552849e-04
## 5     -0.16488474    2.67880645  -1.95969785  8.705404e-01  1.371478e-02
## 6     -8.34558653    3.05630985   9.91788635  1.126825e-03  3.779444e-02
## 7      2.55126609   -4.61160043   2.27404727  1.784850e-02  1.225469e-04
## 8              NA            NA           NA  1.000000e+00  1.000000e+00
## 9              NA            NA           NA  1.000000e+00  1.000000e+00
## 10     0.62337608   -0.46606857  -0.11394052  5.408567e-01  6.467569e-01
## 11     0.81324000    0.56518762  -0.10325061  4.240752e-01  5.771919e-01
## 12     1.88252169    0.24787901  -0.60242229  7.194396e-02  8.063365e-01
## 13     3.81382490   -3.38521138   1.44991651  7.979334e-04  2.351823e-03
## 14             NA            NA           NA  1.000000e+00  1.000000e+00
## 15     1.18124719   -2.61791081   1.63356529  2.507100e-01  1.607231e-02
## 16   -13.45893457   12.05826856   0.64974332  5.217466e-09  1.971760e-08
## 17     4.68716802   -4.08148106  -1.29985850  1.126655e-04  4.944419e-04
## 18     3.30972329   -1.00221200  -0.80508730  2.836206e-03  3.258437e-01
## 19    -4.21247463    8.72235471   4.93726581  1.483807e-01  7.266986e-02
## 20    -0.31501784   -0.50684376   0.93768481  7.586433e-01  6.222701e-01
## 21     3.34979010   -1.42580377  -2.22787402  2.568192e-03  1.662929e-01
## 22     3.26817394   -3.42161594   1.39331663  3.254843e-03  2.235518e-03
## 23    -8.10370118    4.61367627   2.22335914  4.640001e-04  5.768345e-03
## 24     4.85748198   -2.80653931  -1.79709090  5.967464e-05  9.779452e-03
## 25     0.75138086    1.04127532  -1.42932745  4.621352e-01  3.115327e-01
## 26    -0.92007322    0.96605276   2.54588447  4.546640e-01  4.359395e-01
## 27     1.88142799   -0.41541622  -0.67156284  7.321307e-02  6.818610e-01
## 28    -4.48936680  -12.26404809  12.42595861  4.620519e-02  6.583052e-03
## 29    -2.52507783    4.16882859   0.88729603  1.275190e-01  5.300676e-02
## 30     0.05293440   -0.61161238   1.58320958  9.583368e-01  5.480434e-01
## 31     3.59732422   -0.27409114  -3.46913974  1.694255e-03  7.866926e-01
## 32     1.65663854   -1.14260400  -0.39185207  1.124592e-01  2.660633e-01
## 33     3.20284941   -3.98193719   1.22410930  4.467068e-03  7.338706e-04
## 34    -6.45938786    5.99743140  -5.84420979  2.313859e-02  2.669342e-02
## 35     3.26602520   -4.08420249   0.37905543  4.064604e-03  6.321433e-04
## 36             NA            NA           NA  1.000000e+00  1.000000e+00
## 37    -6.00747617    3.21643294  -1.54093327  5.382324e-04  1.472607e-02
## 38     0.74095160    0.06823158  -0.24275673  4.669280e-01  9.462468e-01
## 39    -0.34820280    2.97722889  -1.92298316  7.311542e-01  7.185864e-03
## 40     0.37546967    1.31944946  -1.03105242  7.119550e-01  2.045142e-01
## 41     1.78339711    1.12255107  -4.72319520  1.123665e-01  2.941860e-01
## 42     2.37893950   -3.25660204  -0.51928878  3.868239e-02  8.623582e-03
## 43     0.42890235   -0.59939802  -1.28092683  6.793100e-01  5.654920e-01
## 44     0.85299804   -1.50531941  -0.30784197  4.184652e-01  1.706593e-01
## 45     2.82357322   -1.71167257  -0.21949341  9.399405e-03  9.985198e-02
## 46             NA            NA           NA  1.000000e+00  1.000000e+00
## 47    -4.70345257   -2.03728641   3.60287328  1.000000e+00  1.000000e+00
## 48    -6.14700101    0.68658893   4.78850516  2.748981e-04  5.117480e-01
## 49     3.23494416   -2.20619106  -0.48123032  4.358685e-03  3.987993e-02
## 50    -0.20358154    0.72181352  -0.10776161  8.414166e-01  4.815106e-01
## 51     0.20678651   -3.57170473   2.19035173  8.393803e-01  3.411649e-03
## 52             NA            NA           NA  1.000000e+00  1.000000e+00
## 53     4.19487731   -4.89524727  -2.55855927  4.061343e-03  1.762656e-03
## 54    -0.76006026   -0.30817476   3.43777336  4.814997e-01  7.703747e-01
## 55     1.94685409   -0.49215905  -2.01634368  6.732912e-02  6.285556e-01
## 56    -1.09450190    0.30019161   3.08156210  2.899386e-01  7.678962e-01
## 57    -6.53376780    2.41422380  -0.02160176  2.263242e-02  1.371429e-01
## 58    -2.09961949    2.35846722   2.14810844  5.012776e-02  2.986031e-02
## 59    -3.67906553    4.46813167  -0.40297981  1.430753e-02  6.590645e-03
## 60     0.63107652    2.06144511  -1.22250860  5.355091e-01  5.321272e-02
## 61    -8.53472327   -6.20559368  10.12045199  1.345205e-02  2.499808e-02
## 62     3.62436219   -1.38273330   0.72902178  1.103937e-02  2.160083e-01
## 63    -3.64420471    0.58249057   3.02260090  3.859614e-03  5.719827e-01
## 64     0.29933207   -2.19487248   1.97903870  7.683153e-01  4.234952e-02
## 65     1.71056240    0.16184092  -1.10972683  1.179503e-01  8.746544e-01
## 66    -3.20409071    3.73355587   1.00187036  7.574817e-03  2.855339e-03
## 67     1.80963467   -0.51480084  -3.92542214  1.203359e-01  6.250998e-01
## 68     7.03428691    0.16994957  -6.35225827  4.123990e-04  8.706353e-01
## 69    -2.17087366    4.12672748  -5.11281088  8.205249e-02  9.113976e-03
## 70    -0.32672890    0.96951317   2.54809352  7.549698e-01  3.697315e-01
## 71   -10.67090331   -1.13448915  23.03493531  1.759156e-03  3.390395e-01
## 72    -2.07745151    2.36975699   0.13837993  5.322877e-02  2.989729e-02
## 73   -15.55476896    8.75260496   8.58755051  4.087138e-02  7.242087e-02
## 74             NA            NA           NA  1.000000e+00  1.000000e+00
## 75     1.17051119   -0.31270086  -0.54635733  2.945428e-01  7.671297e-01
## 76     2.64724887    1.61866116  -1.37269648  3.307351e-02  1.495527e-01
## 77     1.23166781   -2.21231535  -0.41690882  2.416644e-01  4.708421e-02
## 78             NA            NA           NA  1.000000e+00  1.000000e+00
## 79     1.52338652   -0.29446566  -0.43090340  1.460475e-01  7.719675e-01
## 80    -3.10121084    1.02246719  -1.01745791  2.108387e-02  3.460033e-01
## 81    -0.96159041    1.26669619   0.56151581  3.515013e-01  2.245803e-01
## 82    -0.50301542    0.74146018   0.77472736  6.222586e-01  4.698648e-01
## 83    -9.76183533    7.68474018  11.22068572  2.284063e-03  4.578493e-03
## 84   -12.70116093    9.11001365   3.86043410  5.001977e-02  6.960266e-02
## 85    -2.05444019    4.00074529   0.47684173  8.571355e-02  7.112937e-03
## 86    29.81968157   -5.64815908  -9.07493054  2.134098e-02  1.115567e-01
## 87    -6.20747532    3.01313000   1.27703074  8.063653e-04  2.360509e-02
## 88    -0.76876396    1.85460254   1.18762275  4.671734e-01  1.060484e-01
## 89    -6.82548497    4.08566431   0.59137361  2.409302e-03  1.502844e-02
## 90    -0.15157634    1.97263059  -3.76993756  8.828644e-01  8.000284e-02
## 91     0.08891887   -0.18510610  -0.54162856  9.309020e-01  8.568462e-01
## 92    -0.79143255    1.27111911  -0.21203168  4.645626e-01  2.596098e-01
## 93     0.80005985   -1.10546084  -2.07644546  4.380524e-01  2.889992e-01
## 94   -11.74048787    3.90336094   5.63054683  7.176832e-03  5.980570e-02
## 95   -10.51281760    5.61855758   4.87574040  1.344315e-04  2.471883e-03
## 96     8.84638034   -3.22713634  -3.91720395  9.015071e-04  3.206253e-02
## 97             NA            NA           NA  1.000000e+00  1.000000e+00
## 98    -2.96282961   12.72820049 -10.20527458  5.940885e-02  1.046170e-03
## 99     0.45172792    0.33392769   0.40919804  6.575283e-01  7.427708e-01
## 100   -8.30067845    4.49586928   4.52408217  7.632707e-02  1.393328e-01
## 101   -8.84124641    4.23891780   4.44280758  3.075940e-04  8.178074e-03
## 102            NA            NA           NA  1.000000e+00  1.000000e+00
## 103   -3.93821820    0.77612207   2.17696254  1.698520e-02  4.810113e-01
## 104   -5.17964689   -3.84821404   3.63700526  1.396886e-02  3.097751e-02
## 105   -5.09694865    1.34982663   5.14448192  1.404266e-03  2.190905e-01
## 106    2.02210798    1.32676122  -4.41548949  2.923771e-01  4.111765e-01
## 107            NA            NA           NA  1.000000e+00  1.000000e+00
## 108            NA            NA           NA  1.000000e+00  1.000000e+00
## 109    0.30339167   -6.39431556  -0.77247408  7.814147e-01  7.746638e-03
## 110   -3.01878670    1.98968866   2.14779441  1.168101e-02  7.206518e-02
## 111   -7.52187564    3.75615572   6.12713657  4.869976e-03  3.297679e-02
## 112   -6.78805651    4.47670844   4.80438535  8.012599e-05  1.539706e-03
## 113    3.74060709    0.94064019           NA  1.663026e-01  5.194667e-01
## 114   -0.58214356   -0.47007325   2.32217558  6.013320e-01  6.703636e-01
## 115  -11.34678600   11.34944749   1.31366005  1.000000e+00  1.000000e+00
## 116            NA            NA           NA  1.000000e+00  1.000000e+00
## 117            NA            NA           NA  1.000000e+00  1.000000e+00
## 118   -0.80958031   -0.11143279   0.90217383  4.448251e-01  9.144011e-01
## 119   -5.96865406    1.73773815   2.63020878  9.410387e-03  1.806483e-01
## 120   -4.22745673    3.56911004   2.31450821  3.900396e-03  9.107753e-03
## 121   -8.87093696    5.97875828   5.07183026  4.685432e-05  5.538823e-04
## 122   22.32312947    7.32424863 -11.65622359  2.849934e-02  8.638533e-02
## 123            NA            NA           NA  1.000000e+00  1.000000e+00
## 124            NA            NA           NA  1.000000e+00  1.000000e+00
## 125    4.09853440    3.22321114  -5.17938614  5.469285e-02  8.426661e-02
## 126   -6.97858420   -1.14476443   4.32205879  1.992206e-02  3.708267e-01
## 127   -3.43278825    0.83466125   3.05885491  1.804587e-01  5.572173e-01
## 128            NA            NA           NA  1.000000e+00  1.000000e+00
## 129   -4.65233507    3.37694241  -1.13364527  4.322796e-02  7.761835e-02
## 130  -10.31988597    4.49192513   5.26597355  1.000000e+00  1.000000e+00
## 131            NA            NA           NA  1.000000e+00  1.000000e+00
## 132            NA            NA           NA  1.000000e+00  1.000000e+00
## 133   -0.97446287   -0.48281284   0.21026778  3.622903e-01  6.439620e-01
## 134            NA            NA           NA  1.000000e+00  1.000000e+00
## 135    5.45179594    1.15935277           NA  1.212272e-02  3.302064e-01
## 136   -4.97984282    4.54118947   0.91716470  1.000000e+00  1.000000e+00
## 137            NA            NA           NA  1.000000e+00  1.000000e+00
## 138   -4.17086744    4.41328972   1.18040461  3.118391e-03  2.246518e-03
## 139            NA            NA           NA  1.000000e+00  1.000000e+00
## 140   -1.37025332    3.61818768   5.80629344  4.013519e-01  1.716646e-01
## 141   -4.88632964    7.07356627           NA  1.285114e-01  8.940735e-02
## 142   -7.99301332    7.64706358   1.66164434  7.923534e-02  8.278049e-02
## 143   -4.33620349   -3.36666945   5.64011555  7.455563e-03  1.996208e-02
## 144            NA            NA           NA  1.000000e+00  1.000000e+00
## 145   -8.75118474   10.07573637  -3.59365602  1.280738e-02  9.707039e-03
## 146  -27.88985895   13.26833075   4.94426638  2.281644e-02  4.788986e-02
## 147   -0.80664077    2.34994106           NA  5.678771e-01  2.561313e-01
## 148            NA            NA           NA  1.000000e+00  1.000000e+00
## 149   -3.88711473    5.90459087   4.18916369  1.773392e-02  4.117149e-03
## 150   -5.67369458    4.97053177   4.34860171  1.110648e-01  1.263916e-01
## 151            NA            NA           NA  1.000000e+00  1.000000e+00
## 152    0.67397299   -1.00338493  -2.55539658  6.224566e-01  4.989244e-01
## 153   -6.97226964   11.34394858   1.16055974  6.054614e-03  1.469468e-03
## 154   -2.65908488   -1.10231552   9.05323870  1.171010e-01  3.852356e-01
## 155   -2.90812809    3.69576619   0.65575650  6.209089e-02  3.437873e-02
## 156    0.35699305    1.63257343   0.80471327  7.552455e-01  2.441543e-01
## 157  -16.81534341    4.23084974   7.93289147  3.781492e-02  1.477592e-01
## 158   -2.57640692   -2.28877761   3.13151748  4.964792e-02  7.075735e-02
## 159            NA            NA           NA  1.000000e+00  1.000000e+00
## 160            NA            NA           NA  1.000000e+00  1.000000e+00
## 161   -1.06601089   -4.08978785  -1.22770433  3.980681e-01  5.490810e-02
## 162   -1.80581917    2.04068294   0.86042809  9.414394e-02  6.213113e-02
## 163            NA            NA           NA  1.000000e+00  1.000000e+00
## 164   -9.18346098    1.09630252   9.97680152  6.905036e-02  4.707747e-01
## 165    0.19034008   -1.36363606   1.60399726  8.522247e-01  1.977165e-01
## 166   -1.54600955    1.30203938   1.86066207  2.621421e-01  3.226722e-01
## 167            NA            NA           NA  1.000000e+00  1.000000e+00
## 168   -6.27037067    2.82769546   4.04277110  1.006804e-01  2.163987e-01
## 169    3.22352995    2.28009024  -0.08071858  1.914987e-01  2.631250e-01
## 170   -4.26178115    4.62423700   1.08587327  1.467245e-01  1.355825e-01
## 171   -4.27317739   -0.66466820   5.66261064  1.291747e-02  5.426224e-01
## 172            NA            NA           NA  1.000000e+00  1.000000e+00
## 173   -3.56831140    4.85414037   1.53897050  2.341110e-02  8.314461e-03
## 174    1.02713248    3.02612994   2.97311768  4.914796e-01  2.031822e-01
## 175            NA            NA           NA  1.000000e+00  1.000000e+00
## 176            NA            NA           NA  1.000000e+00  1.000000e+00
## 177            NA            NA           NA  1.000000e+00  1.000000e+00
## 178            NA            NA           NA  1.000000e+00  1.000000e+00
## 179            NA            NA           NA  1.000000e+00  1.000000e+00
## 180   -5.15268770    5.98620438           NA  3.566192e-02  2.678965e-02
## 181   -1.98878980    4.91866629   1.23830029  2.966010e-01  1.276890e-01
##      p_typeafter q_(Intercept) q_fractionRNA q_typeafter diff_(Intercept)
## 1   0.0118864660  1.201315e-01  1.000000e+00  1.00000000            FALSE
## 2   0.0601274605  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 3   0.3045070652  3.630096e-05  1.484012e-04  1.00000000             TRUE
## 4   0.1483065408  1.136264e-03  1.146749e-01  1.00000000             TRUE
## 5   0.0628207746  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 6   0.0005802286  1.847994e-01  1.000000e+00  0.10444114            FALSE
## 7   0.0326118075  1.000000e+00  2.193590e-02  1.00000000            FALSE
## 8   1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 9   1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 10  0.9105462206  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 11  0.9186220039  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 12  0.5525429963  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 13  0.1595123758  1.332549e-01  3.927544e-01  1.00000000            FALSE
## 14  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 15  0.1172543165  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 16  0.5271727973  9.443613e-07  3.568885e-06  1.00000000             TRUE
## 17  0.2071047043  1.971646e-02  8.801065e-02  1.00000000             TRUE
## 18  0.4283602005  4.452844e-01  1.000000e+00  1.00000000            FALSE
## 19  0.1272207415  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 20  0.3685384245  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 21  0.0351211354  4.057743e-01  1.000000e+00  1.00000000            FALSE
## 22  0.1762952027  5.045006e-01  3.778025e-01  1.00000000            FALSE
## 23  0.0768032495  7.888002e-02  9.287035e-01  1.00000000            FALSE
## 24  0.0849177536  1.056241e-02  1.000000e+00  1.00000000             TRUE
## 25  0.1700365862  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 26  0.1258184256  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 27  0.5088530781  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 28  0.0064142514  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 29  0.4685319913  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 30  0.1298773268  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 31  0.0022932161  2.744693e-01  1.000000e+00  0.39213996            FALSE
## 32  0.6991157459  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 33  0.2351359505  6.655932e-01  1.276935e-01  1.00000000            FALSE
## 34  0.0280523202  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 35  0.7088495914  6.173241e-01  1.112572e-01  1.00000000            FALSE
## 36  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 37  0.1672326428  9.096127e-02  1.000000e+00  1.00000000            FALSE
## 38  0.8105491713  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 39  0.0681483524  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 40  0.3169628736  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 41  0.0014959885  1.000000e+00  1.000000e+00  0.26030200            FALSE
## 42  0.6148573997  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 43  0.2361025770  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 44  0.7660648862  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 45  0.8281218138  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 46  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 47  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 48  0.0013756107  4.755737e-02  1.000000e+00  0.24210749             TRUE
## 49  0.6358483655  6.538027e-01  1.000000e+00  1.00000000            FALSE
## 50  0.9156129915  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 51  0.0473283404  1.000000e+00  5.595104e-01  1.00000000            FALSE
## 52  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 53  0.0376306236  6.173241e-01  2.996514e-01  1.00000000            FALSE
## 54  0.0184797353  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 55  0.0589369127  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 56  0.0071506930  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 57  0.9847270286  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 58  0.0455669207  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 59  0.7036162224  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 60  0.2364643716  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 61  0.0096226820  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 62  0.4934510127  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 63  0.0116016731  5.943806e-01  1.000000e+00  1.00000000            FALSE
## 64  0.0642438565  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 65  0.2930893030  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 66  0.3361811629  1.000000e+00  4.711309e-01  1.00000000            FALSE
## 67  0.0077534597  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 68  0.0007135444  7.052023e-02  1.000000e+00  0.12772445            FALSE
## 69  0.0037296449  1.000000e+00  1.000000e+00  0.63030998            FALSE
## 70  0.0436025042  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 71  0.0001792136  2.832242e-01  1.000000e+00  0.03243767            FALSE
## 72  0.8915661056  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 73  0.0738004890  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 74  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 75  0.6083061250  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 76  0.2122074788  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 77  0.6841061889  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 78  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 79  0.6719525948  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 80  0.3481943786  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 81  0.5827388822  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 82  0.4505412410  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 83  0.0015175116  3.654501e-01  7.417159e-01  0.26252951            FALSE
## 84  0.1613620863  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 85  0.6503453575  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 86  0.0698695854  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 87  0.2487826616  1.338566e-01  1.000000e+00  1.00000000            FALSE
## 88  0.2737132476  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 89  0.5860724454  3.830790e-01  1.000000e+00  1.00000000            FALSE
## 90  0.0044168832  1.000000e+00  1.000000e+00  0.74203637            FALSE
## 91  0.5999332708  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 92  0.8404546824  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 93  0.0582466731  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 94  0.0301246456  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 95  0.0045699340  2.339109e-02  4.103325e-01  0.76317898             TRUE
## 96  0.0172884005  1.487487e-01  1.000000e+00  1.00000000            FALSE
## 97  1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 98  0.0020053253  1.000000e+00  1.809875e-01  0.34491595            FALSE
## 99  0.6878197157  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 100 0.1384911679  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 101 0.0067472409  5.290617e-02  1.000000e+00  1.00000000            FALSE
## 102 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 103 0.0950667800  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 104 0.0358168959  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 105 0.0013321137  2.288953e-01  1.000000e+00  0.23578412            FALSE
## 106 0.1417869111  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 107 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 108 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 109 0.4961098938  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 110 0.0548569607  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 111 0.0087406634  7.207564e-01  1.000000e+00  1.00000000            FALSE
## 112 0.0009677876  1.410217e-02  2.632897e-01  0.17226619             TRUE
## 113 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 114 0.1028795884  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 115 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 116 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 117 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 118 0.3969444233  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 119 0.0783110855  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 120 0.0538281208  5.967605e-01  1.000000e+00  1.00000000            FALSE
## 121 0.0014441399  8.340070e-03  9.803716e-02  0.25272447             TRUE
## 122 0.0544828919  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 123 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 124 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 125 0.0353144087  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 126 0.0495847749  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 127 0.2011509935  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 128 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 129 0.3745400983  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 130 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 131 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 132 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 133 0.8394480566  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 134 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 135 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 136 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 137 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 138 0.2717440080  4.864690e-01  3.778025e-01  1.00000000            FALSE
## 139 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 140 0.1085778712  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 141 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 142 0.3448899663  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 143 0.0024304573  1.000000e+00  1.000000e+00  0.41317774            FALSE
## 144 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 145 0.0694619387  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 146 0.1270453577  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 147 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 148 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 149 0.0138158255  1.000000e+00  6.710952e-01  1.00000000            FALSE
## 150 0.1438948016  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 151 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 152 0.2374647982  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 153 0.3297834779  8.900283e-01  2.527484e-01  1.00000000            FALSE
## 154 0.0119820580  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 155 0.5587839195  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 156 0.5054411527  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 157 0.0798295931  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 158 0.0259130277  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 159 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 160 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 161 0.3444437922  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 162 0.4051425985  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 163 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 164 0.0635975954  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 165 0.1346942632  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 166 0.2038605601  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 167 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 168 0.1543725211  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 169 0.9487241238  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 170 0.4738058075  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 171 0.0047949518  1.000000e+00  1.000000e+00  0.79596200            FALSE
## 172 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 173 0.1986436126  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 174 0.2065580564  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 175 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 176 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 177 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 178 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 179 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 180 1.0000000000  1.000000e+00  1.000000e+00  1.00000000            FALSE
## 181 0.4324767296  1.000000e+00  1.000000e+00  1.00000000            FALSE
##     diff_fractionRNA diff_typeafter passed_ss_(Intercept) passed_ss_fractionRNA
## 1              FALSE          FALSE                  TRUE                 FALSE
## 2              FALSE          FALSE                  TRUE                 FALSE
## 3               TRUE          FALSE                  TRUE                  TRUE
## 4              FALSE          FALSE                 FALSE                 FALSE
## 5              FALSE          FALSE                  TRUE                 FALSE
## 6              FALSE          FALSE                  TRUE                  TRUE
## 7               TRUE          FALSE                 FALSE                  TRUE
## 8              FALSE          FALSE                 FALSE                 FALSE
## 9              FALSE          FALSE                  TRUE                  TRUE
## 10             FALSE          FALSE                 FALSE                 FALSE
## 11             FALSE          FALSE                  TRUE                  TRUE
## 12             FALSE          FALSE                  TRUE                  TRUE
## 13             FALSE          FALSE                 FALSE                 FALSE
## 14             FALSE          FALSE                 FALSE                 FALSE
## 15             FALSE          FALSE                  TRUE                  TRUE
## 16              TRUE          FALSE                  TRUE                  TRUE
## 17             FALSE          FALSE                  TRUE                 FALSE
## 18             FALSE          FALSE                 FALSE                  TRUE
## 19             FALSE          FALSE                  TRUE                  TRUE
## 20             FALSE          FALSE                  TRUE                  TRUE
## 21             FALSE          FALSE                 FALSE                  TRUE
## 22             FALSE          FALSE                 FALSE                 FALSE
## 23             FALSE          FALSE                 FALSE                 FALSE
## 24             FALSE          FALSE                 FALSE                 FALSE
## 25             FALSE          FALSE                  TRUE                  TRUE
## 26             FALSE          FALSE                  TRUE                  TRUE
## 27             FALSE          FALSE                  TRUE                  TRUE
## 28             FALSE          FALSE                  TRUE                  TRUE
## 29             FALSE          FALSE                  TRUE                  TRUE
## 30             FALSE          FALSE                  TRUE                  TRUE
## 31             FALSE          FALSE                 FALSE                  TRUE
## 32             FALSE          FALSE                 FALSE                 FALSE
## 33             FALSE          FALSE                 FALSE                 FALSE
## 34             FALSE          FALSE                  TRUE                  TRUE
## 35             FALSE          FALSE                 FALSE                 FALSE
## 36             FALSE          FALSE                  TRUE                  TRUE
## 37             FALSE          FALSE                  TRUE                  TRUE
## 38             FALSE          FALSE                  TRUE                  TRUE
## 39             FALSE          FALSE                  TRUE                  TRUE
## 40             FALSE          FALSE                  TRUE                  TRUE
## 41             FALSE          FALSE                  TRUE                  TRUE
## 42             FALSE          FALSE                  TRUE                  TRUE
## 43             FALSE          FALSE                 FALSE                  TRUE
## 44             FALSE          FALSE                  TRUE                  TRUE
## 45             FALSE          FALSE                 FALSE                 FALSE
## 46             FALSE          FALSE                  TRUE                  TRUE
## 47             FALSE          FALSE                  TRUE                  TRUE
## 48             FALSE          FALSE                 FALSE                  TRUE
## 49             FALSE          FALSE                 FALSE                 FALSE
## 50             FALSE          FALSE                  TRUE                  TRUE
## 51             FALSE          FALSE                 FALSE                 FALSE
## 52             FALSE          FALSE                 FALSE                 FALSE
## 53             FALSE          FALSE                  TRUE                  TRUE
## 54             FALSE          FALSE                  TRUE                  TRUE
## 55             FALSE          FALSE                 FALSE                  TRUE
## 56             FALSE          FALSE                  TRUE                 FALSE
## 57             FALSE          FALSE                  TRUE                  TRUE
## 58             FALSE          FALSE                  TRUE                  TRUE
## 59             FALSE          FALSE                  TRUE                 FALSE
## 60             FALSE          FALSE                  TRUE                  TRUE
## 61             FALSE          FALSE                  TRUE                  TRUE
## 62             FALSE          FALSE                  TRUE                  TRUE
## 63             FALSE          FALSE                  TRUE                  TRUE
## 64             FALSE          FALSE                  TRUE                 FALSE
## 65             FALSE          FALSE                  TRUE                  TRUE
## 66             FALSE          FALSE                  TRUE                  TRUE
## 67             FALSE          FALSE                  TRUE                  TRUE
## 68             FALSE          FALSE                  TRUE                  TRUE
## 69             FALSE          FALSE                  TRUE                 FALSE
## 70             FALSE          FALSE                  TRUE                  TRUE
## 71             FALSE           TRUE                  TRUE                  TRUE
## 72             FALSE          FALSE                  TRUE                 FALSE
## 73             FALSE          FALSE                  TRUE                  TRUE
## 74             FALSE          FALSE                  TRUE                  TRUE
## 75             FALSE          FALSE                  TRUE                  TRUE
## 76             FALSE          FALSE                  TRUE                  TRUE
## 77             FALSE          FALSE                  TRUE                  TRUE
## 78             FALSE          FALSE                  TRUE                  TRUE
## 79             FALSE          FALSE                 FALSE                 FALSE
## 80             FALSE          FALSE                  TRUE                  TRUE
## 81             FALSE          FALSE                  TRUE                  TRUE
## 82             FALSE          FALSE                  TRUE                  TRUE
## 83             FALSE          FALSE                  TRUE                  TRUE
## 84             FALSE          FALSE                  TRUE                  TRUE
## 85             FALSE          FALSE                  TRUE                  TRUE
## 86             FALSE          FALSE                  TRUE                  TRUE
## 87             FALSE          FALSE                  TRUE                  TRUE
## 88             FALSE          FALSE                  TRUE                  TRUE
## 89             FALSE          FALSE                  TRUE                  TRUE
## 90             FALSE          FALSE                  TRUE                  TRUE
## 91             FALSE          FALSE                 FALSE                  TRUE
## 92             FALSE          FALSE                  TRUE                  TRUE
## 93             FALSE          FALSE                 FALSE                 FALSE
## 94             FALSE          FALSE                  TRUE                  TRUE
## 95             FALSE          FALSE                 FALSE                  TRUE
## 96             FALSE          FALSE                  TRUE                  TRUE
## 97             FALSE          FALSE                  TRUE                  TRUE
## 98             FALSE          FALSE                  TRUE                  TRUE
## 99             FALSE          FALSE                  TRUE                 FALSE
## 100            FALSE          FALSE                  TRUE                  TRUE
## 101            FALSE          FALSE                  TRUE                  TRUE
## 102            FALSE          FALSE                  TRUE                  TRUE
## 103            FALSE          FALSE                  TRUE                  TRUE
## 104            FALSE          FALSE                  TRUE                  TRUE
## 105            FALSE          FALSE                  TRUE                  TRUE
## 106            FALSE          FALSE                  TRUE                  TRUE
## 107            FALSE          FALSE                 FALSE                 FALSE
## 108            FALSE          FALSE                  TRUE                 FALSE
## 109            FALSE          FALSE                  TRUE                  TRUE
## 110            FALSE          FALSE                  TRUE                  TRUE
## 111            FALSE          FALSE                  TRUE                  TRUE
## 112            FALSE          FALSE                 FALSE                  TRUE
## 113            FALSE          FALSE                  TRUE                  TRUE
## 114            FALSE          FALSE                  TRUE                  TRUE
## 115            FALSE          FALSE                  TRUE                  TRUE
## 116            FALSE          FALSE                  TRUE                  TRUE
## 117            FALSE          FALSE                  TRUE                 FALSE
## 118            FALSE          FALSE                  TRUE                 FALSE
## 119            FALSE          FALSE                  TRUE                  TRUE
## 120            FALSE          FALSE                  TRUE                  TRUE
## 121            FALSE          FALSE                 FALSE                  TRUE
## 122            FALSE          FALSE                  TRUE                  TRUE
## 123            FALSE          FALSE                  TRUE                 FALSE
## 124            FALSE          FALSE                  TRUE                 FALSE
## 125            FALSE          FALSE                  TRUE                  TRUE
## 126            FALSE          FALSE                  TRUE                  TRUE
## 127            FALSE          FALSE                  TRUE                  TRUE
## 128            FALSE          FALSE                 FALSE                  TRUE
## 129            FALSE          FALSE                  TRUE                  TRUE
## 130            FALSE          FALSE                  TRUE                  TRUE
## 131            FALSE          FALSE                 FALSE                  TRUE
## 132            FALSE          FALSE                  TRUE                 FALSE
## 133            FALSE          FALSE                  TRUE                  TRUE
## 134            FALSE          FALSE                  TRUE                  TRUE
## 135            FALSE          FALSE                  TRUE                  TRUE
## 136            FALSE          FALSE                  TRUE                  TRUE
## 137            FALSE          FALSE                  TRUE                  TRUE
## 138            FALSE          FALSE                  TRUE                  TRUE
## 139            FALSE          FALSE                  TRUE                  TRUE
## 140            FALSE          FALSE                  TRUE                  TRUE
## 141            FALSE          FALSE                  TRUE                  TRUE
## 142            FALSE          FALSE                  TRUE                  TRUE
## 143            FALSE          FALSE                  TRUE                  TRUE
## 144            FALSE          FALSE                 FALSE                 FALSE
## 145            FALSE          FALSE                  TRUE                  TRUE
## 146            FALSE          FALSE                  TRUE                 FALSE
## 147            FALSE          FALSE                  TRUE                  TRUE
## 148            FALSE          FALSE                  TRUE                 FALSE
## 149            FALSE          FALSE                  TRUE                  TRUE
## 150            FALSE          FALSE                  TRUE                  TRUE
## 151            FALSE          FALSE                  TRUE                  TRUE
## 152            FALSE          FALSE                  TRUE                  TRUE
## 153            FALSE          FALSE                  TRUE                  TRUE
## 154            FALSE          FALSE                  TRUE                  TRUE
## 155            FALSE          FALSE                  TRUE                  TRUE
## 156            FALSE          FALSE                  TRUE                  TRUE
## 157            FALSE          FALSE                  TRUE                  TRUE
## 158            FALSE          FALSE                  TRUE                  TRUE
## 159            FALSE          FALSE                  TRUE                  TRUE
## 160            FALSE          FALSE                 FALSE                 FALSE
## 161            FALSE          FALSE                  TRUE                  TRUE
## 162            FALSE          FALSE                  TRUE                  TRUE
## 163            FALSE          FALSE                  TRUE                  TRUE
## 164            FALSE          FALSE                  TRUE                  TRUE
## 165            FALSE          FALSE                  TRUE                  TRUE
## 166            FALSE          FALSE                  TRUE                  TRUE
## 167            FALSE          FALSE                  TRUE                  TRUE
## 168            FALSE          FALSE                  TRUE                  TRUE
## 169            FALSE          FALSE                  TRUE                  TRUE
## 170            FALSE          FALSE                  TRUE                  TRUE
## 171            FALSE          FALSE                  TRUE                  TRUE
## 172            FALSE          FALSE                  TRUE                  TRUE
## 173            FALSE          FALSE                  TRUE                  TRUE
## 174            FALSE          FALSE                  TRUE                  TRUE
## 175            FALSE          FALSE                  TRUE                  TRUE
## 176            FALSE          FALSE                  TRUE                  TRUE
## 177            FALSE          FALSE                  TRUE                  TRUE
## 178            FALSE          FALSE                 FALSE                  TRUE
## 179            FALSE          FALSE                 FALSE                 FALSE
## 180            FALSE          FALSE                  TRUE                  TRUE
## 181            FALSE          FALSE                  TRUE                  TRUE
##     passed_ss_typeafter
## 1                  TRUE
## 2                  TRUE
## 3                  TRUE
## 4                  TRUE
## 5                 FALSE
## 6                  TRUE
## 7                  TRUE
## 8                  TRUE
## 9                  TRUE
## 10                 TRUE
## 11                 TRUE
## 12                 TRUE
## 13                 TRUE
## 14                 TRUE
## 15                 TRUE
## 16                 TRUE
## 17                 TRUE
## 18                 TRUE
## 19                 TRUE
## 20                 TRUE
## 21                FALSE
## 22                 TRUE
## 23                FALSE
## 24                 TRUE
## 25                FALSE
## 26                 TRUE
## 27                 TRUE
## 28                 TRUE
## 29                 TRUE
## 30                 TRUE
## 31                 TRUE
## 32                 TRUE
## 33                 TRUE
## 34                 TRUE
## 35                 TRUE
## 36                 TRUE
## 37                 TRUE
## 38                 TRUE
## 39                FALSE
## 40                 TRUE
## 41                FALSE
## 42                 TRUE
## 43                FALSE
## 44                 TRUE
## 45                 TRUE
## 46                 TRUE
## 47                 TRUE
## 48                FALSE
## 49                 TRUE
## 50                 TRUE
## 51                 TRUE
## 52                 TRUE
## 53                 TRUE
## 54                 TRUE
## 55                FALSE
## 56                 TRUE
## 57                 TRUE
## 58                FALSE
## 59                 TRUE
## 60                 TRUE
## 61                 TRUE
## 62                FALSE
## 63                 TRUE
## 64                 TRUE
## 65                 TRUE
## 66                 TRUE
## 67                 TRUE
## 68                 TRUE
## 69                 TRUE
## 70                 TRUE
## 71                FALSE
## 72                 TRUE
## 73                 TRUE
## 74                 TRUE
## 75                 TRUE
## 76                 TRUE
## 77                 TRUE
## 78                 TRUE
## 79                 TRUE
## 80                 TRUE
## 81                 TRUE
## 82                 TRUE
## 83                 TRUE
## 84                FALSE
## 85                 TRUE
## 86                 TRUE
## 87                 TRUE
## 88                 TRUE
## 89                 TRUE
## 90                 TRUE
## 91                FALSE
## 92                 TRUE
## 93                FALSE
## 94                FALSE
## 95                 TRUE
## 96                 TRUE
## 97                 TRUE
## 98                 TRUE
## 99                 TRUE
## 100                TRUE
## 101                TRUE
## 102                TRUE
## 103                TRUE
## 104               FALSE
## 105                TRUE
## 106                TRUE
## 107               FALSE
## 108                TRUE
## 109                TRUE
## 110                TRUE
## 111                TRUE
## 112                TRUE
## 113                TRUE
## 114                TRUE
## 115                TRUE
## 116               FALSE
## 117               FALSE
## 118                TRUE
## 119                TRUE
## 120                TRUE
## 121                TRUE
## 122                TRUE
## 123                TRUE
## 124                TRUE
## 125                TRUE
## 126                TRUE
## 127                TRUE
## 128               FALSE
## 129                TRUE
## 130                TRUE
## 131                TRUE
## 132                TRUE
## 133                TRUE
## 134                TRUE
## 135               FALSE
## 136                TRUE
## 137               FALSE
## 138                TRUE
## 139                TRUE
## 140                TRUE
## 141                TRUE
## 142                TRUE
## 143                TRUE
## 144                TRUE
## 145                TRUE
## 146                TRUE
## 147                TRUE
## 148                TRUE
## 149                TRUE
## 150                TRUE
## 151                TRUE
## 152                TRUE
## 153                TRUE
## 154                TRUE
## 155                TRUE
## 156                TRUE
## 157                TRUE
## 158                TRUE
## 159                TRUE
## 160                TRUE
## 161                TRUE
## 162                TRUE
## 163                TRUE
## 164                TRUE
## 165                TRUE
## 166                TRUE
## 167               FALSE
## 168               FALSE
## 169                TRUE
## 170                TRUE
## 171                TRUE
## 172               FALSE
## 173                TRUE
## 174                TRUE
## 175                TRUE
## 176               FALSE
## 177                TRUE
## 178                TRUE
## 179                TRUE
## 180                TRUE
## 181                TRUE
```

``` r
write.table(res_prim,"ANCOM-BC2_all.txt",sep="\t",col.names=NA)
# nothing significantly different after treatment relative to before (column "diff_typeafter"=TRUE AND column "passed_ss_typeafter"=TRUE); One ASV was likely a false positive (column "diff_typeafter"=TRUE AND column "passed_ss_typeafter"=FALSE)
```


Differential abundance testing using ANCOM-BC: DNA fraction only
https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html



``` r
# Start by creating phyloseq object using full dataset (low abundance ASVs are not removed).
otu <- read.table("silva_nochloronomito_otu_table_ps1.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps1.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps1.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps1 <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 6137 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 6137 taxa by 6 taxonomic ranks ]
```

``` r
# prune empty rows
ps_qf <- prune_taxa(taxa_sums(ps1) > 1, ps1) 
ps_qf
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1855 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 1855 taxa by 6 taxonomic ranks ]
```

``` r
# Subset and test DNA fraction
ps_dna <- subset_samples(ps_qf, fraction == "DNA")
ps_dna 
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1855 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 1855 taxa by 6 taxonomic ranks ]
```

``` r
# prune empty rows
ps_dna_qf <- prune_taxa(taxa_sums(ps_dna) > 1, ps_dna) 
ps_dna_qf
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1171 taxa and 18 samples ]
## sample_data() Sample Data:       [ 18 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 1171 taxa by 6 taxonomic ranks ]
```

``` r
sample_data(ps_dna_qf)$type<-factor(sample_data(ps_dna_qf)$type,levels=c("before","after"))
# You can verify the change by checking:
levels(sample_data(ps_dna_qf)$type)
```

```
## [1] "before" "after"
```

``` r
set.seed(123)
output = ancombc2(data = ps_dna_qf, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "type", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "type", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, mdfdr_control = NULL, 
                  trend_control = NULL)
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## Sensitivity analysis for pseudo-count addition to 0s: ...
```

```
## ANCOM-BC2 primary results ...
```

```
## Warning in pt(abs(W), df = dof, lower.tail = FALSE): NaNs produced
```

``` r
res_prim_dna = output$res
print(res_prim_dna)
```

```
##                                  taxon lfc_(Intercept) lfc_typeafter
## 1                          Marinifilum              NA            NA
## 2                               Vibrio    -0.136915288   0.327867404
## 3                            Clade III     0.470962013  -0.879641127
## 4                             Clade Ia     0.393520872  -0.686038275
## 5                       Proteobacteria     0.385181293  -0.934610115
## 6                     Rhodobacteraceae    -0.851637820   1.831259492
## 7              Candidatus Actinomarina     0.221991744  -0.257215453
## 8                        Tenacibaculum    -1.024984340   2.212537413
## 9                       Cryomorphaceae     0.102427457  -0.100548829
## 10                Synechococcus CC9902    -0.023762618   0.357170452
## 11                          Litoricola     0.212677049  -0.233928716
## 12             AEGEAN-169 marine group     0.085777421   0.083320354
## 13                   Blastocatellaceae    -0.192564268   0.746423152
## 14                         Pseudomonas    -0.531314871   1.064657559
## 15                            Clade II     0.529359284  -1.025634305
## 16                              HIMB11     0.211227418  -0.230304638
## 17                       Roseibacillus    -0.552107631   0.978481975
## 18                          Rubritalea    -0.364189001   0.676568252
## 19                        SAR116 clade     0.317172295  -0.495166832
## 20                         SAR86 clade     0.139925837  -0.052050686
## 21                   Halodesulfovibrio              NA            NA
## 22                    NS5 marine group     0.530405764  -1.028250505
## 23                   Campylobacterales     0.364625086  -0.919362485
## 24                       Thalassotalea    -0.829718135   1.365802536
## 25                    OM60(NOR5) clade     0.064296583   0.012495136
## 26                    Enterobacterales    -0.927324658   2.398611351
## 27                Pseudoteredinibacter     0.155728282   0.163496975
## 28                   Flavobacteriaceae    -0.229013002   0.569550057
## 29                             S25-593     0.465118210  -1.116082936
## 30                    NS4 marine group     0.259947881  -0.352105796
## 31                NS11-12 marine group     0.105705348   0.033500535
## 32  Urania-1B-19 marine sediment group     0.078626717  -0.014313699
## 33                     Sandaracinaceae     0.688550728  -1.065233784
## 34         Candidatus Puniceispirillum    -0.269577639  -0.581228166
## 35                    NS9 marine group     0.070529188  -0.151466155
## 36                            Bacteria     0.675923542  -1.689564417
## 37                   Cyclobacteriaceae     0.219099848  -0.539933161
## 38                     Pseudomonadales     0.663133920  -2.111096593
## 39                             MBMPE27     0.339082983  -0.567116019
## 40                          Fabibacter              NA            NA
## 41                      Propionigenium     0.068423310  -0.560346628
## 42                         Fokiniaceae     0.248113663  -0.512313427
## 43                            Balneola     0.297847551  -0.446854971
## 44                   Pseudoalteromonas    -2.019902013   2.380751945
## 45       Marinimicrobia (SAR406 clade)     0.269204042  -0.465622268
## 46                 Alphaproteobacteria    -0.062072046  -0.015845751
## 47                      Rhodopirellula     0.015758479   0.096669219
## 48                     Francisellaceae     0.209440183  -0.440415558
## 49                             Woeseia    -0.469179541   0.975124887
## 50                           PS1 clade     0.327447978  -1.010674360
## 51                      Bradymonadales    -0.272532286   0.927523952
## 52                      Staphylococcus    -0.421201602   1.215617544
## 53                     Pedosphaeraceae     0.194974832  -0.256537428
## 54                      Photobacterium    -0.727466172   1.933674699
## 55                         Alteromonas              NA            NA
## 56             Candidatus Amoebophilus    -0.670346299   1.161149228
## 57                    Coraliomargarita    -0.203080835   0.730854581
## 58                      Saprospiraceae     0.122801121  -0.259326699
## 59                       Cohaesibacter    -0.018049846   0.044739077
## 60                            Serratia              NA            NA
## 61           Candidatus Cryptoprodotis     0.542027279  -1.198111865
## 62                       Rickettsiales     1.207250717  -1.789290444
## 63                      Comamonadaceae              NA            NA
## 64                        Lentisphaera     0.083793124   0.239536596
## 65                              DEV008              NA            NA
## 66                    Desulfocapsaceae    -0.118924497  -0.211405037
## 67       Marine Methylotrophic Group 3     0.947174388  -1.033345317
## 68                     Geminicoccaceae     0.383465322  -0.496937400
## 69                     Marine Group II     0.346406425  -1.057649158
## 70                     Actinomarinales              NA            NA
## 71                            R76-B128     0.411302805  -0.678051555
## 72                     SCGC AAA164-E04    -0.318376546  -0.254702826
## 73                 Gammaproteobacteria    -0.211880861   0.626648430
## 74                               NB1-j    -0.909472331   2.193095749
## 75                          Lawsonella    -0.371970106   0.760068069
## 76                     Chitinophagales     0.062894489  -0.377630521
## 77                         Alcanivorax    -0.252406457  -0.300666018
## 78                        Chromatiales              NA            NA
## 79  [Caedibacter] taeniospiralis group     0.136753374  -1.616589471
## 80                     JGI 0000069-P22     0.019342143  -0.353028796
## 81                             Tistlia              NA            NA
## 82                              BD7-11    -0.104714624  -0.240258215
## 83                          OM27 clade     0.213349449  -0.658436560
## 84                       Bacteroidales    -0.569171272   0.981730485
## 85                         Stappiaceae    -0.517332027   0.738624948
## 86                  Candidatus Megaira     0.986515818  -1.499883350
## 87                              Rothia    -1.027055855   2.183637252
## 88                    Phormidesmiaceae              NA            NA
## 89              Synechococcus PCC-7336    -0.280079164   0.647766130
## 90                      Oligoflexaceae     0.252450943  -0.337798561
## 91           Bacteroidetes VC2.1 Bac22              NA            NA
## 92                       Enhydrobacter    -0.739027411   0.933827902
## 93                        Crocinitomix    -0.400845573   0.283442363
## 94                       Pirellulaceae    -0.590689315   0.601287952
## 95                             Ekhidna              NA            NA
## 96                         Bacteroidia    -0.349373379   0.801985822
## 97                          Halieaceae     0.370191929  -1.056037205
## 98                    Reichenbachiella    -1.434852059   1.514856952
## 99                       Acinetobacter    -0.228952826   0.252835963
## 100                         Lentimonas    -0.169242069  -0.395955378
## 101                              OM190    -0.270669618   0.316236414
## 102                               PB19    -0.572064367   0.990060800
## 103                       Luminiphilus              NA            NA
## 104                           Gven-F17    -0.082754116   0.507643677
## 105     Methylobacterium-Methylorubrum    -0.205057638  -0.055135566
## 106                        KI89A clade              NA            NA
## 107                   Draconibacterium              NA            NA
## 108                 Bacteriovoracaceae    -0.090277859   0.090682365
## 109       SAR324 clade(Marine group B)    -0.183809965  -0.339507211
## 110                          Lewinella    -0.260078496   0.397366490
## 111                         Aureispira    -0.587601141   0.987972091
## 112                     Kiloniellaceae              NA            NA
## 113                     Endozoicomonas     0.073307698  -0.391744348
## 114                  Thalassobaculales     0.625885489  -1.108310925
## 115                        Portibacter    -0.314137103   0.668201384
## 116                    Cloacibacterium              NA            NA
## 117                   NS7 marine group              NA            NA
## 118                      Lentisphaeria    -0.555365082   0.688050233
## 119                   Terasakiellaceae    -0.204791649  -0.510455923
## 120                   Latescibacterota     0.811460455  -0.987415084
## 121                        Turneriella     0.119375546  -0.403059332
## 122                           PRD18C08    -1.063883386   0.772919607
## 123                    Rhodothermaceae              NA            NA
## 124                         Nitrospira              NA            NA
## 125                        Gimesiaceae     0.535446540  -0.678555125
## 126                           Coxiella    -0.237508284  -0.337223076
## 127                           Massilia    -0.641694311  -0.071458543
## 128                  Rubinisphaeraceae              NA            NA
## 129                 Thiohalorhabdaceae    -1.074193030   1.291992912
## 130                       Roseimarinus              NA            NA
## 131                            BIrii41    -1.140958726   1.150565263
## 132                    Hyphomonadaceae    -0.132198806   0.009432929
## 133                    Arcobacteraceae              NA            NA
## 134                   Cellvibrionaceae    -0.165629023   0.025728574
## 135                       Acholeplasma              NA            NA
## 136        Clostridium sensu stricto 1    -0.070547076   1.021784276
## 137                  Sedimenticolaceae              NA            NA
## 138                          Weissella              NA            NA
## 139                 Puniceispirillales              NA            NA
## 140                       Cytophagales              NA            NA
## 141                            CL500-3              NA            NA
## 142           BD2-11 terrestrial group     1.871736374  -2.392267797
## 143                      Spirochaeta 2    -0.026739739   0.404121170
## 144          Candidatus Nitrosopumilus     0.183872141  -0.429325479
## 145                       Pla3 lineage    -0.274842945   0.946927974
## 146                        BD1-7 clade              NA            NA
## 147                       Pla4 lineage    -0.704326308   0.490437868
## 148                   Chryseobacterium              NA            NA
## 149                    Lachnospiraceae    -0.031654162  -0.748426720
## 150                  Helicobacteraceae     1.707340189  -2.779174553
## 151                             LCP-80              NA            NA
## 152                          P.palmC41              NA            NA
## 153                   Flavobacteriales              NA            NA
## 154                     Nannocystaceae              NA            NA
## 155                    Micavibrionales    -1.284799762   1.713206377
## 156                      Luteolibacter              NA            NA
## 157                          vadinHA49    -0.678463096   1.050794635
## 158                           WCHB1-41    -0.041541057  -0.692207767
## 159                        Gemmataceae              NA            NA
## 160                        Lishizhenia              NA            NA
## 161                   Ardenticatenales              NA            NA
## 162                      SUP05 cluster    -1.140958726   1.371868449
## 163                             DEV007    -0.444408501   0.437734324
## 164                   Micavibrionaceae              NA            NA
## 165                        Ponticaulis              NA            NA
## 166                     SCGC AAA011-D5    -0.546428825   0.627770524
## 167                      Marinoscillum              NA            NA
## 168            Candidatus Hepatoplasma    -0.523339396   1.382409565
## 169                     Unknown Family    -0.691390332   1.123919021
## 170                    Blastopirellula     0.001734839  -0.388910074
## 171                       Flammeovirga     0.018676908  -0.079468728
## 172                    Calditrichaceae    -0.733361790   1.424608668
## 173                Xenococcus PCC-7305    -1.050899898   1.630575001
## 174                                A4b    -0.349820664   0.005693578
## 175                    Amoebophilaceae     0.023232491  -0.125504255
## 176                           AT-s2-59              NA            NA
## 177                    Woesearchaeales    -0.063799762   0.025817292
## 178                    Spirochaetaceae    -0.293936467   0.141342135
## 179                   Planctomycetales    -0.388561955  -0.152269569
## 180                         Haliangium              NA            NA
## 181                Sediminispirochaeta    -0.460326560   0.718542243
## 182                 Desulfobacteraceae              NA            NA
## 183                        SS1-B-02-17    -0.368607254  -0.025787810
## 184                   Rhodospirillales    -0.522513250   1.011557927
## 185                   Izemoplasmatales              NA            NA
## 186                         Actibacter              NA            NA
## 187                    Aestuariibacter              NA            NA
## 188                              SZB85              NA            NA
## 189                        Subgroup 10              NA            NA
## 190                            SBR1031              NA            NA
## 191                     Flavobacterium              NA            NA
## 192                          Pirellula     0.239624286  -0.209255327
## 193                   Carboxylicivirga              NA            NA
## 194                   Saccharibacillus              NA            NA
## 195                       Pir4 lineage    -0.357752717   0.295744705
## 196                       SAR202 clade     0.811460455  -1.025348894
## 197                            P3OB-42    -0.678860606   0.252890011
## 198                    Kapabacteriales              NA            NA
## 199                     Marinifilaceae              NA            NA
## 200        Absconditabacteriales (SR1)     0.182297143   0.441996973
## 201                     Rickettsiaceae              NA            NA
## 202                     Cyanobacteriia              NA            NA
## 203                          0319-6G20              NA            NA
## 204                       Zixibacteria              NA            NA
## 205                      Desulfobacter              NA            NA
## 206                   Margulisbacteria              NA            NA
## 207                 Caenarcaniphilales              NA            NA
## 208                       FS142-36B-02    -0.158759318   0.200283690
## 209                  Lentimicrobiaceae              NA            NA
## 210                              WPS-2              NA            NA
## 211                    Anaerolineaceae    -0.275216424   0.061327985
## 212                             MidBa8              NA            NA
## 213                   Defluviicoccales              NA            NA
## 214                   Phycisphaeraceae    -0.419057460   0.796000009
## 215                Bacteroidetes BD2-2              NA            NA
## 216                          Sumerlaea              NA            NA
## 217             Candidatus Omnitrophus    -0.134552259  -0.042298064
## 218          Candidatus Falkowbacteria              NA            NA
## 219                    Gracilibacteria              NA            NA
## 220                       Bdellovibrio              NA            NA
## 221                     Thermoplasmata    -0.072483870   0.102852828
## 222            AKAU3564 sediment group    -1.140958726   0.965004097
## 223                     Microcystaceae              NA            NA
## 224                      Desulfovibrio              NA            NA
## 225                             SM2D12    -0.704326308   0.713932844
## 226                            Clade I              NA            NA
## 227                  NS10 marine group              NA            NA
## 228                      Thalassospira              NA            NA
## 229                Gastranaerophilales    -0.068382423   0.124804056
## 230                     Leptospiraceae    -0.477948978   0.466793093
## 231                       Peredibacter     1.014193009  -0.541217678
## 232                    SCGC AAA286-E23    -0.275216424   0.559782002
## 233                              BD2-3              NA            NA
##     se_(Intercept) se_typeafter W_(Intercept) W_typeafter p_(Intercept)
## 1               NA           NA            NA          NA  1.000000e+00
## 2     1.087775e-01    0.2793378 -1.258673e+00  1.17373098  2.398167e-01
## 3     1.547486e-01    0.2680259  3.043400e+00 -3.28192623  9.418831e-03
## 4     1.615531e-01    0.2784175  2.435861e+00 -2.46406338  2.999913e-02
## 5     3.477148e-01    0.4543736  1.107751e+00 -2.05691991  2.916039e-01
## 6     5.783464e-02    0.2810349 -1.472539e+01  6.51612863  4.580095e-03
## 7     2.233229e-01    0.3183107  9.940393e-01 -0.80806423  3.383467e-01
## 8     8.737780e-03    0.2745817 -1.173049e+02  8.05784658  7.266419e-05
## 9     3.303602e-01    0.4101624  3.100478e-01 -0.24514393  7.618431e-01
## 10    1.588499e-01    0.2671691 -1.495916e-01  1.33687031  8.833827e-01
## 11    1.265688e-01    0.2938070  1.680328e+00 -0.79619846  1.167492e-01
## 12    1.299800e-01    0.2659234  6.599277e-01  0.31332466  5.208214e-01
## 13    4.306470e-01    0.5457108 -4.471511e-01  1.36779975  6.634425e-01
## 14    6.588205e-02    0.2969045 -8.064638e+00  3.58585870  1.283975e-03
## 15    1.638770e-01    0.2869574  3.230223e+00 -3.57416937  6.573783e-03
## 16    1.363169e-01    0.2643185  1.549532e+00 -0.87131499  1.452481e-01
## 17    2.061403e-01    0.3522578 -2.678310e+00  2.77774400  2.274902e-01
## 18    1.844791e-01    0.7554494 -1.974148e+00  0.89558376  8.380422e-02
## 19    2.619151e-01    0.3602042  1.210974e+00 -1.37468377  2.474610e-01
## 20    1.459218e-01    0.2833293  9.589100e-01 -0.18371090  3.550996e-01
## 21              NA           NA            NA          NA  1.000000e+00
## 22    1.893554e-01    0.3853611  2.801113e+00 -2.66827822  1.499774e-02
## 23    2.200052e-01    0.3836920  1.657348e+00 -2.39609480  1.318277e-01
## 24    5.267126e-17    0.3695636 -1.575277e+16  3.69571687  4.041320e-17
## 25    1.575708e-01    0.2969502  4.080489e-01  0.04207822  6.904243e-01
## 26    1.116097e-01    0.2981663 -8.308639e+00  8.04454341  1.417838e-02
## 27    5.267126e-17    0.4642724  2.956608e+15  0.35215740  1.143964e-31
## 28    8.343672e-02    0.3400217 -2.744751e+00  1.67504052  2.267067e-02
## 29    1.672751e-01    0.3432034  2.780558e+00 -3.25195813  1.663317e-02
## 30    2.195544e-01    0.3742302  1.183979e+00 -0.94088029  2.576184e-01
## 31    1.201911e-01    0.2549085  8.794776e-01  0.13142180  3.951090e-01
## 32    1.148623e-01    0.3260993  6.845302e-01 -0.04389369  5.066464e-01
## 33    6.523258e-02    0.2846971  1.055532e+01 -3.74163848  1.816347e-03
## 34    2.111661e-02    0.2699563 -1.276614e+01 -2.15304576  1.037005e-03
## 35    2.822781e-01    0.4259119  2.498571e-01 -0.35562794  8.073025e-01
## 36    2.693157e-01    0.3859219  2.509781e+00 -4.37799528  2.899774e-02
## 37    1.848886e-01    0.3840185  1.185037e+00 -1.40600814  2.634005e-01
## 38    2.068410e-01    0.3398803  3.206008e+00 -6.21129456  3.271302e-02
## 39    2.917435e-01    0.4245727  1.162264e+00 -1.33573361  2.892644e-01
## 40              NA           NA            NA          NA  1.000000e+00
## 41    3.154559e-01    0.4153967  2.169029e-01 -1.34894340  8.368566e-01
## 42    3.836631e-01    0.5166334  6.466967e-01 -0.99163829  5.463384e-01
## 43    1.642101e-01    0.3530153  1.813820e+00 -1.26582323  9.284502e-02
## 44    5.267126e-17    0.3270813 -3.834922e+16  7.27877726  3.910219e-50
## 45    1.261152e-01    0.3185216  2.134588e+00 -1.46182316  5.411071e-02
## 46    2.989408e-01    0.4082815 -2.076399e-01 -0.03881085  8.401325e-01
## 47    1.832836e-01    0.3129487  8.597865e-02  0.30889796  9.330284e-01
## 48    2.327881e-01    0.3604082  8.997032e-01 -1.22199090  4.191459e-01
## 49    8.706938e-02    0.2891466 -5.388571e+00  3.37242417  3.275642e-02
## 50    1.872899e-01    0.3945646  1.748349e+00 -2.56149282  1.109770e-01
## 51    1.979421e-01    0.3522014 -1.376829e+00  2.63350448  1.937089e-01
## 52    2.286918e-01    0.4747276 -1.841787e+00  2.56066340  9.863105e-02
## 53    1.165381e-01    0.3355510  1.673056e+00 -0.76452578  1.224843e-01
## 54    1.833676e-02    0.2805849 -3.967256e+01  6.89158407  1.604346e-02
## 55              NA           NA            NA          NA  1.000000e+00
## 56    8.063376e-02    0.5817290 -8.313470e+00  1.99603120  1.642872e-04
## 57    1.767210e-01    0.3454468 -1.149161e+00  2.11567940  2.748635e-01
## 58    9.803526e-02    0.3323175  1.252622e+00 -0.78035826  2.505675e-01
## 59    7.107342e-02    0.2998935 -2.539606e-01  0.14918323  8.052329e-01
## 60              NA           NA            NA          NA  1.000000e+00
## 61    2.679559e-01    0.3881233  2.022822e+00 -3.08693630  9.901322e-02
## 62    8.411507e-02    0.3402465  1.435237e+01 -5.25880675  7.330969e-04
## 63              NA           NA            NA          NA  1.000000e+00
## 64    8.942258e-02    0.3237112  9.370466e-01  0.73997005  4.178608e-01
## 65              NA           NA            NA          NA  1.000000e+00
## 66    6.398916e-02    0.3328590 -1.858510e+00 -0.63511879  2.041981e-01
## 67    1.869111e-01    0.3545638  5.067513e+00 -2.91441252  1.483619e-02
## 68    4.391070e-17    0.2853018  8.732844e+15 -1.74179549  1.000000e+00
## 69    2.715553e-01    0.4545202  1.275639e+00 -2.32695763  2.378747e-01
## 70              NA           NA            NA          NA  1.000000e+00
## 71    1.428262e-01    0.3669040  2.879744e+00 -1.84803517  1.383717e-02
## 72    7.356707e-02    0.3087858 -4.327704e+00 -0.82485287  2.274774e-02
## 73    7.836912e-02    0.2951855 -2.703627e+00  2.12289719  2.692241e-02
## 74    8.249863e-02    0.2876825 -1.102409e+01  7.62332029  8.128197e-03
## 75    4.391070e-17    0.2853018 -8.471058e+15  2.66408432  1.000000e+00
## 76    2.121435e-01    0.3256579  2.964715e-01 -1.15959286  7.744210e-01
## 77    6.516503e-02    0.2773526 -3.873342e+00 -1.08405682  3.045897e-02
## 78              NA           NA            NA          NA  1.000000e+00
## 79    1.767547e-01    0.3165206  7.736900e-01 -5.10737482  4.740871e-01
## 80    3.113733e-01    0.4085929  6.211882e-02 -0.86401110  9.524858e-01
## 81              NA           NA            NA          NA  1.000000e+00
## 82    2.389588e-01    0.3656417 -4.382120e-01 -0.65708656  6.908680e-01
## 83    1.815264e-01    0.3592115  1.175308e+00 -1.83300539  2.670964e-01
## 84    4.391070e-17    0.2965886 -1.296202e+16  3.31007492  4.911425e-17
## 85    3.194788e-02    0.2764646 -1.619300e+01  2.67167971  3.792012e-03
## 86    6.496334e-02    0.3025812  1.518573e+01 -4.95696113  4.308388e-03
## 87    9.987633e-02    0.2936157 -1.028328e+01  7.43705782  9.324580e-03
## 88              NA           NA            NA          NA  1.000000e+00
## 89    4.391070e-17    0.2853018 -6.378380e+15  2.27045926  1.000000e+00
## 90    1.361987e-01    0.3870556  1.853548e+00 -0.87273908  9.079440e-02
## 91              NA           NA            NA          NA  1.000000e+00
## 92    4.085207e-02    0.2777483 -1.809033e+01  3.36213710  3.041739e-03
## 93    8.310497e-02    0.2825564 -4.823365e+00  1.00313542  1.698209e-02
## 94    7.970419e-02    0.3239339 -7.411020e+00  1.85620542  5.082514e-03
## 95              NA           NA            NA          NA  1.000000e+00
## 96    4.161315e-02    0.3604150 -8.395745e+00  2.22517319  3.544425e-03
## 97    9.471734e-03    0.2801035  3.908386e+01 -3.77016846  1.628501e-02
## 98    4.391070e-17    0.2853018 -3.267659e+16  5.30966477  1.000000e+00
## 99    1.943880e-01    0.3790420 -1.177814e+00  0.66703941  2.918725e-01
## 100   4.273823e-02    0.2780595 -3.959969e+00 -1.42399532  5.825350e-02
## 101   1.037574e-01    0.3249896 -2.608677e+00  0.97306625  4.019164e-02
## 102   4.326772e-02    0.3720637 -1.322150e+01  2.66099807  4.422245e-05
## 103             NA           NA            NA          NA  1.000000e+00
## 104   4.531083e-02    0.2839114 -1.826365e+00  1.78803579  3.189142e-01
## 105   4.391070e-17    0.2853018 -4.669878e+15 -0.19325348  1.000000e+00
## 106             NA           NA            NA          NA  1.000000e+00
## 107             NA           NA            NA          NA  1.000000e+00
## 108   1.429581e-01    0.3253224 -6.314989e-01  0.27874613  5.509994e-01
## 109   1.354701e-01    0.3170355 -1.356831e+00 -1.07088078  2.463548e-01
## 110   6.872681e-02    0.2889960 -3.784237e+00  1.37498963  1.644698e-01
## 111   4.599314e-02    0.3614681 -1.277584e+01  2.73322044  5.227119e-05
## 112             NA           NA            NA          NA  1.000000e+00
## 113   2.190302e-01    0.3569380  3.346921e-01 -1.09751373  7.514410e-01
## 114   4.986609e-02    0.2847442  1.255132e+01 -3.89230439  5.061441e-02
## 115   6.170929e-02    0.2730383 -5.090596e+00  2.44728092  7.028532e-03
## 116             NA           NA            NA          NA  1.000000e+00
## 117             NA           NA            NA          NA  1.000000e+00
## 118   4.372183e-02    0.2782271 -1.270224e+01  2.47298042  6.140798e-03
## 119   1.045751e-02    0.2801419 -1.958321e+01 -1.82213370  3.248023e-02
## 120   4.391070e-17    0.2853018  1.847979e+16 -3.46094929  1.000000e+00
## 121   1.945888e-02    0.2806677  6.134760e+00 -1.43607323  1.028678e-01
## 122   4.391070e-17    0.2853018 -2.422834e+16  2.70912973  1.000000e+00
## 123             NA           NA            NA          NA  1.000000e+00
## 124             NA           NA            NA          NA  1.000000e+00
## 125   3.859284e-17    0.2858391  1.387424e+16 -2.37390587  4.588501e-17
## 126   8.314463e-02    0.2825692 -2.856568e+00 -1.19341771  6.475896e-02
## 127   4.391070e-17    0.2853018 -1.461362e+16 -0.25046649  1.000000e+00
## 128             NA           NA            NA          NA  1.000000e+00
## 129   4.391070e-17    0.2853018 -2.446312e+16  4.52851290  1.000000e+00
## 130             NA           NA            NA          NA  1.000000e+00
## 131   4.391070e-17    0.2852965 -2.598361e+16  4.03287560  2.450082e-17
## 132   1.237743e-01    0.2983218 -1.068063e+00  0.03161998  3.638046e-01
## 133             NA           NA            NA          NA  1.000000e+00
## 134   1.172680e-01    0.3098965 -1.412397e+00  0.08302313  2.169404e-01
## 135             NA           NA            NA          NA  1.000000e+00
## 136   4.662030e-02    0.2841429 -1.513226e+00  3.59602276  3.717590e-01
## 137             NA           NA            NA          NA  1.000000e+00
## 138             NA           NA            NA          NA  1.000000e+00
## 139             NA           NA            NA          NA  1.000000e+00
## 140             NA           NA            NA          NA  1.000000e+00
## 141             NA           NA            NA          NA  1.000000e+00
## 142   3.859284e-17    0.2800387  4.849957e+16 -8.54263290  1.312630e-17
## 143   1.164489e-01    0.3082266 -2.296264e-01  1.31111697  8.260099e-01
## 144   2.262416e-01    0.3618401  8.127249e-01 -1.18650612  4.619781e-01
## 145   1.369808e-01    0.3045043 -2.006435e+00  3.10973616  1.384604e-01
## 146             NA           NA            NA          NA  1.000000e+00
## 147   4.391070e-17    0.2853018 -1.603997e+16  1.71901424  1.000000e+00
## 148             NA           NA            NA          NA  1.000000e+00
## 149   1.126944e-01    0.3036010 -2.808850e-01 -2.46516546  8.256749e-01
## 150   4.391070e-17    0.2853018  3.888210e+16 -9.74117405  1.000000e+00
## 151             NA           NA            NA          NA  1.000000e+00
## 152             NA           NA            NA          NA  1.000000e+00
## 153             NA           NA            NA          NA  1.000000e+00
## 154             NA           NA            NA          NA  1.000000e+00
## 155   4.391070e-17    0.2853018 -2.925938e+16  6.00489144  1.000000e+00
## 156             NA           NA            NA          NA  1.000000e+00
## 157   8.138451e-02    0.2918645 -8.336514e+00  3.60028260  3.617980e-03
## 158   7.104358e-02    0.2843273 -5.847264e-01 -2.43454574  6.179079e-01
## 159             NA           NA            NA          NA  1.000000e+00
## 160             NA           NA            NA          NA  1.000000e+00
## 161             NA           NA            NA          NA  1.000000e+00
## 162   4.391070e-17    0.3005411 -2.598361e+16  4.56466146  2.450082e-17
## 163   9.645205e-02    0.3000070 -4.607559e+00  1.45908042  1.922702e-02
## 164             NA           NA            NA          NA  1.000000e+00
## 165             NA           NA            NA          NA  1.000000e+00
## 166   5.369128e-02    0.2746390 -1.017724e+01  2.28580226  2.021570e-03
## 167             NA           NA            NA          NA  1.000000e+00
## 168   4.531083e-02    0.2839114 -1.154999e+01  4.86915901  5.498155e-02
## 169   3.371406e-02    0.2803007 -2.050748e+01  4.00968979  2.369353e-03
## 170   1.384166e-02    0.2748123  1.253346e-01 -1.41518417  9.117210e-01
## 171   4.391070e-17    0.2853018  4.253384e+14 -0.27854267  1.000000e+00
## 172   4.391070e-17    0.2853018 -1.670121e+16  4.99333911  1.000000e+00
## 173   4.391070e-17    0.2853018 -2.393266e+16  5.71526350  1.000000e+00
## 174   4.341104e-02    0.2945017 -8.058334e+00  0.01933293  3.991788e-03
## 175   9.964421e-02    0.3066989  2.331545e-01 -0.40921001  8.208578e-01
## 176             NA           NA            NA          NA  1.000000e+00
## 177   1.687806e-01    0.3007114 -3.780041e-01  0.08585405  7.152601e-01
## 178   9.062158e-02    0.2954864 -3.243559e+00  0.47833720  1.903856e-01
## 179   1.696019e-01    0.3308246 -2.291024e+00 -0.46027279  2.620067e-01
## 180             NA           NA            NA          NA  1.000000e+00
## 181   4.391070e-17    0.3017765 -1.048324e+16  2.38104119  6.072737e-17
## 182             NA           NA            NA          NA  1.000000e+00
## 183   9.201342e-02    0.2959529 -4.006016e+00 -0.08713486  1.557333e-01
## 184   1.368533e-01    0.3175026 -3.818055e+00  3.18598326  3.161493e-02
## 185             NA           NA            NA          NA  1.000000e+00
## 186             NA           NA            NA          NA  1.000000e+00
## 187             NA           NA            NA          NA  1.000000e+00
## 188             NA           NA            NA          NA  1.000000e+00
## 189             NA           NA            NA          NA  1.000000e+00
## 190             NA           NA            NA          NA  1.000000e+00
## 191             NA           NA            NA          NA  1.000000e+00
## 192   4.391070e-17    0.2853018  5.457081e+15 -0.73345251  1.000000e+00
## 193             NA           NA            NA          NA  1.000000e+00
## 194             NA           NA            NA          NA  1.000000e+00
## 195   4.391070e-17    0.2853018 -8.147278e+15  1.03660299  1.000000e+00
## 196   4.391070e-17    0.2853018  1.847979e+16 -3.59390957  1.000000e+00
## 197   4.391070e-17    0.2799671 -1.546003e+16  0.90328478  4.117844e-17
## 198             NA           NA            NA          NA  1.000000e+00
## 199             NA           NA            NA          NA  1.000000e+00
## 200   5.915411e-02    0.2866781  3.081732e+00  1.54178834  1.997542e-01
## 201             NA           NA            NA          NA  1.000000e+00
## 202             NA           NA            NA          NA  1.000000e+00
## 203             NA           NA            NA          NA  1.000000e+00
## 204             NA           NA            NA          NA  1.000000e+00
## 205             NA           NA            NA          NA  1.000000e+00
## 206             NA           NA            NA          NA  1.000000e+00
## 207             NA           NA            NA          NA  1.000000e+00
## 208   4.391070e-17    0.2853018 -3.615504e+15  0.70200639  1.000000e+00
## 209             NA           NA            NA          NA  1.000000e+00
## 210             NA           NA            NA          NA  1.000000e+00
## 211   4.391070e-17    0.2853018 -6.267639e+15  0.21495828  1.000000e+00
## 212             NA           NA            NA          NA  1.000000e+00
## 213             NA           NA            NA          NA  1.000000e+00
## 214   4.391070e-17    0.2853018 -9.543401e+15  2.79002793  1.000000e+00
## 215             NA           NA            NA          NA  1.000000e+00
## 216             NA           NA            NA          NA  1.000000e+00
## 217   9.452036e-02    0.2756237 -1.423527e+00 -0.15346306  2.044417e-01
## 218             NA           NA            NA          NA  1.000000e+00
## 219             NA           NA            NA          NA  1.000000e+00
## 220             NA           NA            NA          NA  1.000000e+00
## 221   4.391070e-17    0.2853018 -1.650711e+15  0.36050535  1.000000e+00
## 222   4.391070e-17    0.2853018 -2.598361e+16  3.38239743  1.000000e+00
## 223             NA           NA            NA          NA  1.000000e+00
## 224             NA           NA            NA          NA  1.000000e+00
## 225   4.391070e-17    0.2852965 -1.603997e+16  2.50242419  3.968959e-17
## 226             NA           NA            NA          NA  1.000000e+00
## 227             NA           NA            NA          NA  1.000000e+00
## 228             NA           NA            NA          NA  1.000000e+00
## 229   4.695327e-03    0.2799711 -1.456393e+01  0.44577475  4.364358e-02
## 230   4.391070e-17    0.2853018 -1.088457e+16  1.63613788  1.000000e+00
## 231   4.391070e-17    0.2853018  2.309671e+16 -1.89700053  1.000000e+00
## 232   4.391070e-17    0.2853018 -6.267639e+15  1.96206960  1.000000e+00
## 233             NA           NA            NA          NA  1.000000e+00
##     p_typeafter q_(Intercept) q_typeafter diff_(Intercept) diff_typeafter
## 1   1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 2   0.270626748  1.000000e+00   1.0000000            FALSE          FALSE
## 3   0.005951117  1.000000e+00   1.0000000            FALSE          FALSE
## 4   0.028450711  1.000000e+00   1.0000000            FALSE          FALSE
## 5   0.064205446  1.000000e+00   1.0000000            FALSE          FALSE
## 6   0.022750965  9.434995e-01   1.0000000            FALSE          FALSE
## 7   0.433592003  1.000000e+00   1.0000000            FALSE          FALSE
## 8   0.015054550  1.598612e-02   1.0000000             TRUE          FALSE
## 9   0.810487076  1.000000e+00   1.0000000            FALSE          FALSE
## 10  0.204190034  1.000000e+00   1.0000000            FALSE          FALSE
## 11  0.440214683  1.000000e+00   1.0000000            FALSE          FALSE
## 12  0.759004394  1.000000e+00   1.0000000            FALSE          FALSE
## 13  0.198668241  1.000000e+00   1.0000000            FALSE          FALSE
## 14  0.023046809  2.773385e-01   1.0000000            FALSE          FALSE
## 15  0.003395593  1.000000e+00   0.7877775            FALSE          FALSE
## 16  0.399387618  1.000000e+00   1.0000000            FALSE          FALSE
## 17  0.219989982  1.000000e+00   1.0000000            FALSE          FALSE
## 18  0.396621633  1.000000e+00   1.0000000            FALSE          FALSE
## 19  0.192462255  1.000000e+00   1.0000000            FALSE          FALSE
## 20  0.857074754  1.000000e+00   1.0000000            FALSE          FALSE
## 21  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 22  0.019327579  1.000000e+00   1.0000000            FALSE          FALSE
## 23  0.040154161  1.000000e+00   1.0000000            FALSE          FALSE
## 24  0.168230344  9.173797e-15   1.0000000             TRUE          FALSE
## 25  0.967128376  1.000000e+00   1.0000000            FALSE          FALSE
## 26  0.015103259  1.000000e+00   1.0000000            FALSE          FALSE
## 27  0.758365995  2.653997e-29   1.0000000             TRUE          FALSE
## 28  0.128247617  1.000000e+00   1.0000000            FALSE          FALSE
## 29  0.006931218  1.000000e+00   1.0000000            FALSE          FALSE
## 30  0.363921878  1.000000e+00   1.0000000            FALSE          FALSE
## 31  0.897453766  1.000000e+00   1.0000000            FALSE          FALSE
## 32  0.965711091  1.000000e+00   1.0000000            FALSE          FALSE
## 33  0.033307022  3.905146e-01   1.0000000            FALSE          FALSE
## 34  0.120361727  2.250302e-01   1.0000000            FALSE          FALSE
## 35  0.728845147  1.000000e+00   1.0000000            FALSE          FALSE
## 36  0.001102780  1.000000e+00   0.2569476            FALSE          FALSE
## 37  0.190024994  1.000000e+00   1.0000000            FALSE          FALSE
## 38  0.003418771  1.000000e+00   0.7897360            FALSE          FALSE
## 39  0.230069638  1.000000e+00   1.0000000            FALSE          FALSE
## 40  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 41  0.235224748  1.000000e+00   1.0000000            FALSE          FALSE
## 42  0.366906640  1.000000e+00   1.0000000            FALSE          FALSE
## 43  0.227791667  1.000000e+00   1.0000000            FALSE          FALSE
## 44  0.005352361  9.110809e-48   1.0000000             TRUE          FALSE
## 45  0.169478490  1.000000e+00   1.0000000            FALSE          FALSE
## 46  0.969888468  1.000000e+00   1.0000000            FALSE          FALSE
## 47  0.763170911  1.000000e+00   1.0000000            FALSE          FALSE
## 48  0.288797236  1.000000e+00   1.0000000            FALSE          FALSE
## 49  0.077802806  1.000000e+00   1.0000000            FALSE          FALSE
## 50  0.028300620  1.000000e+00   1.0000000            FALSE          FALSE
## 51  0.021833991  1.000000e+00   1.0000000            FALSE          FALSE
## 52  0.030653746  1.000000e+00   1.0000000            FALSE          FALSE
## 53  0.460640547  1.000000e+00   1.0000000            FALSE          FALSE
## 54  0.091736138  1.000000e+00   1.0000000            FALSE          FALSE
## 55  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 56  0.092936020  3.597889e-02   1.0000000             TRUE          FALSE
## 57  0.058003003  1.000000e+00   1.0000000            FALSE          FALSE
## 58  0.460748360  1.000000e+00   1.0000000            FALSE          FALSE
## 59  0.884698443  1.000000e+00   1.0000000            FALSE          FALSE
## 60  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 61  0.027254058  1.000000e+00   1.0000000            FALSE          FALSE
## 62  0.013396080  1.598151e-01   1.0000000            FALSE          FALSE
## 63  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 64  0.512965203  1.000000e+00   1.0000000            FALSE          FALSE
## 65  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 66  0.590320612  1.000000e+00   1.0000000            FALSE          FALSE
## 67  0.061775157  1.000000e+00   1.0000000            FALSE          FALSE
## 68  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 69  0.048390970  1.000000e+00   1.0000000            FALSE          FALSE
## 70  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 71  0.089381275  1.000000e+00   1.0000000            FALSE          FALSE
## 72  0.469924598  1.000000e+00   1.0000000            FALSE          FALSE
## 73  0.066524514  1.000000e+00   1.0000000            FALSE          FALSE
## 74  0.016775488  1.000000e+00   1.0000000            FALSE          FALSE
## 75  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 76  0.279651214  1.000000e+00   1.0000000            FALSE          FALSE
## 77  0.357684124  1.000000e+00   1.0000000            FALSE          FALSE
## 78  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 79  0.003746774  1.000000e+00   0.8617580            FALSE          FALSE
## 80  0.420774989  1.000000e+00   1.0000000            FALSE          FALSE
## 81  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 82  0.558036350  1.000000e+00   1.0000000            FALSE          FALSE
## 83  0.096699273  1.000000e+00   1.0000000            FALSE          FALSE
## 84  0.186777610  1.100159e-14   1.0000000             TRUE          FALSE
## 85  0.116184028  7.925304e-01   1.0000000            FALSE          FALSE
## 86  0.038370507  8.918364e-01   1.0000000            FALSE          FALSE
## 87  0.017603964  1.000000e+00   1.0000000            FALSE          FALSE
## 88  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 89  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 90  0.401454886  1.000000e+00   1.0000000            FALSE          FALSE
## 91  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 92  0.078225143  6.448487e-01   1.0000000            FALSE          FALSE
## 93  0.389707764  1.000000e+00   1.0000000            FALSE          FALSE
## 94  0.160432921  1.000000e+00   1.0000000            FALSE          FALSE
## 95  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 96  0.112500284  7.478736e-01   1.0000000            FALSE          FALSE
## 97  0.165056470  1.000000e+00   1.0000000            FALSE          FALSE
## 98  1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 99  0.534290066  1.000000e+00   1.0000000            FALSE          FALSE
## 100 0.290460421  1.000000e+00   1.0000000            FALSE          FALSE
## 101 0.368100070  1.000000e+00   1.0000000            FALSE          FALSE
## 102 0.044830158  9.817383e-03   1.0000000             TRUE          FALSE
## 103 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 104 0.324634480  1.000000e+00   1.0000000            FALSE          FALSE
## 105 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 106 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 107 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 108 0.789797672  1.000000e+00   1.0000000            FALSE          FALSE
## 109 0.344537636  1.000000e+00   1.0000000            FALSE          FALSE
## 110 0.400306433  1.000000e+00   1.0000000            FALSE          FALSE
## 111 0.041122347  1.155193e-02   1.0000000             TRUE          FALSE
## 112 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 113 0.322437374  1.000000e+00   1.0000000            FALSE          FALSE
## 114 0.160096065  1.000000e+00   1.0000000            FALSE          FALSE
## 115 0.070651862  1.000000e+00   1.0000000            FALSE          FALSE
## 116 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 117 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 118 0.131920500  1.000000e+00   1.0000000            FALSE          FALSE
## 119 0.319536610  1.000000e+00   1.0000000            FALSE          FALSE
## 120 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 121 0.387235181  1.000000e+00   1.0000000            FALSE          FALSE
## 122 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 123 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 124 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 125 0.253812202  1.032413e-14   1.0000000             TRUE          FALSE
## 126 0.318479024  1.000000e+00   1.0000000            FALSE          FALSE
## 127 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 128 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 129 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 130 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 131 0.154736580  5.635188e-15   1.0000000             TRUE          FALSE
## 132 0.976761150  1.000000e+00   1.0000000            FALSE          FALSE
## 133 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 134 0.937054478  1.000000e+00   1.0000000            FALSE          FALSE
## 135 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 136 0.172671683  1.000000e+00   1.0000000            FALSE          FALSE
## 137 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 138 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 139 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 140 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 141 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 142 0.074185042  3.032175e-15   1.0000000             TRUE          FALSE
## 143 0.237758339  1.000000e+00   1.0000000            FALSE          FALSE
## 144 0.301078343  1.000000e+00   1.0000000            FALSE          FALSE
## 145 0.052892388  1.000000e+00   1.0000000            FALSE          FALSE
## 146 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 147 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 148 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 149 0.245333855  1.000000e+00   1.0000000            FALSE          FALSE
## 150 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 151 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 152 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 153 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 154 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 155 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 156 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 157 0.036754868  7.597757e-01   1.0000000            FALSE          FALSE
## 158 0.135304588  1.000000e+00   1.0000000            FALSE          FALSE
## 159 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 160 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 161 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 162 0.137297982  5.635188e-15   1.0000000             TRUE          FALSE
## 163 0.240638933  1.000000e+00   1.0000000            FALSE          FALSE
## 164 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 165 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 166 0.106365982  4.326160e-01   1.0000000            FALSE          FALSE
## 167 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 168 0.128952269  1.000000e+00   1.0000000            FALSE          FALSE
## 169 0.056938009  5.046722e-01   1.0000000            FALSE          FALSE
## 170 0.292650693  1.000000e+00   1.0000000            FALSE          FALSE
## 171 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 172 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 173 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 174 0.985789447  8.302919e-01   1.0000000            FALSE          FALSE
## 175 0.691952539  1.000000e+00   1.0000000            FALSE          FALSE
## 176 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 177 0.933692271  1.000000e+00   1.0000000            FALSE          FALSE
## 178 0.715960843  1.000000e+00   1.0000000            FALSE          FALSE
## 179 0.725385232  1.000000e+00   1.0000000            FALSE          FALSE
## 180 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 181 0.253129364  1.354220e-14   1.0000000             TRUE          FALSE
## 182 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 183 0.944667980  1.000000e+00   1.0000000            FALSE          FALSE
## 184 0.049864455  1.000000e+00   1.0000000            FALSE          FALSE
## 185 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 186 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 187 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 188 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 189 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 190 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 191 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 192 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 193 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 194 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 195 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 196 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 197 0.532321966  9.306328e-15   1.0000000             TRUE          FALSE
## 198 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 199 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 200 0.366303773  1.000000e+00   1.0000000            FALSE          FALSE
## 201 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 202 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 203 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 204 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 205 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 206 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 207 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 208 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 209 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 210 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 211 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 212 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 213 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 214 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 215 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 216 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 217 0.883064416  1.000000e+00   1.0000000            FALSE          FALSE
## 218 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 219 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 220 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 221 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 222 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 223 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 224 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 225 0.242025194  9.049226e-15   1.0000000             TRUE          FALSE
## 226 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 227 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 228 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 229 0.733043267  1.000000e+00   1.0000000            FALSE          FALSE
## 230 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 231 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 232 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
## 233 1.000000000  1.000000e+00   1.0000000            FALSE          FALSE
##     passed_ss_(Intercept) passed_ss_typeafter
## 1                    TRUE                TRUE
## 2                    TRUE                TRUE
## 3                   FALSE               FALSE
## 4                   FALSE               FALSE
## 5                    TRUE               FALSE
## 6                    TRUE                TRUE
## 7                    TRUE                TRUE
## 8                   FALSE                TRUE
## 9                    TRUE                TRUE
## 10                   TRUE                TRUE
## 11                   TRUE                TRUE
## 12                   TRUE                TRUE
## 13                   TRUE                TRUE
## 14                   TRUE                TRUE
## 15                  FALSE               FALSE
## 16                   TRUE                TRUE
## 17                   TRUE                TRUE
## 18                   TRUE                TRUE
## 19                   TRUE                TRUE
## 20                   TRUE                TRUE
## 21                   TRUE               FALSE
## 22                  FALSE               FALSE
## 23                   TRUE               FALSE
## 24                  FALSE                TRUE
## 25                   TRUE                TRUE
## 26                   TRUE                TRUE
## 27                  FALSE                TRUE
## 28                   TRUE                TRUE
## 29                  FALSE               FALSE
## 30                   TRUE                TRUE
## 31                   TRUE                TRUE
## 32                   TRUE                TRUE
## 33                   TRUE                TRUE
## 34                   TRUE                TRUE
## 35                   TRUE                TRUE
## 36                   TRUE               FALSE
## 37                   TRUE                TRUE
## 38                   TRUE                TRUE
## 39                   TRUE                TRUE
## 40                   TRUE                TRUE
## 41                   TRUE                TRUE
## 42                   TRUE                TRUE
## 43                   TRUE                TRUE
## 44                  FALSE               FALSE
## 45                   TRUE                TRUE
## 46                   TRUE                TRUE
## 47                   TRUE                TRUE
## 48                   TRUE                TRUE
## 49                   TRUE                TRUE
## 50                   TRUE                TRUE
## 51                   TRUE               FALSE
## 52                   TRUE                TRUE
## 53                   TRUE                TRUE
## 54                   TRUE                TRUE
## 55                  FALSE               FALSE
## 56                  FALSE                TRUE
## 57                   TRUE                TRUE
## 58                   TRUE                TRUE
## 59                   TRUE                TRUE
## 60                   TRUE                TRUE
## 61                   TRUE                TRUE
## 62                   TRUE                TRUE
## 63                   TRUE                TRUE
## 64                   TRUE                TRUE
## 65                   TRUE                TRUE
## 66                   TRUE                TRUE
## 67                   TRUE                TRUE
## 68                   TRUE                TRUE
## 69                   TRUE                TRUE
## 70                   TRUE                TRUE
## 71                   TRUE                TRUE
## 72                   TRUE                TRUE
## 73                   TRUE                TRUE
## 74                   TRUE                TRUE
## 75                   TRUE                TRUE
## 76                   TRUE               FALSE
## 77                   TRUE                TRUE
## 78                   TRUE                TRUE
## 79                   TRUE               FALSE
## 80                   TRUE               FALSE
## 81                   TRUE                TRUE
## 82                   TRUE                TRUE
## 83                   TRUE                TRUE
## 84                  FALSE                TRUE
## 85                   TRUE                TRUE
## 86                   TRUE                TRUE
## 87                   TRUE                TRUE
## 88                   TRUE                TRUE
## 89                   TRUE                TRUE
## 90                   TRUE                TRUE
## 91                   TRUE                TRUE
## 92                   TRUE                TRUE
## 93                   TRUE                TRUE
## 94                   TRUE                TRUE
## 95                   TRUE                TRUE
## 96                   TRUE                TRUE
## 97                   TRUE                TRUE
## 98                   TRUE                TRUE
## 99                   TRUE                TRUE
## 100                  TRUE                TRUE
## 101                  TRUE                TRUE
## 102                 FALSE               FALSE
## 103                  TRUE                TRUE
## 104                  TRUE                TRUE
## 105                  TRUE                TRUE
## 106                  TRUE                TRUE
## 107                  TRUE                TRUE
## 108                  TRUE                TRUE
## 109                  TRUE                TRUE
## 110                  TRUE                TRUE
## 111                 FALSE                TRUE
## 112                  TRUE                TRUE
## 113                  TRUE                TRUE
## 114                  TRUE                TRUE
## 115                  TRUE                TRUE
## 116                  TRUE                TRUE
## 117                  TRUE                TRUE
## 118                  TRUE                TRUE
## 119                  TRUE                TRUE
## 120                  TRUE                TRUE
## 121                  TRUE                TRUE
## 122                  TRUE                TRUE
## 123                  TRUE                TRUE
## 124                  TRUE                TRUE
## 125                 FALSE                TRUE
## 126                  TRUE                TRUE
## 127                  TRUE                TRUE
## 128                  TRUE                TRUE
## 129                  TRUE                TRUE
## 130                  TRUE                TRUE
## 131                 FALSE                TRUE
## 132                  TRUE                TRUE
## 133                  TRUE                TRUE
## 134                  TRUE                TRUE
## 135                  TRUE                TRUE
## 136                  TRUE                TRUE
## 137                  TRUE                TRUE
## 138                  TRUE                TRUE
## 139                  TRUE                TRUE
## 140                  TRUE                TRUE
## 141                  TRUE                TRUE
## 142                 FALSE                TRUE
## 143                  TRUE                TRUE
## 144                  TRUE                TRUE
## 145                  TRUE                TRUE
## 146                  TRUE                TRUE
## 147                  TRUE                TRUE
## 148                  TRUE                TRUE
## 149                  TRUE                TRUE
## 150                  TRUE                TRUE
## 151                  TRUE                TRUE
## 152                  TRUE                TRUE
## 153                  TRUE                TRUE
## 154                  TRUE                TRUE
## 155                  TRUE                TRUE
## 156                  TRUE                TRUE
## 157                  TRUE                TRUE
## 158                  TRUE                TRUE
## 159                  TRUE                TRUE
## 160                  TRUE                TRUE
## 161                  TRUE                TRUE
## 162                 FALSE                TRUE
## 163                  TRUE                TRUE
## 164                  TRUE                TRUE
## 165                  TRUE                TRUE
## 166                  TRUE                TRUE
## 167                  TRUE                TRUE
## 168                  TRUE                TRUE
## 169                  TRUE                TRUE
## 170                  TRUE                TRUE
## 171                  TRUE                TRUE
## 172                  TRUE                TRUE
## 173                  TRUE                TRUE
## 174                  TRUE                TRUE
## 175                  TRUE                TRUE
## 176                  TRUE                TRUE
## 177                  TRUE                TRUE
## 178                  TRUE                TRUE
## 179                  TRUE                TRUE
## 180                  TRUE               FALSE
## 181                 FALSE                TRUE
## 182                  TRUE                TRUE
## 183                  TRUE                TRUE
## 184                  TRUE                TRUE
## 185                  TRUE                TRUE
## 186                  TRUE               FALSE
## 187                  TRUE                TRUE
## 188                  TRUE                TRUE
## 189                  TRUE                TRUE
## 190                  TRUE                TRUE
## 191                  TRUE                TRUE
## 192                  TRUE                TRUE
## 193                  TRUE               FALSE
## 194                  TRUE                TRUE
## 195                  TRUE                TRUE
## 196                  TRUE                TRUE
## 197                 FALSE                TRUE
## 198                  TRUE                TRUE
## 199                  TRUE                TRUE
## 200                  TRUE                TRUE
## 201                  TRUE                TRUE
## 202                  TRUE                TRUE
## 203                  TRUE                TRUE
## 204                  TRUE                TRUE
## 205                  TRUE               FALSE
## 206                  TRUE                TRUE
## 207                  TRUE                TRUE
## 208                  TRUE                TRUE
## 209                  TRUE                TRUE
## 210                  TRUE                TRUE
## 211                  TRUE                TRUE
## 212                  TRUE                TRUE
## 213                  TRUE               FALSE
## 214                  TRUE                TRUE
## 215                  TRUE                TRUE
## 216                  TRUE                TRUE
## 217                  TRUE                TRUE
## 218                  TRUE                TRUE
## 219                  TRUE                TRUE
## 220                  TRUE                TRUE
## 221                  TRUE                TRUE
## 222                  TRUE                TRUE
## 223                  TRUE                TRUE
## 224                  TRUE               FALSE
## 225                 FALSE                TRUE
## 226                  TRUE                TRUE
## 227                  TRUE                TRUE
## 228                  TRUE                TRUE
## 229                  TRUE                TRUE
## 230                  TRUE                TRUE
## 231                  TRUE                TRUE
## 232                  TRUE                TRUE
## 233                  TRUE                TRUE
```

``` r
write.table(res_prim_dna,"ANCOM-BC2_DNA.txt",sep="\t",col.names=NA)
# nothing significantly different after treatment relative to before (all taxa are FALSE for column "diff_typeafter")
```


Differential abundance testing using ANCOM-BC: RNA fraction only
https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html



``` r
# Start by creating phyloseq object using full dataset (low abundance ASVs are not removed).

otu <- read.table("silva_nochloronomito_otu_table_ps1.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("silva_nochloronomito_taxa_table_ps1.txt",sep="\t",header=TRUE,row.names=1)
samples<-read.table("metadata_ps1.txt",sep="\t",header=T,row.names=1)
OTU = otu_table(otu, taxa_are_rows=FALSE)
taxon<-as.matrix(taxon)
TAX = tax_table(taxon)
sampledata = sample_data(samples)
ps1 <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(taxon))
ps1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 6137 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 6137 taxa by 6 taxonomic ranks ]
```

``` r
# prune empty rows
ps_qf <- prune_taxa(taxa_sums(ps1) > 1, ps1) 
ps_qf
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1855 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 1855 taxa by 6 taxonomic ranks ]
```

``` r
# Subset and test RNA fraction
ps_qf
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1855 taxa and 35 samples ]
## sample_data() Sample Data:       [ 35 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 1855 taxa by 6 taxonomic ranks ]
```

``` r
ps_rna <- subset_samples(ps_qf, fraction == "RNA")
ps_rna 
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1855 taxa and 17 samples ]
## sample_data() Sample Data:       [ 17 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 1855 taxa by 6 taxonomic ranks ]
```

``` r
# prune empty rows
ps_rna_qf <- prune_taxa(taxa_sums(ps_rna) > 1, ps_rna) 
ps_rna_qf
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 889 taxa and 17 samples ]
## sample_data() Sample Data:       [ 17 samples by 4 sample variables ]
## tax_table()   Taxonomy Table:    [ 889 taxa by 6 taxonomic ranks ]
```

``` r
sample_data(ps_rna_qf)$type<-factor(sample_data(ps_rna_qf)$type,levels=c("before","after"))
# You can verify the change by checking:
levels(sample_data(ps_rna_qf)$type)
```

```
## [1] "before" "after"
```

``` r
set.seed(123)
output = ancombc2(data = ps_rna_qf, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "type", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "type", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = FALSE, 
                  dunnet = FALSE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, mdfdr_control = NULL, 
                  trend_control = NULL)
```

```
## Warning: The group variable has < 3 categories 
## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
```

```
## Obtaining initial estimates ...
```

```
## Estimating sample-specific biases ...
```

```
## Sensitivity analysis for pseudo-count addition to 0s: ...
```

```
## ANCOM-BC2 primary results ...
```

```
## Warning in pt(abs(W), df = dof, lower.tail = FALSE): NaNs produced
```

``` r
res_prim_rna = output$res
print(res_prim_rna)
```

```
##                                                  taxon lfc_(Intercept)
## 1                                          Marinifilum    -0.549124078
## 2                                               Vibrio    -0.046543481
## 3                                            Clade III     0.076042305
## 4                                             Clade Ia     0.325461660
## 5                                       Proteobacteria     0.579587946
## 6                              Candidatus Actinomarina    -0.763628778
## 7                                         Sphingopyxis    -0.189960738
## 8                                     Alteromonadaceae    -0.745966996
## 9                                 Synechococcus CC9902     0.536967616
## 10                                Symphothece PCC-7002     0.348123275
## 11                                          Litoricola     0.306726654
## 12                             AEGEAN-169 marine group    -0.119257155
## 13                                             Delftia     0.613440432
## 14                                   Blastocatellaceae    -0.100270234
## 15                                         Pseudomonas     0.636224063
## 16                                            Clade II    -0.076093207
## 17                                              Kordia    -0.711862871
## 18                                              HIMB11     0.329135474
## 19                                        SAR116 clade     0.725040206
## 20                                         SAR86 clade    -0.089168447
## 21                                   Halodesulfovibrio    -0.182205503
## 22                                    NS5 marine group     0.342275005
## 23                                   Campylobacterales     0.212236654
## 24                                    OM60(NOR5) clade     0.392967029
## 25                                             S25-593     0.623448022
## 26                                      Cryomorphaceae     0.073191313
## 27                                    NS4 marine group     0.013240209
## 28                                NS11-12 marine group    -0.112297000
## 29                                     Aestuariibacter     0.203392192
## 30                  Urania-1B-19 marine sediment group     0.006098683
## 31                                   Flavobacteriaceae    -0.152627791
## 32                         Candidatus Puniceispirillum    -0.175862887
## 33                                    NS9 marine group     0.230245021
## 34                                            Bacteria     0.144784856
## 35                                   Cyclobacteriaceae     0.133505793
## 36                                     Pseudomonadales     0.110476853
## 37                                             MBMPE27    -0.053941445
## 38                                         Fokiniaceae    -0.216471090
## 39                                            Balneola     0.073211077
## 40                                         KI89A clade              NA
## 41                              Synechococcus PCC-7336              NA
## 42                                   Pseudoalteromonas    -0.183164427
## 43                                       Thalassotalea     0.831922321
## 44                       Marinimicrobia (SAR406 clade)     0.115969818
## 45                                 Alphaproteobacteria     0.116062213
## 46                                      Rhodopirellula    -0.488770666
## 47  Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium     0.233002947
## 48                                     Francisellaceae     0.626040088
## 49                                           PS1 clade     0.124355346
## 50                                      Bradymonadales    -0.017038132
## 51                                    Terasakiellaceae    -0.699203967
## 52                                      Staphylococcus     0.202364217
## 53                                         Azotobacter     0.284630305
## 54                                     Pedosphaeraceae     0.298792657
## 55                                         Alteromonas     0.165992394
## 56                             Candidatus Amoebophilus    -0.417017529
## 57                                    Coraliomargarita     0.027681313
## 58                                       Cohaesibacter              NA
## 59                           Candidatus Cryptoprodotis              NA
## 60                                       Rickettsiales     0.611174349
## 61                                      Comamonadaceae     0.722503330
## 62                                             Ekhidna    -1.428422682
## 63                                          Aureitalea     1.063229600
## 64                                       Acinetobacter     0.249221676
## 65                                 Xenococcus PCC-7305              NA
## 66                                      Cyanobacteriia              NA
## 67                                    Desulfocapsaceae     0.291364360
## 68                       Marine Methylotrophic Group 3     0.382712544
## 69                                     Marine Group II    -0.868086717
## 70                                              Dyella    -1.379211491
## 71                                      Saprospiraceae     0.424029593
## 72                                      Unknown Family    -0.242294885
## 73                                            R76-B128    -0.330230482
## 74                                     SCGC AAA164-E04    -0.092190498
## 75                               Acrophormium PCC-7375              NA
## 76                                 Gammaproteobacteria     0.271361570
## 77                                      Propionigenium              NA
## 78                                     Chitinophagales    -0.258924158
## 79                                          Lawsonella              NA
## 80                                    Rhodobacteraceae              NA
## 81                                     Hyphomonadaceae     0.282616177
## 82                                     Corynebacterium     1.582107681
## 83                                         Alcanivorax    -0.583213088
## 84                                      Endozoicomonas    -0.240677547
## 85                                        Chromatiales    -0.518298266
## 86                  [Caedibacter] taeniospiralis group     0.135592619
## 87                                     JGI 0000069-P22              NA
## 88                                               NB1-j              NA
## 89                                          OM27 clade              NA
## 90                                       Bacteroidales              NA
## 91                                         Stappiaceae              NA
## 92                                  Candidatus Megaira     0.290326120
## 93                                    Phormidesmiaceae     2.279478962
## 94                                    Oscillatoriaceae              NA
## 95                                        Pla4 lineage     0.055778090
## 96                                      Oligoflexaceae    -0.300645054
## 97                                          Ferrimonas              NA
## 98                                       Enhydrobacter    -0.106558452
## 99                                Trichodesmium IMS101    -0.346169801
## 100                                     Microcystaceae              NA
## 101                                         Acidovorax              NA
## 102                                       Crocinitomix    -0.321774719
## 103                                         Rubritalea     0.537838919
## 104                                      Halarcobacter              NA
## 105                                 Cyanobium PCC-6307    -0.793853442
## 106                                      Aquabacterium              NA
## 107                                     Photobacterium              NA
## 108                                           Reinekea              NA
## 109                                         Halieaceae              NA
## 110                                   Unknown Family_1              NA
## 111                                               PB19              NA
## 112                                          Lewinella              NA
## 113                     Methylobacterium-Methylorubrum              NA
## 114                                      Marinoscillum              NA
## 115                                     Psychrosphaera              NA
## 116                                       FS142-36B-02    -0.525082036
## 117                       SAR324 clade(Marine group B)    -0.563722882
## 118                                             BD7-11              NA
## 119                                         Aureispira              NA
## 120                                       Lentisphaera    -0.301702999
## 121                                      Streptococcus    -0.679738752
## 122                                        OM182 clade    -0.755573291
## 123                                    Cloacibacterium     0.692042260
## 124                                         Shewanella     0.027203365
## 125                                        Turneriella              NA
## 126                                     Niveispirillum     0.442456797
## 127                                           Coxiella    -0.169046207
## 128                                        Sphingobium              NA
## 129                                           Azospira     0.701900883
## 130                                           PRD18C08              NA
## 131                                        Marinomonas              NA
## 132                                     Flavobacterium              NA
## 133                                     Desulfosarcina    -0.094406640
## 134                                   Latescibacterota    -1.738924738
## 135                                   Cellvibrionaceae    -0.471238406
## 136                                       Acholeplasma              NA
## 137                                       Cytophagales              NA
## 138                                       Chlamydiales     0.213327218
## 139                                          vadinHA49    -0.850772085
## 140                                        BD1-7 clade     0.735654927
## 141                                      Lactobacillus    -0.443936258
## 142                                      Spirochaetota              NA
## 143                                         Aliicoccus              NA
## 144                             Sva0081 sediment group              NA
## 145                                       Sphingomonas              NA
## 146                                   Sphingobacterium    -0.220358240
## 147                                       Pla3 lineage    -0.066361907
## 148                                     Rubritaleaceae    -0.902493424
## 149                                   Flammeovirgaceae     0.020915726
## 150                                   Chryseobacterium              NA
## 151                                        Caulobacter              NA
## 152                                            Woeseia     0.126532781
## 153                                                A4b              NA
## 154                                      Colwelliaceae    -1.483511926
## 155                                        Bacteroidia    -0.592366486
## 156                                           Vicingus              NA
## 157                                            Dietzia    -0.359391579
## 158                                          Algivirga    -0.254689075
## 159                                   Stenotrophomonas    -1.882765774
## 160                                   Flavobacteriales    -0.594921909
## 161                                              HOC36              NA
## 162                                         Caldithrix              NA
## 163                                       Terasakiella              NA
## 164                                           WCHB1-41    -0.346169801
## 165                                   Cyanobacteriales              NA
## 166                                     Agaribacterium     0.280211683
## 167                                    Amoebophilaceae    -0.387186155
## 168                                             KD4-96              NA
## 169                                 Verrucomicrobiales              NA
## 170                                 Desulfobacteraceae    -0.727239827
## 171                                           Zoogloea              NA
## 172                                         Lentimonas              NA
## 173                                       Peredibacter    -0.244699379
## 174                                    Planctomycetota    -1.189618594
## 175                                       Anaerococcus              NA
## 176                                        Kordiimonas     0.797425886
## 177                                              OM190    -0.205288320
## 178                                             DEV007              NA
## 179                                 Bacteriovoracaceae    -0.073379386
## 180                                    Woesearchaeales    -0.165500127
## 181                                   Rhodospirillales              NA
## 182                                          Thiovulum     0.464937734
## 183                                Bacteroidetes BD2-2     0.253315193
## 184                                                TG3              NA
## 185                                      Brevundimonas              NA
## 186                                 Spongiibacteraceae              NA
## 187                                   Ardenticatenales     0.201877745
## 188                                   Margulisbacteria              NA
## 189                            Planktothricoides SR001    -1.333459630
## 190                                       Cyanobiaceae    -1.189618594
## 191                                    Spirochaetaceae              NA
##     lfc_typeafter se_(Intercept) se_typeafter W_(Intercept) W_typeafter
## 1      0.48782484   7.569282e-02    0.2713301 -7.254639e+00  1.79790179
## 2      0.09078514   2.011506e-01    0.3683711 -2.313862e-01  0.24645023
## 3     -0.13118871   1.888716e-01    0.4543234  4.026138e-01 -0.28875625
## 4     -0.63002742   2.940896e-01    0.4897129  1.106675e+00 -1.28652394
## 5     -1.33579105   5.373798e-01    0.7696023  1.078544e+00 -1.73569009
## 6      1.67292147   3.653993e-01    0.4809101 -2.089847e+00  3.47865736
## 7      0.48369446   9.195517e-01    1.0500937 -2.065797e-01  0.46062028
## 8      0.46299998   9.283754e-17    0.2754077 -8.035187e+15  1.68114382
## 9     -1.07623705   2.979510e-01    0.6409397  1.802201e+00 -1.67915499
## 10    -0.64688024   9.283754e-17    0.2696754  3.749811e+15 -2.39873622
## 11    -0.65503371   3.567333e-01    0.5114996  8.598207e-01 -1.28061430
## 12     0.33050336   4.561711e-01    0.5675909 -2.614308e-01  0.58229151
## 13    -0.68209323   2.969246e-01    0.5180698  2.065980e+00 -1.31660476
## 14     0.24499391   4.707630e-01    0.7313595 -2.129952e-01  0.33498424
## 15    -0.72006595   2.189513e-01    0.4717684  2.905779e+00 -1.52631225
## 16     0.05620834   2.798407e-01    0.3917332 -2.719162e-01  0.14348627
## 17    -0.20277423   9.283754e-17    0.2754077 -7.667834e+15 -0.73626923
## 18    -0.64101400   2.359392e-01    0.3312958  1.395001e+00 -1.93486898
## 19    -1.49880759   4.128430e-01    0.5240141  1.756213e+00 -2.86024311
## 20     0.26840279   2.315325e-01    0.3721541 -3.851227e-01  0.72121404
## 21     0.26790455   1.224515e-01    0.4202012 -1.487980e+00  0.63756258
## 22    -0.67473331   2.986205e-01    0.5756799  1.146187e+00 -1.17206334
## 23    -0.67227764   4.753298e-01    0.6283494  4.465040e-01 -1.06991045
## 24    -1.10699618   4.500637e-01    0.5519206  8.731365e-01 -2.00571652
## 25    -1.42137442   2.022589e-01    0.4612189  3.082426e+00 -3.08177820
## 26    -0.60366190   4.685301e-01    0.5487053  1.562147e-01 -1.10015687
## 27    -0.33332816   2.332470e-01    0.3528604  5.676476e-02 -0.94464605
## 28     0.14179390   2.819483e-01    0.4121654 -3.982893e-01  0.34402189
## 29    -1.61612564   2.968284e-02    0.2714830  6.852181e+00 -5.95295371
## 30    -0.37836216   1.181045e-01    0.3036121  5.163803e-02 -1.24620271
## 31     0.32557195   2.485832e-01    0.8414503 -6.139908e-01  0.38691763
## 32    -0.67651048   1.100261e-01    0.3056338 -1.598374e+00 -2.21346717
## 33    -0.54143815   3.602958e-01    0.4837080  6.390444e-01 -1.11934917
## 34    -0.22093987   5.475116e-01    0.6191938  2.644416e-01 -0.35681859
## 35    -0.51298364   1.631396e-01    0.3640591  8.183531e-01 -1.40906707
## 36    -0.92325016   2.863585e-01    0.3952504  3.857991e-01 -2.33586120
## 37     0.05098669   3.388333e-01    0.4372346 -1.591976e-01  0.11661176
## 38    -0.01515143   4.046364e-01    0.4957672 -5.349768e-01 -0.03056158
## 39    -0.12130807   2.362298e-01    0.3688944  3.099147e-01 -0.32884223
## 40             NA             NA           NA            NA          NA
## 41             NA             NA           NA            NA          NA
## 42     0.13101662   2.394257e-01    0.3706545 -7.650157e-01  0.35347365
## 43    -0.92846448   9.283754e-17    0.2754077  8.961055e+15 -3.37123627
## 44    -0.21530398   2.249454e-01    0.4007591  5.155465e-01 -0.53724039
## 45    -0.66700270   2.267480e-01    0.4543469  5.118554e-01 -1.46804737
## 46     1.40836602   2.262744e-02    0.2707109 -2.160079e+01  5.20247207
## 47    -0.28513435   4.414360e-01    0.6367849  5.278295e-01 -0.44777183
## 48    -1.74908948   1.790461e-01    0.3244278  3.496531e+00 -5.39130624
## 49    -0.55756252   2.926646e-01    0.3894582  4.249074e-01 -1.43163641
## 50    -0.45922477   9.106192e-02    0.2790406 -1.871049e-01 -1.64572733
## 51     0.02250310   9.283754e-17    0.2754077 -7.531479e+15  0.08170833
## 52    -0.36160242   2.097660e-01    0.3474736  9.647143e-01 -1.04066158
## 53    -0.41494150   2.123501e-01    0.3521533  1.340382e+00 -1.17829802
## 54    -1.01298213   1.524974e-01    0.3535765  1.959330e+00 -2.86495927
## 55    -0.20585452   1.705493e-01    0.3400813  9.732812e-01 -0.60530967
## 56     0.59518714   4.020491e-01    0.4881670 -1.037230e+00  1.21922857
## 57    -0.52327655   3.230692e-02    0.3379037  8.568232e-01 -1.54859668
## 58             NA             NA           NA            NA          NA
## 59             NA             NA           NA            NA          NA
## 60    -2.92735101   4.419248e-01    0.5281422  1.382983e+00 -5.54273250
## 61    -2.08198367   4.270063e-01    0.6133283  1.692020e+00 -3.39456623
## 62     4.26109999   4.322592e-02    0.2677111 -3.304551e+01 15.91678300
## 63    -1.63058586   9.283754e-17    0.2754077  1.145258e+16 -5.92062521
## 64    -0.46786744   5.698627e-01    0.6416479  4.373364e-01 -0.72916542
## 65             NA             NA           NA            NA          NA
## 66             NA             NA           NA            NA          NA
## 67    -0.65819685   7.004595e-02    0.2759112  4.159617e+00 -2.38553860
## 68    -0.25663604   3.038104e-01    0.4155395  1.259708e+00 -0.61759722
## 69     0.81420297   1.948415e-01    0.3426911 -4.455348e+00  2.37590950
## 70     3.22230178   9.283754e-17    0.2754077 -1.485618e+16 11.70011442
## 71    -1.03639630   1.474374e-01    0.3064416  2.875997e+00 -3.38203577
## 72     0.42148535   9.283754e-17    0.2754077 -2.609880e+15  1.53040502
## 73     0.41860659   1.147357e-02    0.3053652 -2.878185e+01  1.37083937
## 74    -0.84281986   6.963099e-03    0.2675744 -1.323987e+01 -3.14985268
## 75             NA             NA           NA            NA          NA
## 76    -0.89812550   3.910563e-01    0.4726420  6.939194e-01 -1.90022350
## 77             NA             NA           NA            NA          NA
## 78     0.27196664   3.378671e-01    0.4618308 -7.663490e-01  0.58888809
## 79             NA             NA           NA            NA          NA
## 80             NA             NA           NA            NA          NA
## 81    -0.33941739   1.103013e-16    0.3014877  2.562221e+15 -1.12580849
## 82    -2.22444308   9.283754e-17    0.2728058  1.704168e+16 -8.15394247
## 83     0.28680320   9.354009e-02    0.2818638 -6.234900e+00  1.01752416
## 84     1.01077113   1.661713e-01    0.3216012 -1.448370e+00  3.14293312
## 85    -0.24953885   5.229664e-02    0.2753198 -9.910738e+00 -0.90635985
## 86    -0.99326171   2.546272e-01    0.3712271  5.325143e-01 -2.67561754
## 87             NA             NA           NA            NA          NA
## 88             NA             NA           NA            NA          NA
## 89             NA             NA           NA            NA          NA
## 90             NA             NA           NA            NA          NA
## 91             NA             NA           NA            NA          NA
## 92    -0.56239555   1.031222e-01    0.2909735  2.815360e+00 -1.93280656
## 93    -3.37042836   1.377987e-16    0.3772604  1.654210e+16 -8.93395768
## 94             NA             NA           NA            NA          NA
## 95     0.73128797   9.283754e-17    0.2754077  6.008140e+14  2.65529226
## 96     0.29118912   1.272074e-01    0.2948233 -2.363424e+00  0.98767326
## 97             NA             NA           NA            NA          NA
## 98     0.42881564   1.757693e-01    0.3224592 -6.062403e-01  1.32982916
## 99     0.32811968   9.283754e-17    0.2754077 -3.728770e+15  1.19139612
## 100            NA             NA           NA            NA          NA
## 101            NA             NA           NA            NA          NA
## 102    0.02873438   9.283754e-17    0.2754077 -3.465998e+15  0.10433396
## 103   -0.80569749   1.107233e-02    0.2642646  4.857506e+01 -3.04882820
## 104            NA             NA           NA            NA          NA
## 105    1.72282303   4.943673e-02    0.2689310 -1.605797e+01  6.40619084
## 106            NA             NA           NA            NA          NA
## 107            NA             NA           NA            NA          NA
## 108            NA             NA           NA            NA          NA
## 109            NA             NA           NA            NA          NA
## 110            NA             NA           NA            NA          NA
## 111            NA             NA           NA            NA          NA
## 112            NA             NA           NA            NA          NA
## 113            NA             NA           NA            NA          NA
## 114            NA             NA           NA            NA          NA
## 115            NA             NA           NA            NA          NA
## 116    0.31969242   1.551382e-01    0.3168072 -3.384608e+00  1.00910731
## 117    1.41769324   5.199325e-02    0.2694783 -1.084223e+01  5.26088159
## 118            NA             NA           NA            NA          NA
## 119            NA             NA           NA            NA          NA
## 120    0.57449331   1.020413e-16    0.2947186 -2.956674e+15  1.94929450
## 121    1.82335630   9.283754e-17    0.2754077 -7.321809e+15  6.62057089
## 122    0.46228622   1.987168e-01    0.3411601 -3.802261e+00  1.35504195
## 123   -1.46418090   1.517489e-01    0.3087565  4.560442e+00 -4.74218622
## 124   -0.11751115   9.283754e-17    0.2754077  2.930212e+14 -0.42668068
## 125            NA             NA           NA            NA          NA
## 126   -0.82161812   9.283754e-17    0.2818758  4.765925e+15 -2.91482362
## 127   -0.17806875   6.522764e-02    0.2887557 -2.591635e+00 -0.61667615
## 128            NA             NA           NA            NA          NA
## 129   -0.78089234   9.283754e-17    0.2754077  7.560529e+15 -2.83540474
## 130            NA             NA           NA            NA          NA
## 131            NA             NA           NA            NA          NA
## 132            NA             NA           NA            NA          NA
## 133   -0.17540664   9.283754e-17    0.2754077 -1.016902e+15 -0.63689806
## 134    2.10194465   9.283754e-17    0.2754077 -1.873084e+16  7.63211969
## 135   -0.05288313   4.968642e-03    0.2656354 -9.484249e+01 -0.19908165
## 136            NA             NA           NA            NA          NA
## 137            NA             NA           NA            NA          NA
## 138    0.01615652   9.283754e-17    0.2728058  2.297855e+15  0.05922351
## 139    0.68420078   9.576277e-02    0.2881512 -8.884163e+00  2.37445031
## 140   -1.01833675   1.020413e-16    0.2942359  7.209381e+15 -3.46095355
## 141    0.62591330   6.158063e-02    0.2774782 -7.209024e+00  2.25572038
## 142            NA             NA           NA            NA          NA
## 143            NA             NA           NA            NA          NA
## 144            NA             NA           NA            NA          NA
## 145            NA             NA           NA            NA          NA
## 146    1.22857832   1.259769e-01    0.3008112 -1.749195e+00  4.08421780
## 147    0.11136717   9.283754e-17    0.2754077 -7.148176e+14  0.40437200
## 148    1.71994838   9.283754e-17    0.2754077 -9.721212e+15  6.24509876
## 149    0.49469930   3.034460e-02    0.2715659  6.892734e-01  1.82165460
## 150            NA             NA           NA            NA          NA
## 151            NA             NA           NA            NA          NA
## 152    0.11102853   1.405994e-02    0.2944212  8.999525e+00  0.37710775
## 153            NA             NA           NA            NA          NA
## 154    1.59111902   9.283754e-17    0.2754077 -1.597966e+16  5.77732189
## 155    1.88301864   3.079168e-01    0.4121981 -1.923787e+00  4.56823729
## 156            NA             NA           NA            NA          NA
## 157    0.90627388   9.283754e-17    0.2754077 -3.871188e+15  3.29066263
## 158   -0.05035795   1.249219e-16    0.3759416 -2.038786e+15 -0.13395152
## 159    2.38962672   9.283754e-17    0.2754077 -2.028022e+16  8.67668764
## 160    1.52179358   1.174094e-01    0.2969399 -5.067071e+00  5.12492065
## 161            NA             NA           NA            NA          NA
## 162            NA             NA           NA            NA          NA
## 163            NA             NA           NA            NA          NA
## 164    0.67509050   9.283754e-17    0.2754077 -3.728770e+15  2.45124033
## 165            NA             NA           NA            NA          NA
## 166   -0.56193569   9.283754e-17    0.2754077  3.018301e+15 -2.04037746
## 167    0.14852690   1.230819e-01    0.2890757 -3.145760e+00  0.51379939
## 168            NA             NA           NA            NA          NA
## 169            NA             NA           NA            NA          NA
## 170    1.63814918   9.283754e-17    0.2754077 -7.833467e+15  5.94808747
## 171            NA             NA           NA            NA          NA
## 172            NA             NA           NA            NA          NA
## 173    0.12517884   9.283754e-17    0.2754077 -2.635780e+15  0.45452190
## 174    1.00333236   9.283754e-17    0.2754077 -1.281398e+16  3.64308007
## 175            NA             NA           NA            NA          NA
## 176   -1.06039024   9.283754e-17    0.2754077  8.589477e+15 -3.85025610
## 177    0.15680036   2.217664e-01    0.3444502 -9.256964e-01  0.45521919
## 178            NA             NA           NA            NA          NA
## 179   -0.53173811   9.283754e-17    0.2754077 -7.904064e+14 -1.93073062
## 180    0.33014409   1.037259e-01    0.3334624 -1.595553e+00  0.99004895
## 181            NA             NA           NA            NA          NA
## 182   -0.39952147   9.283754e-17    0.2754077  5.008079e+15 -1.45065462
## 183   -0.01630120   6.021163e-02    0.2891817  4.207081e+00 -0.05637012
## 184            NA             NA           NA            NA          NA
## 185            NA             NA           NA            NA          NA
## 186            NA             NA           NA            NA          NA
## 187   -0.22528426   9.283754e-17    0.2825092  2.174527e+15 -0.79744042
## 188            NA             NA           NA            NA          NA
## 189    1.04041929   9.283754e-17    0.2754077 -1.436337e+16  3.77774199
## 190    1.34619563   9.283754e-17    0.2754077 -1.281398e+16  4.88800987
## 191            NA             NA           NA            NA          NA
##     p_(Intercept) p_typeafter q_(Intercept) q_typeafter diff_(Intercept)
## 1    5.403641e-03 0.170036870  9.510409e-01   1.0000000            FALSE
## 2    8.216806e-01 0.810318927  1.000000e+00   1.0000000            FALSE
## 3    6.957044e-01 0.778665087  1.000000e+00   1.0000000            FALSE
## 4    2.943452e-01 0.227245625  1.000000e+00   1.0000000            FALSE
## 5    3.061189e-01 0.113270191  1.000000e+00   1.0000000            FALSE
## 6    6.619974e-02 0.006953249  1.000000e+00   1.0000000            FALSE
## 7    8.401121e-01 0.654043748  1.000000e+00   1.0000000            FALSE
## 8    1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 9    1.016870e-01 0.124045774  1.000000e+00   1.0000000            FALSE
## 10   1.697738e-16 0.251450902  3.140815e-14   1.0000000             TRUE
## 11   4.100272e-01 0.229233423  1.000000e+00   1.0000000            FALSE
## 12   7.985906e-01 0.572112075  1.000000e+00   1.0000000            FALSE
## 13   7.268163e-02 0.224436701  1.000000e+00   1.0000000            FALSE
## 14   8.360771e-01 0.745312386  1.000000e+00   1.0000000            FALSE
## 15   1.971625e-02 0.165446661  1.000000e+00   1.0000000            FALSE
## 16   7.925722e-01 0.889454574  1.000000e+00   1.0000000            FALSE
## 17   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 18   1.905430e-01 0.079128974  1.000000e+00   1.0000000            FALSE
## 19   1.068177e-01 0.015509896  1.000000e+00   1.0000000            FALSE
## 20   7.082111e-01 0.487296838  1.000000e+00   1.0000000            FALSE
## 21   2.109753e-01 0.558425379  1.000000e+00   1.0000000            FALSE
## 22   2.784045e-01 0.268337965  1.000000e+00   1.0000000            FALSE
## 23   6.670723e-01 0.315875037  1.000000e+00   1.0000000            FALSE
## 24   4.052820e-01 0.075850368  1.000000e+00   1.0000000            FALSE
## 25   1.506142e-02 0.015076221  1.000000e+00   1.0000000            FALSE
## 26   8.819751e-01 0.321388662  1.000000e+00   1.0000000            FALSE
## 27   9.563189e-01 0.376294991  1.000000e+00   1.0000000            FALSE
## 28   7.042030e-01 0.742565884  1.000000e+00   1.0000000            FALSE
## 29   9.225633e-02 0.105952612  1.000000e+00   1.0000000            FALSE
## 30   9.604934e-01 0.259139637  1.000000e+00   1.0000000            FALSE
## 31   5.544219e-01 0.707806979  1.000000e+00   1.0000000            FALSE
## 32   2.082531e-01 0.113733288  1.000000e+00   1.0000000            FALSE
## 33   5.387159e-01 0.291974007  1.000000e+00   1.0000000            FALSE
## 34   7.973958e-01 0.729450483  1.000000e+00   1.0000000            FALSE
## 35   4.444335e-01 0.208475141  1.000000e+00   1.0000000            FALSE
## 36   7.253775e-01 0.101603765  1.000000e+00   1.0000000            FALSE
## 37   8.836271e-01 0.914536078  1.000000e+00   1.0000000            FALSE
## 38   6.461837e-01 0.978394746  1.000000e+00   1.0000000            FALSE
## 39   7.629879e-01 0.749057976  1.000000e+00   1.0000000            FALSE
## 40   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 41   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 42   4.869046e-01 0.741576435  1.000000e+00   1.0000000            FALSE
## 43   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 44   6.246093e-01 0.610430720  1.000000e+00   1.0000000            FALSE
## 45   6.305501e-01 0.202020471  1.000000e+00   1.0000000            FALSE
## 46   2.945103e-02 0.120894200  1.000000e+00   1.0000000            FALSE
## 47   6.103821e-01 0.664901575  1.000000e+00   1.0000000            FALSE
## 48   7.295650e-02 0.032724799  1.000000e+00   1.0000000            FALSE
## 49   6.836601e-01 0.195344233  1.000000e+00   1.0000000            FALSE
## 50   8.635171e-01 0.198371732  1.000000e+00   1.0000000            FALSE
## 51   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 52   3.629402e-01 0.328458827  1.000000e+00   1.0000000            FALSE
## 53   2.377994e-01 0.291696132  1.000000e+00   1.0000000            FALSE
## 54   9.090938e-02 0.024165560  1.000000e+00   1.0000000            FALSE
## 55   4.021986e-01 0.587695535  1.000000e+00   1.0000000            FALSE
## 56   3.582130e-01 0.289736049  1.000000e+00   1.0000000            FALSE
## 57   4.306838e-01 0.182160632  1.000000e+00   1.0000000            FALSE
## 58   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 59   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 60   3.008317e-01 0.031042403  1.000000e+00   1.0000000            FALSE
## 61   1.658994e-01 0.027414281  1.000000e+00   1.0000000            FALSE
## 62   9.144900e-04 0.003923985  1.627792e-01   0.7494811            FALSE
## 63   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 64   6.703296e-01 0.481138110  1.000000e+00   1.0000000            FALSE
## 65   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 66   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 67   5.322335e-02 0.139797209  1.000000e+00   1.0000000            FALSE
## 68   2.968603e-01 0.580553892  1.000000e+00   1.0000000            FALSE
## 69   2.104648e-02 0.097979138  1.000000e+00   1.0000000            FALSE
## 70   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 71   1.026240e-01 0.077411161  1.000000e+00   1.0000000            FALSE
## 72   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 73   8.673404e-06 0.242310429  1.569886e-03   1.0000000             TRUE
## 74   5.656345e-03 0.087729816  9.898604e-01   1.0000000            FALSE
## 75   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 76   5.136878e-01 0.106136932  1.000000e+00   1.0000000            FALSE
## 77   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 78   4.725323e-01 0.577418201  1.000000e+00   1.0000000            FALSE
## 79   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 80   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 81   1.523235e-31 0.377183093  2.894147e-29   1.0000000             TRUE
## 82   3.735663e-17 0.077687151  7.060403e-15   1.0000000             TRUE
## 83   2.477228e-02 0.415963548  1.000000e+00   1.0000000            FALSE
## 84   3.846934e-01 0.196108045  1.000000e+00   1.0000000            FALSE
## 85   6.401868e-02 0.531245568  1.000000e+00   1.0000000            FALSE
## 86   6.312883e-01 0.075329211  1.000000e+00   1.0000000            FALSE
## 87   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 88   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 89   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 90   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 91   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 92   2.172750e-01 0.303958306  1.000000e+00   1.0000000            FALSE
## 93   3.848483e-17 0.070963048  7.235148e-15   1.0000000             TRUE
## 94   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 95   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 96   7.737341e-02 0.379225814  1.000000e+00   1.0000000            FALSE
## 97   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 98   6.059990e-01 0.314961900  1.000000e+00   1.0000000            FALSE
## 99   1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 100  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 101  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 102  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 103  4.235428e-04 0.092841778  7.581417e-02   1.0000000            FALSE
## 104  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 105  3.855684e-03 0.023510972  6.824560e-01   1.0000000            FALSE
## 106  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 107  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 108  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 109  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 110  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 111  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 112  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 113  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 114  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 115  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 116  7.730681e-02 0.419160255  1.000000e+00   1.0000000            FALSE
## 117  8.399698e-03 0.034283965  1.000000e+00   1.0000000            FALSE
## 118  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 119  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 120  1.143913e-31 0.190582003  2.184874e-29   1.0000000             TRUE
## 121  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 122  1.637242e-01 0.404741688  1.000000e+00   1.0000000            FALSE
## 123  4.487090e-02 0.041705450  1.000000e+00   1.0000000            FALSE
## 124  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 125  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 126  1.335774e-16 0.210397314  2.484539e-14   1.0000000             TRUE
## 127  8.095943e-02 0.581087017  1.000000e+00   1.0000000            FALSE
## 128  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 129  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 130  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 131  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 132  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 133  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 134  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 135  1.111531e-04 0.860602444  2.000756e-02   1.0000000             TRUE
## 136  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 137  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 138  2.770496e-16 0.962341130  5.097712e-14   1.0000000             TRUE
## 139  7.135747e-02 0.253759977  1.000000e+00   1.0000000            FALSE
## 140  8.830436e-17 0.179066667  1.651292e-14   1.0000000             TRUE
## 141  8.774878e-02 0.265650456  1.000000e+00   1.0000000            FALSE
## 142  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 143  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 144  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 145  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 146  3.306249e-01 0.152865759  1.000000e+00   1.0000000            FALSE
## 147  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 148  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 149  6.158060e-01 0.319607225  1.000000e+00   1.0000000            FALSE
## 150  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 151  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 152  1.212292e-02 0.742347472  1.000000e+00   1.0000000            FALSE
## 153  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 154  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 155  1.500784e-01 0.019676853  1.000000e+00   1.0000000            FALSE
## 156  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 157  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 158  3.122543e-16 0.915228428  5.683029e-14   1.0000000             TRUE
## 159  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 160  1.240446e-01 0.122678966  1.000000e+00   1.0000000            FALSE
## 161  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 162  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 163  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 164  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 165  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 166  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 167  5.143371e-02 0.642841431  1.000000e+00   1.0000000            FALSE
## 168  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 169  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 170  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 171  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 172  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 173  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 174  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 175  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 176  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 177  4.069973e-01 0.672565467  1.000000e+00   1.0000000            FALSE
## 178  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 179  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 180  2.088589e-01 0.395137454  1.000000e+00   1.0000000            FALSE
## 181  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 182  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 183  5.212093e-02 0.960171935  1.000000e+00   1.0000000            FALSE
## 184  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 185  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 186  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 187  2.927624e-16 0.571441402  5.357552e-14   1.0000000             TRUE
## 188  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 189  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 190  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
## 191  1.000000e+00 1.000000000  1.000000e+00   1.0000000            FALSE
##     diff_typeafter passed_ss_(Intercept) passed_ss_typeafter
## 1            FALSE                  TRUE                TRUE
## 2            FALSE                  TRUE                TRUE
## 3            FALSE                  TRUE                TRUE
## 4            FALSE                  TRUE                TRUE
## 5            FALSE                  TRUE                TRUE
## 6            FALSE                  TRUE                TRUE
## 7            FALSE                  TRUE                TRUE
## 8            FALSE                  TRUE                TRUE
## 9            FALSE                  TRUE                TRUE
## 10           FALSE                 FALSE                TRUE
## 11           FALSE                  TRUE                TRUE
## 12           FALSE                  TRUE                TRUE
## 13           FALSE                  TRUE                TRUE
## 14           FALSE                  TRUE                TRUE
## 15           FALSE                  TRUE                TRUE
## 16           FALSE                  TRUE                TRUE
## 17           FALSE                  TRUE                TRUE
## 18           FALSE                  TRUE               FALSE
## 19           FALSE                  TRUE               FALSE
## 20           FALSE                  TRUE                TRUE
## 21           FALSE                  TRUE                TRUE
## 22           FALSE                  TRUE                TRUE
## 23           FALSE                  TRUE                TRUE
## 24           FALSE                  TRUE               FALSE
## 25           FALSE                  TRUE                TRUE
## 26           FALSE                  TRUE                TRUE
## 27           FALSE                  TRUE                TRUE
## 28           FALSE                  TRUE                TRUE
## 29           FALSE                  TRUE                TRUE
## 30           FALSE                  TRUE                TRUE
## 31           FALSE                  TRUE                TRUE
## 32           FALSE                  TRUE                TRUE
## 33           FALSE                  TRUE                TRUE
## 34           FALSE                  TRUE                TRUE
## 35           FALSE                  TRUE                TRUE
## 36           FALSE                  TRUE                TRUE
## 37           FALSE                  TRUE                TRUE
## 38           FALSE                  TRUE                TRUE
## 39           FALSE                  TRUE                TRUE
## 40           FALSE                  TRUE               FALSE
## 41           FALSE                  TRUE                TRUE
## 42           FALSE                  TRUE                TRUE
## 43           FALSE                  TRUE                TRUE
## 44           FALSE                  TRUE                TRUE
## 45           FALSE                  TRUE                TRUE
## 46           FALSE                  TRUE                TRUE
## 47           FALSE                  TRUE                TRUE
## 48           FALSE                  TRUE                TRUE
## 49           FALSE                  TRUE                TRUE
## 50           FALSE                  TRUE                TRUE
## 51           FALSE                  TRUE                TRUE
## 52           FALSE                  TRUE                TRUE
## 53           FALSE                  TRUE                TRUE
## 54           FALSE                  TRUE               FALSE
## 55           FALSE                  TRUE                TRUE
## 56           FALSE                  TRUE                TRUE
## 57           FALSE                  TRUE                TRUE
## 58           FALSE                  TRUE               FALSE
## 59           FALSE                  TRUE                TRUE
## 60           FALSE                  TRUE                TRUE
## 61           FALSE                  TRUE                TRUE
## 62           FALSE                  TRUE                TRUE
## 63           FALSE                  TRUE                TRUE
## 64           FALSE                  TRUE                TRUE
## 65           FALSE                  TRUE                TRUE
## 66           FALSE                  TRUE                TRUE
## 67           FALSE                  TRUE                TRUE
## 68           FALSE                  TRUE                TRUE
## 69           FALSE                  TRUE                TRUE
## 70           FALSE                  TRUE                TRUE
## 71           FALSE                  TRUE                TRUE
## 72           FALSE                  TRUE                TRUE
## 73           FALSE                 FALSE                TRUE
## 74           FALSE                  TRUE                TRUE
## 75           FALSE                  TRUE                TRUE
## 76           FALSE                  TRUE                TRUE
## 77           FALSE                  TRUE               FALSE
## 78           FALSE                  TRUE                TRUE
## 79           FALSE                  TRUE                TRUE
## 80           FALSE                  TRUE                TRUE
## 81           FALSE                 FALSE                TRUE
## 82           FALSE                 FALSE                TRUE
## 83           FALSE                  TRUE                TRUE
## 84           FALSE                  TRUE                TRUE
## 85           FALSE                  TRUE                TRUE
## 86           FALSE                  TRUE                TRUE
## 87           FALSE                 FALSE               FALSE
## 88           FALSE                  TRUE                TRUE
## 89           FALSE                 FALSE               FALSE
## 90           FALSE                  TRUE                TRUE
## 91           FALSE                 FALSE               FALSE
## 92           FALSE                  TRUE                TRUE
## 93           FALSE                 FALSE                TRUE
## 94           FALSE                  TRUE                TRUE
## 95           FALSE                  TRUE                TRUE
## 96           FALSE                  TRUE                TRUE
## 97           FALSE                  TRUE                TRUE
## 98           FALSE                  TRUE                TRUE
## 99           FALSE                  TRUE                TRUE
## 100          FALSE                  TRUE                TRUE
## 101          FALSE                  TRUE                TRUE
## 102          FALSE                  TRUE                TRUE
## 103          FALSE                  TRUE                TRUE
## 104          FALSE                  TRUE               FALSE
## 105          FALSE                  TRUE                TRUE
## 106          FALSE                  TRUE                TRUE
## 107          FALSE                  TRUE                TRUE
## 108          FALSE                  TRUE                TRUE
## 109          FALSE                  TRUE                TRUE
## 110          FALSE                  TRUE                TRUE
## 111          FALSE                 FALSE               FALSE
## 112          FALSE                  TRUE                TRUE
## 113          FALSE                  TRUE                TRUE
## 114          FALSE                  TRUE               FALSE
## 115          FALSE                 FALSE               FALSE
## 116          FALSE                  TRUE                TRUE
## 117          FALSE                  TRUE                TRUE
## 118          FALSE                  TRUE                TRUE
## 119          FALSE                  TRUE               FALSE
## 120          FALSE                 FALSE                TRUE
## 121          FALSE                  TRUE                TRUE
## 122          FALSE                  TRUE                TRUE
## 123          FALSE                  TRUE                TRUE
## 124          FALSE                  TRUE                TRUE
## 125          FALSE                  TRUE                TRUE
## 126          FALSE                 FALSE                TRUE
## 127          FALSE                  TRUE                TRUE
## 128          FALSE                  TRUE               FALSE
## 129          FALSE                  TRUE                TRUE
## 130          FALSE                  TRUE                TRUE
## 131          FALSE                  TRUE                TRUE
## 132          FALSE                  TRUE               FALSE
## 133          FALSE                  TRUE                TRUE
## 134          FALSE                  TRUE                TRUE
## 135          FALSE                 FALSE                TRUE
## 136          FALSE                  TRUE                TRUE
## 137          FALSE                  TRUE                TRUE
## 138          FALSE                 FALSE                TRUE
## 139          FALSE                  TRUE                TRUE
## 140          FALSE                 FALSE                TRUE
## 141          FALSE                  TRUE                TRUE
## 142          FALSE                  TRUE                TRUE
## 143          FALSE                  TRUE                TRUE
## 144          FALSE                  TRUE                TRUE
## 145          FALSE                  TRUE                TRUE
## 146          FALSE                  TRUE                TRUE
## 147          FALSE                  TRUE                TRUE
## 148          FALSE                  TRUE                TRUE
## 149          FALSE                  TRUE                TRUE
## 150          FALSE                  TRUE                TRUE
## 151          FALSE                  TRUE                TRUE
## 152          FALSE                  TRUE                TRUE
## 153          FALSE                  TRUE                TRUE
## 154          FALSE                  TRUE                TRUE
## 155          FALSE                  TRUE                TRUE
## 156          FALSE                  TRUE                TRUE
## 157          FALSE                  TRUE                TRUE
## 158          FALSE                 FALSE                TRUE
## 159          FALSE                  TRUE                TRUE
## 160          FALSE                  TRUE                TRUE
## 161          FALSE                  TRUE                TRUE
## 162          FALSE                  TRUE                TRUE
## 163          FALSE                  TRUE                TRUE
## 164          FALSE                  TRUE                TRUE
## 165          FALSE                  TRUE                TRUE
## 166          FALSE                  TRUE                TRUE
## 167          FALSE                  TRUE                TRUE
## 168          FALSE                  TRUE                TRUE
## 169          FALSE                  TRUE                TRUE
## 170          FALSE                  TRUE                TRUE
## 171          FALSE                  TRUE                TRUE
## 172          FALSE                  TRUE                TRUE
## 173          FALSE                  TRUE                TRUE
## 174          FALSE                  TRUE                TRUE
## 175          FALSE                  TRUE                TRUE
## 176          FALSE                  TRUE                TRUE
## 177          FALSE                  TRUE                TRUE
## 178          FALSE                 FALSE               FALSE
## 179          FALSE                  TRUE                TRUE
## 180          FALSE                  TRUE                TRUE
## 181          FALSE                  TRUE                TRUE
## 182          FALSE                  TRUE                TRUE
## 183          FALSE                  TRUE                TRUE
## 184          FALSE                  TRUE                TRUE
## 185          FALSE                  TRUE                TRUE
## 186          FALSE                  TRUE                TRUE
## 187          FALSE                 FALSE                TRUE
## 188          FALSE                  TRUE                TRUE
## 189          FALSE                  TRUE                TRUE
## 190          FALSE                  TRUE                TRUE
## 191          FALSE                  TRUE                TRUE
```

``` r
write.table(res_prim_rna,"ANCOM-BC2_RNA.txt",sep="\t",col.names=NA)
# nothing significantly different after treatment relative to before (all taxa are FALSE for column "diff_typeafter")
```


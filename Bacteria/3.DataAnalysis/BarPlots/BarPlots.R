                               ##Bar Plots##

##load phyloseq object for bacteria, (BACTERIAFINALWOOD.RDS)##
BACTERIAFINALWOOD <- readRDS("../../2.PhyloseqObjectPrep/BACTERIAFINALWOOD.RDS")

##rename phyloseq object before proceeding further to Bar Plot analysis##
BACTERIAMerged <- BACTERIAFINALWOOD 

##merge samples according to GTDs condition/variety##
BACTERIAMerged <- merge_samples(BACTERIAFINALWOOD, "GTDs_condition_variety") 

##merging ccheck sample_data after merging GTDs condition/variety##
sample_data(BACTERIAMerged)

##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAMerged100 <- transform_sample_counts(BACTERIAMerged, function(OTU) 100*OTU/sum(OTU))
BACTERIAMerged100 <- prune_taxa(taxa_sums(BACTERIAMerged100)>0,BACTERIAMerged100)

##agglomerate in Phylum level##
BACTERIAMerged100GlomPhylum <-tax_glom(BACTERIAMerged100, taxrank = "Phylum")

##subset Phylum to the desirable taxa for plotting,in our case Proteobacteria,Firmicutes,Actinobacteriota,Bacteroidota,Acidobacteriota,Chloroflexi,Abditibacteriota,Verrucomicrobiota,Planctomycetota. Average > 2% calculated in Excel##
BACTERIAMerged100GlomPhylumSelection <- subset_taxa(BACTERIAMerged100GlomPhylum, Phylum %in% (c("Proteobacteria","Firmicutes","Actinobacteriota","Bacteroidota","Acidobacteriota","Chloroflexi","Abditibacteriota","Verrucomicrobiota","Planctomycetota")))

##Plot Bars for bacteria in selected taxa, first bar plot figure of the manuscript## 
plot_bar(BACTERIAMerged100GlomPhylumSelection, x="Sample", fill="Phylum", title = "") + geom_col()

##moreover prepare phyloseq object for Proteobacteria/Classes which is the second bar plot figure of manuscript##
BACTERIAMerged100Proteobacteria <- subset_taxa(BACTERIAMerged100, Phylum %in% (c("Proteobacteria")))

##remove NA, uncharacterized Classes in Proteobacteria##
BACTERIAMerged100Proteobacteria <- <- subset_taxa(BACTERIAMerged100Proteobacteria, !is.na(Phylum) & !Phylum %in% (c("", "uncharacterized")))

##Plot Bars for Proteobacteria/Classes, Alpha- and Gamma- in our study## 
plot_bar(psbacteriaWOODRawAll100Proteobacteria, x="Sample", fill="Class", title = "") + geom_col()

##PERMANOVA######
##load package vegan and adonis function##
library(vegan)

##load phyloseq object for bacteria, (BACTERIAWOOD.RDS)##
BACTERIAFINALWOOD <- readRDS("../../2.PhyloseqObjectPrep/BACTERIAFINALWOOD.RDS")

##PERMANOVA analysis must be applied to all vine varieties collectively and for each variety separately## 
##ALL VARIETIES COLLECTIVELY##
##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAFINALWOOD100 <- transform_sample_counts(BACTERIAFINALWOOD, function(OTU) 100*OTU/sum(OTU))

##All Varieties##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaBACTERIAFINAL <- adonis(BACTERIAFINALWOOD100@otu_table ~ GTDs_condition + variety + vineyard, method = "bray", data = data.frame(BACTERIAFINALWOOD100@sam_data))                              

##All Varieties##
##export the statistics## 
mypermanovaBACTERIAFINAL                              


##EACH VARIETY SEPARATELY##
##Subset DataSet to each variety and perform PERMANOVA analysis##         
##Agiorgitiko##
BACTERIAWOODAgiorgitiko<- subset_samples(BACTERIAFINALWOOD, variety=="Agiorgitiko")
BACTERIAWOODAgiorgitiko <- prune_taxa(taxa_sums(BACTERIAWOODAgiorgitiko)>0,BACTERIAWOODAgiorgitiko)

##Agiorgitiko##
##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAWOODAgiorgitiko100 <- transform_sample_counts(BACTERIAWOODAgiorgitiko, function(OTU) 100*OTU/sum(OTU))

##Agiorgitiko##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaAgiorgitiko <- adonis(BACTERIAWOODAgiorgitiko100@otu_table ~ GTDs_condition + vineyard, method = "bray", data = data.frame(BACTERIAWOODAgiorgitiko100@sam_data))
##Agiorgitiko##
##export the statistics##
mypermanovaAgiorgitiko

##Xinomavro##
BACTERIAWOODXinomavro<- subset_samples(BACTERIAFINALWOOD, variety=="Xinomavro")
BACTERIAWOODXinomavro <- prune_taxa(taxa_sums(BACTERIAWOODXinomavro)>0,BACTERIAWOODXinomavro)

##Xinomavro
##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAWOODXinomavro100 <- transform_sample_counts(BACTERIAWOODXinomavro, function(OTU) 100*OTU/sum(OTU))

##Xinomavro##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaXinomavro <- adonis(BACTERIAWOODXinomavro100@otu_table ~ GTDs_condition + vineyard, method = "bray", data = data.frame(BACTERIAWOODXinomavro100@sam_data))
##Xinomavro##
##export the statistics##
mypermanovaXinomavro

##Vidiano##
BACTERIAWOODVidiano<- subset_samples(BACTERIAFINALWOOD, variety=="Vidiano")
BACTERIAWOODVidiano <- prune_taxa(taxa_sums(BACTERIAWOODVidiano)>0,BACTERIAWOODVidiano)

##Vidiano
##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAWOODVidiano100 <- transform_sample_counts(BACTERIAWOODVidiano, function(OTU) 100*OTU/sum(OTU))

##Vidiano##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaVidiano <- adonis(BACTERIAWOODVidiano100@otu_table ~ GTDs_condition + vineyard, method = "bray", data = data.frame(BACTERIAWOODVidiano100@sam_data))
##Vidiano##
##export the statistics##
mypermanovaVidiano

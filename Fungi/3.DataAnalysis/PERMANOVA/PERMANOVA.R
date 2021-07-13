                               ##PERMANOVA######
##load package vegan and adonis function##
library(vegan)
      
##load phyloseq object for fungi, (FUNGIFINALWOOD.RDS)##
FUNGIFINALWOOD <- readRDS("../../2.PhyloseqObjectPrep/FUNGIFINALWOOD.RDS")  

##PERMANOVA analysis must be applied to all vine varieties collectively and for each variety separately## 
##ALL VARIETIES COLLECTIVELY##
##transform phyloseq object raw counts to relative abundance (100%)##
FUNGIFINALWOOD100 <- transform_sample_counts(FUNGIFINALWOOD, function(OTU) 100*OTU/sum(OTU))

##All Varieties##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaFUNGIFINAL <- adonis(FUNGIFINALWOOD100@otu_table ~ GTDs_condition + variety + vineyard, method = "bray", data = data.frame(FUNGIFINALWOOD100@sam_data))                              

##All Varieties##
##export the statistics## 
mypermanovaFUNGIFINAL                              
                               
                               
##EACH VARIETY SEPARATELY##
##Subset DataSet to each variety and perform PERMANOVA analysis##         
##Agiorgitiko##
FUNGIWOODAgiorgitiko<- subset_samples(FUNGIFINALWOOD, variety=="Agiorgitiko")
FUNGIWOODAgiorgitiko <- prune_taxa(taxa_sums(FUNGIWOODAgiorgitiko)>0,FUNGIWOODAgiorgitiko)

##Agiorgitiko##
##transform phyloseq object raw counts to relative abundance (100%)##
FUNGIWOODAgiorgitiko100 <- transform_sample_counts(FUNGIWOODAgiorgitiko, function(OTU) 100*OTU/sum(OTU))

##Agiorgitiko##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaAgiorgitiko <- adonis(FUNGIWOODAgiorgitiko100@otu_table ~ GTDs_condition + vineyard, method = "bray", data = data.frame(FUNGIWOODAgiorgitiko100@sam_data))
##Agiorgitiko##
##export the statistics##
mypermanovaAgiorgitiko

##Xinomavro##
FUNGIWOODXinomavro<- subset_samples(FUNGIFINALWOOD, variety=="Xinomavro")
FUNGIWOODXinomavro <- prune_taxa(taxa_sums(FUNGIWOODXinomavro)>0,FUNGIWOODXinomavro)

##Xinomavro
##transform phyloseq object raw counts to relative abundance (100%)##
FUNGIWOODXinomavro100 <- transform_sample_counts(FUNGIWOODXinomavro, function(OTU) 100*OTU/sum(OTU))

##Xinomavro##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaXinomavro <- adonis(FUNGIWOODXinomavro100@otu_table ~ GTDs_condition + vineyard, method = "bray", data = data.frame(FUNGIWOODXinomavro100@sam_data))
##Xinomavro##
##export the statistics##
mypermanovaXinomavro

##Vidiano##
FUNGIWOODVidiano<- subset_samples(FUNGIFINALWOOD, variety=="Vidiano")
FUNGIWOODVidiano <- prune_taxa(taxa_sums(FUNGIWOODVidiano)>0,FUNGIWOODVidiano)

##Vidiano
##transform phyloseq object raw counts to relative abundance (100%)##
FUNGIWOODVidiano100 <- transform_sample_counts(FUNGIWOODVidiano, function(OTU) 100*OTU/sum(OTU))

##Vidiano##
##PERMANOVA using Bray-Curtis dissimilarity index##
mypermanovaVidiano <- adonis(FUNGIWOODVidiano100@otu_table ~ GTDs_condition + vineyard, method = "bray", data = data.frame(FUNGIWOODVidiano100@sam_data))
##Vidiano##
##export the statistics##
mypermanovaVidiano











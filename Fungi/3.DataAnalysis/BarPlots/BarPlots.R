                              ##Bar Plots##

##load phyloseq object for fungi, (FUNGIFINALWOOD.RDS)##
FUNGIFINALWOOD <- readRDS("../../2.PhyloseqObjectPrep/FUNGIFINALWOOD.RDS")

##rename phyloseq object before proceeding further to Bar Plot analysis##
FUNGIMerged <- FUNGIFINALWOOD 

##merge samples according to GTDs condition/variety##
FUNGIMerged <- merge_samples(FUNGIFINALWOOD, "GTDs_condition_variety") 

##check sample_data after merging GTDs condition/variety##
sample_data(FUNGIMerged)

##transform phyloseq object raw counts to relative abundance (100%)##
FUNGIMerged100 <- transform_sample_counts(FUNGIMerged, function(OTU) 100*OTU/sum(OTU))
FUNGIMerged100 <- prune_taxa(taxa_sums(FUNGIMerged100)>0,FUNGIMerged100)

##agglomerate in Class level##
FUNGIMerged100GlomClass <-tax_glom(FUNGIMerged100, taxrank = "Class")

##subset Class to the desirable taxa for plotting,in our case Eurotiomycetes,Dothideomycetes,Sordariomycetes,Agaricomycetes,Leotiomycetes,Lecanoromycetes. Average > 2% calculated in Excel##
FUNGIMerged100GlomClassSelection <- subset_taxa(FUNGIMerged100GlomClass, Class %in% (c("c__Eurotiomycetes","c__Dothideomycetes","c__Sordariomycetes","c__Agaricomycetes","c__Leotiomycetes","c__Lecanoromycetes")))

##Plot Bars for fungi in selected taxa## 
plot_bar(FUNGIMerged100GlomClassSelection, x="Sample", fill="Class", title = "") + geom_col()



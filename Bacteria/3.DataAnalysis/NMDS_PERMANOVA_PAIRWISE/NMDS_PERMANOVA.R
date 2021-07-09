                            ##NMDS Analysis##

##load phyloseq object for bacteria, (BACTERIAFINALWOOD.RDS)##
##NMDS must be applied to all vine varieties collectively and for each variety separately## 

BACTERIAFINALWOOD1 <- readRDS("../../2.PhyloseqObjectPrep/BACTERIAFINALWOOD.RDS")

##ALL VARIETIES COLLECTIVELY##
##remove outlier samples *X8*,*70T*,*29T*, *13K*##
BACTERIAFINALWOOD1 <- prune_samples(!(sample_names(BACTERIAFINALWOOD) %in% (c("X8","70T","29T","13K"))), BACTERIAFINALWOOD)
BACTERIAFINALWOOD1 <- prune_taxa(taxa_sums(BACTERIAFINALWOOD1)>0,BACTERIAFINALWOOD1)

##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAFINALWOOD100 <- transform_sample_counts(BACTERIAFINALWOOD1, function(OTU) 100*OTU/sum(OTU))

##calculate bray-curtis distance##
ord.nmds.bray.BACTERIAFINALWOOD100 <- ordinate(BACTERIAFINALWOOD100, method="NMDS", distance="bray")

##plot the NMDS for all vine varieties collectively (colored per GTDs condition firstly and secondly per variety)##
plot_ordination(BACTERIAFINALWOOD100, ord.nmds.bray.BACTERIAFINALWOOD100, color="GTDs_codition",label = "NA", shape = "variety", title=paste("NMDS (stress ",round(ord.nmds.bray.BACTERIAFINALWOOD100$stress, 2),")", sep = "")) + theme_bw() + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00"))

plot_ordination(BACTERIAFINALWOOD100, ord.nmds.bray.BACTERIAFINALWOOD100, color="variety",label = "NA", shape = "GTDs_condition", title=paste("NMDS (stress ",round(ord.nmds.bray.BACTERIAFINALWOOD100$stress, 2),")", sep = "")) + theme_bw() + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00", "80b1d3"))


##EACH VARIETY SEPARATELY##
##Subset DataSet to each variety and perform NMDS only per GTDs_condition##

##Agiorgitiko##
##subset samples##
BACTERIAWOODAgiorgitiko<- subset_samples(BACTERIAFINALWOOD, variety=="Agiorgitiko")
BACTERIAWOODAgiorgitiko <- prune_taxa(taxa_sums(BACTERIAWOODAgiorgitiko)>0,BACTERIAWOODAgiorgitiko)

##Agiorgitiko
##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAWOODAgiorgitiko100 <- transform_sample_counts(BACTERIAWOODAgiorgitiko, function(OTU) 100*OTU/sum(OTU))

##Agiorgitiko
##calculate bray-curtis distance##
ord.nmds.bray.BACTERIAWOODAgiorgitiko100 <- ordinate(BACTERIAWOODAgiorgitiko100, method="NMDS", distance="bray")

##Agiorgitiko
##plot the NMDS for Agiorgitiko and color GTDs condition##
plot_ordination(BACTERIAWOODAgiorgitiko100, ord.nmds.bray.BACTERIAWOODAgiorgitiko100, color="GTDs_condition",label = "NA", shape = "vineyard", title=paste("NMDS (stress ",round(ord.nmds.bray.BACTERIAWOODAgiorgitiko100$stress, 2),")", sep = "")) + theme_bw() + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00"))

##Xinomavro##
##subset samples##
BACTERIAWOODXinomavro<- subset_samples(BACTERIAFINALWOOD, variety=="Xinomavro")
BACTERIAWOODXinomavro <- prune_taxa(taxa_sums(BACTERIAWOODXinomavro)>0,BACTERIAWOODXinomavro)

##Xinomavro##
##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAWOODXinomavro100 <- transform_sample_counts(BACTERIAWOODXinomavro, function(OTU) 100*OTU/sum(OTU))

##Xinomavro##
##calculate bray-curtis distance##
ord.nmds.bray.BACTERIAWOODXinomavro100 <- ordinate(BACTERIAWOODXinomavro100, method="NMDS", distance="bray")

##Xinomavro##
##plot the NMDS for Xinomavro and color GTDs condition##
plot_ordination(BACTERIAWOODXinomavro100, ord.nmds.bray.BACTERIAWOODXinomavro100, color="GTDs_condition",label = "NA", shape = "terroir", title=paste("NMDS (stress ",round(ord.nmds.bray.BACTERIAWOODXinomavro100$stress, 2),")", sep = "")) + theme_bw() + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00"))

##Vidiano##
##subset samples##
BACTERIAWOODVidiano <- subset_samples(BACTERIAFINALWOOD, variety=="Vidiano")
BACTERIAWOODVidiano <- prune_taxa(taxa_sums(BACTERIAWOODVidiano)>0,BACTERIAWOODVidiano)

##Vidiano##
##transform phyloseq object raw counts to relative abundance (100%)##
BACTERIAWOODVidiano100 <- transform_sample_counts(BACTERIAWOODVidiano, function(OTU) 100*OTU/sum(OTU))

##Vidiano##
##calculate bray-curtis distance##
ord.nmds.bray.BACTERIAWOODVidiano100 <- ordinate(BACTERIAWOODVidiano100, method="NMDS", distance="bray")

##plot the NMDS for Vidiano and color GTDs condition##
plot_ordination(BACTERIAWOODVidiano100, ord.nmds.bray.BACTERIAWOODVidiano100, color="GTDs_condition",label = "NA", shape = "terroir", title=paste("NMDS (stress ",round(ord.nmds.bray.BACTERIAWOODVidiano100$stress, 2),")", sep = "")) + theme_bw() + geom_point(size = 4) + scale_color_manual(values=c("#984EA3","#FF7F00"))


##pairwise permanova to groups of NMDS###

##check if pairwise permanova between variables (variety, GTDs condition, vineyard) differ significally and use p-value to manuscript data##

##load package pairwiseAdonis##
library(pairwiseAdonis)

##ALL VARIETIES COLLECTIVELY##
##variety vs GTDs condition##
mycmpfactorbacteriaAll <- interaction(data.frame(BACTERIAFINALWOOD100@sam_data)$variety, data.frame(BACTERIAFINALWOOD100@sam_data)$GTDs_condition)

mympairwisepermbacteriaAll <- pairwise.adonis(BACTERIAFINALWOOD100@otu_table, mycmpfactorfungiAll, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Agiorgitiko##
##variety vs GTDs condition##
mycmpfactorAgiorgitiko <- interaction(data.frame(BACTERIAWOODAgiorgitiko100@sam_data)$variety, data.frame(BACTERIAWOODAgiorgitiko100@sam_data)$GTDs_condition)

mympairwisepermAgiorgitiko <- pairwise.adonis(BACTERIAWOODAgiorgitiko100@otu_table, mycmpfactorAgiorgitiko, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Agiorgitiko##
##vineyard vs GTDs condition##
mycmpfactorAgiorgitiko2 <- interaction(data.frame(BACTERIAWOODAgiorgitiko100@sam_data)$GTDs_condition, data.frame(BACTERIAWOODAgiorgitiko100@sam_data)$vineyard)

mympairwisepermAgiorgitiko2 <- pairwise.adonis(BACTERIAWOODAgiorgitiko100@otu_table, mycmpfactorAgiorgitiko2, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Agiorgitiko##
##vineyard vs variety###
mycmpfactorAgiorgitiko3 <- interaction(data.frame(BACTERIAWOODAgiorgitiko100@sam_data)$variety, data.frame(BACTERIAWOODAgiorgitiko100@sam_data)$vineyard)

mympairwisepermAgiorgitiko3 <- pairwise.adonis(BACTERIAWOODAgiorgitiko100@otu_table, mycmpfactorAgiorgitiko3, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Xinomavro##
##variety vs GTDs condition##
mycmpfactorXinomavro <- interaction(data.frame(BACTERIAWOODXinomavro100@sam_data)$variety, data.frame(BACTERIAWOODXinomavro100@sam_data)$GTDs_condition)

mympairwisepermXinomavro <- pairwise.adonis(BACTERIAWOODXinomavro100@otu_table, mycmpfactorXinomavro, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Xinomavro##
##vineyard vs GTDs condition##
mycmpfactorXinomavro2 <- interaction(data.frame(BACTERIAWOODXinomavro100@sam_data)$GTDs_condition, data.frame(BACTERIAWOODXinomavro100@sam_data)$vineyard)

mympairwisepermXinomavro2 <- pairwise.adonis(BACTERIAWOODXinomavro100@otu_table, mycmpfactorXinomavro2, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Xinomavro##
##vineyard vs variety###
mycmpfactorXinomavro3 <- interaction(data.frame(BACTERIAWOODXinomavro100@sam_data)$variety, data.frame(BACTERIAWOODXinomavro100@sam_data)$vineyard)

mympairwisepermXinomavro3 <- pairwise.adonis(BACTERIAWOODXinomavro100@otu_table, mycmpfactorXinomavro3, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Vidiano##
##Variety vs GTDs condition##
mycmpfactorVidiano <- interaction(data.frame(BACTERIAWOODVidiano100@sam_data)$variety, data.frame(BACTERIAWOODVidiano100@sam_data)$GTDs_condition)

mympairwisepermVidiano <- pairwise.adonis(BACTERIAWOODVidiano100@otu_table, mycmpfactorVidiano, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Vidiano##
##vineyard vs GTDs condition##
mycmpfactorVidiano2 <- interaction(data.frame(BACTERIAWOODVidiano100@sam_data)$GTDs_condition, data.frame(BACTERIAWOODVidiano100@sam_data)$vineyard)

mympairwisepermVidiano2 <- pairwise.adonis(BACTERIAWOODVidiano100@otu_table, mycmpfactorVidiano2, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

##Vidiano##
##vineyard vs variety###
mycmpfactorVidiano3 <- interaction(data.frame(BACTERIAWOODVidiano100@sam_data)$variety, data.frame(BACTERIAWOODVidiano100@sam_data)$vineyard)

mympairwisepermVidiano3 <- pairwise.adonis(BACTERIAWOODVidiano100@otu_table, mycmpfactorVidiano3, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")


##### Perform the analysis with the transformed abundance matrix and agglomerate to genus ----
library(pime)
library(phyloseq)
ps_orig_pre_pre <- readRDS("../../2.PhyloseqObjectPrep/FUNGIFINALWOOD.RDS")
ps_orig_pre <- ps_orig_pre_pre
otu_table(ps_orig_pre) <- t(otu_table(ps_orig_pre_pre))

# ps_orig <- tax_glom(ps_orig_pre, "Genus")
ps_orig <- ps_orig_pre

#ps_tr <- transform_sample_counts(ps_orig, fun = function(x) {100*x/sum(x)})
#ps_tr <- microbial::normalize(ps_orig, "GTDs_condition", method = "TMM")
set.seed(42)
ps_tr <- rarefy_even_depth(ps_orig, sample.size = summary(sample_sums(ps_orig))[2])
#ps_tr <- rarefy_even_depth(ps_orig, sample.size = min(sample_sums(ps_orig)))

# calculates the overall OOB error of the original dataset (or whatever dataset provided)
ps_orig.oob.error <- pime.oob.error(ps_tr, "GTDs_condition")

# splits the object by treatment level and chooses features that best fit the treatment levels
per_variable_obj <- pime.split.by.variable(ps_tr, "GTDs_condition")

# the following creates a list of the phyloseq objects according to the selection performed previously
prevalences <- pime.prevalence(per_variable_obj)

# set.seed(42)

# calculates the OOB error for choosing the best feature set
best.prev <- pime.best.prevalence(prevalences, "GTDs_condition")


imp <- best.prev$`Importance`$`Prevalence 40`



##### perform the second selection layer with the Kruskal-Wallis test ----

my_ps_sel <- prevalences$`40`

## prepare the percentage table for plotting
ps_orig_ra <- transform_sample_counts(ps_orig, fun = function(x) {100*x/sum(x)})
# keep only the selected taxa
my_ps_sel_ra <- prune_taxa(taxa_names(my_ps_sel),ps_orig_ra)
my_ps_sel_ra <- prune_samples(sample_names(my_ps_sel),my_ps_sel_ra)
my_ps_sel_ra <- prune_taxa(taxa_sums(my_ps_sel_ra) > 0, my_ps_sel_ra) 

otu_table_ps_sel_ra <- data.frame(otu_table(my_ps_sel_ra))

myotus <- row.names(otu_table_ps_sel_ra)

# perform the Kruskal analysis
mystatsout_all <- list()
library(agricolae)
for(myotu in myotus){
  mykrusk <- kruskal(t(otu_table_ps_sel_ra)[,myotu], data.frame(sample_data(my_ps_sel_ra))$GTDs_condition, group = T)
  mystatsout_all[[myotu]][["krusk"]] <- mykrusk
}



mylett_tbl <- array(dim = c(nrow(otu_table_ps_sel_ra),length(levels(factor(data.frame(sample_data(my_ps_sel_ra))$GTDs_condition))) + 1))
row.names(mylett_tbl) <- row.names(otu_table_ps_sel_ra)
colnames(mylett_tbl) <- c(levels(factor(data.frame(sample_data(my_ps_sel_ra))$GTDs_condition)),"Kruskal_Wallis_P_value")

for(myotu in myotus){
  if(all(mystatsout_all[[myotu]]$krusk$groups[colnames(mylett_tbl)[-grep("Kruskal",colnames(mylett_tbl))],"groups"] == "a")){
    mylett_tbl[myotu,-grep("Kruskal",colnames(mylett_tbl))] <- ""
    mylett_tbl[myotu,"Kruskal_Wallis_P_value"] <- mystatsout_all[[myotu]]$krusk$statistics$p.chisq
  } else {
    mylett_tbl[myotu,-grep("Kruskal",colnames(mylett_tbl))] <- mystatsout_all[[myotu]]$krusk$groups[colnames(mylett_tbl)[1:length(levels(factor(data.frame(sample_data(my_ps_sel_ra))$GTDs_condition)))],"groups"]
    mylett_tbl[myotu,"Kruskal_Wallis_P_value"] <- mystatsout_all[[myotu]]$krusk$statistics$p.chisq
  }
}

mylett_tbl_df <- as.data.frame(mylett_tbl, stringsAsFactors = F)
mylett_tbl_df$Kruskal_Wallis_P_value <- as.numeric(mylett_tbl_df$Kruskal_Wallis_P_value)

mylett_tbl_df_sig <- mylett_tbl_df[which(mylett_tbl_df$Kruskal_Wallis_P_value <= 0.05),-grep("Kruskal",colnames(mylett_tbl_df))]

## obtain the significant OTUs
significant_OTUs <- row.names(mylett_tbl_df_sig)
#### prepare the bar plots
# mystatsout_all
# mylett_tbl
# mytptgsTblOrdrd_sig

library(gtools)
# prepare the barplot with the habitat and time index
cairo_pdf(paste("DA.pdf",sep=""), height = 22, width = 20, onefile = T)
par(mfrow = c(11,5))
for(myotu_sig in row.names(mylett_tbl_df_sig)){
  mytestvarord <- colnames(mylett_tbl_df_sig)
  par(mar = c(4,10,4,8))
  par(xpd = T)
  
  barerrplt <- bar.err(mystatsout_all[[myotu_sig]]$krusk$means[mytestvarord[length(mytestvarord):1],], variation="SE", xlim=c(0,max(mystatsout_all[[myotu_sig]]$krusk$means[mytestvarord[length(mytestvarord):1],1] + mystatsout_all[[myotu_sig]]$krusk$means[mytestvarord[length(mytestvarord):1],3])),horiz=TRUE, bar=T, col="grey60", space=0.5
                       , main= paste(myotu_sig," ",data.frame(tax_table(my_ps_sel_ra))[myotu_sig,]$Phylum,": ",data.frame(tax_table(my_ps_sel_ra))[myotu_sig,]$Genus," ",  stars.pval(mystatsout_all[[myotu_sig]]$krusk$statistics$p.chisq), "\n", "P value ", round(mystatsout_all[[myotu_sig]]$krusk$statistics$p.chisq,3), sep = "")
                       , las=1, xlab = "%", cex.lab = 1.5, font.lab = 2)
  
  text(x = mystatsout_all[[myotu_sig]]$krusk$means[mytestvarord[length(mytestvarord):1],1] + mystatsout_all[[myotu_sig]]$krusk$means[mytestvarord[length(mytestvarord):1],3]/sqrt(mystatsout_all[[myotu_sig]]$krusk$means[mytestvarord[length(mytestvarord):1],4] - 1), barerrplt$x,labels = mylett_tbl_df_sig[myotu_sig,ncol(mylett_tbl_df_sig):1], pos = 4, font = 3)
  par(xpd = F)
  
}
dev.off()

#### proceed with deeper analysis ----

# according to the best.prev$`OOB error` output, prevalence 30% provided the lowest OOB error with the most Nseqs... so lets choose this for the followup analysis

# ps_sel_ra <- prevalences$`40`


# run the random forests with the OTUs passing both above tests (random forest selection and Kruskal-Wallis test)
library(randomForest)
ps_choice <- prune_taxa(significant_OTUs,my_ps_sel)
ps_choice <- prune_taxa(taxa_sums(ps_choice) > 0, ps_choice)
ps_choice <- prune_samples(sample_sums(ps_choice) > 0, ps_choice)

mytbl <- t(data.frame(otu_table(ps_choice)))

myntree = 2*round((sqrt(dim(mytbl)[1] * dim(mytbl)[2])/4)/2) + 1
# if these are too low you kan keep a minimum of 1001
if(myntree < 1001){
  myntree_fin <- 1001 
} else {
  myntree_fin <- myntree
}

# the nPerm below should not be effective in the classification forests
# run the random forests
library(randomForest)
RF_GTDs_condition_classify <- randomForest(x=mytbl, y=factor(sample_data(ps_choice)$GTDs_condition), ntree = myntree_fin, importance=TRUE, proximity=TRUE, nPerm = 1000)


# prepare some plots
mydist <- as.dist(RF_GTDs_condition_classify$proximity)

library(ape)

myPCoA <- pcoa(mydist)

plot(myPCoA$vectors[,1:2], type = "p", pch = 21, bg = as.factor(sample_data(ps_choice)$GTDs_condition), bty = "n")
library(vegan)
ordiellipse(myPCoA$vectors[,1:2], groups = as.factor(sample_data(ps_choice)$GTDs_condition), kind = "ehull")
#ordiellipse(myPCoA$vectors[row.names(sample_data(ps_sel_ra)),1:2], groups = as.factor(sample_data(ps_sel_ra)$GTDs_condition), kind = "sd")

### run the RDA or CCA
# run DCA to check if RDA or CCA is more appropriate
ps_sel_ra.ord <- ordinate(ps_choice)
## result
# Call:
#   decorana(veg = veganifyOTU(physeq)) 
# 
# Detrended correspondence analysis with 26 segments.
# Rescaling of axes with 4 iterations.
# 
# DCA1   DCA2   DCA3   DCA4
# Eigenvalues     0.5280 0.2735 0.1618 0.1306
# Decorana values 0.5314 0.2689 0.1522 0.1197
# Axis lengths    2.6880 4.0048 3.4760 2.5982
## CCA seems a legitimate choice

# perform a PERMANOVA whose outcomes to add on the plot title
mypermanova2 <- adonis2(phyloseq::distance(ps_choice,method = "bray") ~  GTDs_condition, data.frame(sample_data(ps_choice)), permutations = 999,
                        sqrt.dist = FALSE, add = FALSE, by = NULL,
                        parallel = 4)

ps_choice.rda <- ordinate(ps_choice, method = "RDA", formula = ~GTDs_condition)
library(ggplot2)
p1 = plot_ordination(ps_choice, ps_choice.rda, type="samples", color="GTDs_condition", title= paste("model: comm ~ ","GTDs_condition", "; Rsq = ",round(100*mypermanova2$R2,2),"% , P = ",round(mypermanova2$`Pr(>F)`,3),sep = ""))
p2 <- p1 + geom_point(size=2) + scale_colour_brewer(type="qual", palette="Set1") + theme_minimal() +stat_ellipse(type = "t", linetype = 2)
cairo_pdf("RDA_rarefied_at_1st_quratile_40perc_prev_genus_glom.pdf", height = 5, width = 8)
print(p2)
dev.off()


### do the same as above for the non canonical
mypermanova2 <- adonis2(phyloseq::distance(ps_choice,method = "bray") ~  GTDs_condition, data.frame(sample_data(ps_choice)), permutations = 999,
                        sqrt.dist = FALSE, add = FALSE, by = NULL,
                        parallel = 4)

ps_choice.nmds <- ordinate(ps_choice, method = "NMDS")
library(ggplot2)
p1 = plot_ordination(ps_choice, ps_choice.nmds, type="samples", color="GTDs_condition", title= paste("model: comm ~ ","GTDs_condition", "; Rsq = ",round(100*mypermanova2$R2,2),"% , P = ",round(mypermanova2$`Pr(>F)`,3), "; Stress ", round(ps_choice.nmds$stress,3),sep = ""))
p2 <- p1 + geom_point(size=2) + scale_colour_brewer(type="qual", palette="Set1") + theme_minimal() +stat_ellipse(type = "t", linetype = 2)
cairo_pdf("NMDS_rarefied_at_1st_quratile_40perc_prev_genus_glom.pdf", height = 5, width = 8)
print(p2)
dev.off()



# #### permanova
# mytbl <- data.frame(otu_table(ps_sel_ra), stringsAsFactors = F)
# mydesign <- data.frame(sample_data(ps_sel_ra))
# library(vegan)
# mypermanova <- adonis2(mytbl ~  pest_grp + dose, mydesign, permutations = 999, method = "bray",
#                        sqrt.dist = FALSE, add = FALSE, by = "terms",
#                        parallel = 4)
# mypermanova2 <- adonis2(mytbl ~  pest_grp + dose, mydesign, permutations = 999, method = "bray",
#                         sqrt.dist = FALSE, add = FALSE, by = NULL,
#                         parallel = 4)
# write.table(data.frame(mypermanova, check.names = F), file = paste("PERMANOVA.txt",sep=""), quote = F, col.names = NA)





#### pairwise PERMANOVA analysis between reactors (before and after treatment grouped bioreactors) ----
library(pairwiseAdonis)


pariwisecmps <- pairwise.adonis(t(data.frame(otu_table(ps_choice))), data.frame(sample_data(ps_choice))$GTDs_condition, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "BH")

### prepare a heatmap with F

pariwisecmps_F <- data.frame(do.call("rbind",strsplit(as.character(pariwisecmps$pairs)," vs ", fixed = T)), F = pariwisecmps$F.Model)
pariwisecmps_P <- data.frame(do.call("rbind",strsplit(as.character(pariwisecmps$pairs)," vs ", fixed = T)), F = pariwisecmps$p.value)

# convert to wide
library(spaa)
mypairwisedist <- list2dist(pariwisecmps_F)
as.matrix(mypairwisedist)
mypairwisepval <- list2dist(pariwisecmps_P)
as.matrix(mypairwisepval)

# create the plot

mypairwisedistfrlpt <- as.matrix(mypairwisedist)
mypairwisepvalfrlpt <- as.matrix(mypairwisepval)
mypairwisepvalfrlpt[is.na(mypairwisepvalfrlpt)] <- 1

mybreaks <- seq(0,max(mypairwisedistfrlpt, na.rm = T), length.out = 120)
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors

cairo_pdf(paste("PERMANOVA_pairwise_rarefied_at_1stQ_40perc_prev_genus_glom.pdf",sep=""), height = 3, width = 4)
pheatmap::pheatmap(mypairwisedistfrlpt, color = mycolors, breaks = mybreaks, cellwidth = 15, cellheight = 15, clustering_method = "complete"
                   #, cluster_rows = F, cluster_cols = F
                   , treeheight_row = 10
                   , treeheight_col = 10
                   , display_numbers = round(as.matrix(mypairwisepvalfrlpt), 2)
                   , fontsize_number = 6
                   , number_color = "black"
)
dev.off()



### prepare the taxa heatmap

sel_import_pre <- RF_GTDs_condition_classify$importance[-which(RF_GTDs_condition_classify$importance[,"MeanDecreaseAccuracy"] <= 0),-grep("Mean",colnames(RF_GTDs_condition_classify$importance))]
sel_import <- vegan::decostand(sel_import_pre,method = "range", MARGIN = 1)

my_tax_names <- data.frame(tax_table(ps_choice))
my_tax_names$name <- paste(row.names(my_tax_names),paste(my_tax_names$Phylum,":",my_tax_names$Genus, sep = ""))

row.names(sel_import) <- my_tax_names[row.names(sel_import),"name"]

# prep also the relative abundance column
ps_tr_ra <- transform_sample_counts(ps_tr, fun = function(x) 100*x/sum(x))
ps_tr_ra_otu_tbl <- data.frame(otu_table(ps_tr_ra), check.names = F)
ps_tr_ra_otu_tbl_sel_max <- apply(ps_tr_ra_otu_tbl[gsub(" .+","",row.names(sel_import)),], MARGIN = 1, max)
ps_tr_ra_otu_tbl_sel_max_df <- data.frame(`max RA %`=ps_tr_ra_otu_tbl_sel_max, check.names = F)
ps_tr_ra_otu_tbl_sel_max_df <- ps_tr_ra_otu_tbl_sel_max_df[gsub(" .+","",row.names(sel_import)),, drop = F]
row.names(ps_tr_ra_otu_tbl_sel_max_df) <- row.names(sel_import)

# mybreaks <- seq(-max(sel_import, na.rm = T),max(sel_import, na.rm = T), length.out = 120)
mybreaks <- seq(0,1, length.out = 120)
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors

cairo_pdf(paste("taxa_rarefied_at_1stQ_40perc_prev_genus_glom.pdf",sep=""), height = 10, width = 12)
pheatmap::pheatmap(sel_import, color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
                   #, cluster_rows = F, cluster_cols = F
                   , treeheight_row = 10
                   , treeheight_col = 10
                   # , display_numbers = round(as.matrix(sel_import), 2)
                   , fontsize_number = 2
                   , number_color = "black"
                   , annotation_row = ps_tr_ra_otu_tbl_sel_max_df
)
dev.off()






test <- ps_tr_ra_otu_tbl[gsub(" .+","",row.names(sel_import)),]
test <- test[gsub(" .+","",row.names(sel_import)),]
row.names(test) <- row.names(sel_import)
test2 <- decostand(test,"range", MARGIN = 1)
mybreaks <- seq(0,1, length.out = 120)
mycolors <- colorRampPalette(RColorBrewer::brewer.pal(n = 6, name = 'RdBu')[6:1])(n = 119) # the nature colors
cairo_pdf(paste("taxa_rarefied_at_1stQ_40perc_prev_genus_glom_all_samples_and_RA.pdf",sep=""), height = 10, width = 15)
pheatmap::pheatmap(test2, color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
                   #, cluster_rows = F, cluster_cols = F
                   , treeheight_row = 10
                   , treeheight_col = 10
                   # , display_numbers = round(as.matrix(sel_import), 2)
                   , fontsize_number = 2
                   , number_color = "black"
                   , annotation_row = ps_tr_ra_otu_tbl_sel_max_df
)
dev.off()

### prepare the core microbiome heatmap
test3 <- stats::aggregate(t(test2[,row.names(sample_data(ps_tr_ra))]), by = list(GTDs_condition = sample_data(ps_tr_ra)$GTDs_condition), mean)

row.names(test3) <- test3$GTDs_condition
test3 <- test3[,-grep("GTDs_condition",colnames(test3))]
test4 <- decostand(t(test3), method = "range", MARGIN = 1)
cairo_pdf(paste("taxa_rarefied_at_1stQ_40perc_prev_genus_glom_all_samples_and_RA_core_microb.pdf",sep=""), height = 10, width = 15)
pheatmap::pheatmap(as.matrix(test4), color = mycolors, breaks = mybreaks, cellwidth = 10, cellheight = 10, clustering_method = "complete"
                   #, cluster_rows = F, cluster_cols = F
                   , treeheight_row = 10
                   , treeheight_col = 10
                   # , display_numbers = round(as.matrix(sel_import), 2)
                   , fontsize_number = 2
                   , number_color = "black"
                   , annotation_row = ps_tr_ra_otu_tbl_sel_max_df[row.names(test4),, drop = F]
)
dev.off()

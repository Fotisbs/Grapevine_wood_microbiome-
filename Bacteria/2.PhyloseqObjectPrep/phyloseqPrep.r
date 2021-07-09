##Quality Control, classification and phyloseq object construction##
##install the necessary packages##
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")

##load the package##
library(dada2); packageVersion("dada2")

mymultthread <- 56

##input file path and storing in variables##
##change to the directory containing the fastq files after unzipping##
path1 <- ("../1.Demultiplex/demux_out/analysis_ready")

##lists files in a path##
list.files(path1)

##set the variables containing all the forward and the reverse paths to the files of interest with the list.files command##
fnFs <- sort(list.files(path1, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path1, pattern="_R2.fastq", full.names = TRUE))

##Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq##
sample.names <- gsub("_.+","",basename(fnFs))

##Plot the per-base qualities## 
plotQualityProfile(c(fnFs[1],fnRs[1]))

pdf("Quality_raw_seqs.pdf",onefile=T)
for (i in c(1:length(fnFs))) {
  print(i)
  plot1 <- plotQualityProfile(fnFs[i])
  print(plot1)
  plot2 <- plotQualityProfile(fnRs[i])
  print(plot2)
}
dev.off() 

##sequence quality filtering and control, error modelling, and dereplication##
##set the file paths where the quality controlled sequences are going to be saved##
filtFs <- file.path("filtered_seqs", paste(sample.names, "_F_filt.fastq.gz", sep = ""))
filtRs <- file.path("filtered_seqs", paste(sample.names, "_R_filt.fastq.gz", sep = ""))

##filter the sequences and save then in the folders provided above and get their statistics in a table##
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 11,
                     maxN=0, maxEE=c(2,2), truncLen = c(230,230), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=mymultthread, matchIDs=FALSE)
head(out)

##learns error rates using a machine learning algorithm##
errF <- learnErrors(filtFs, multithread=mymultthread)
plotErrors(errF, nominalQ=TRUE)
errR <- learnErrors(filtRs, multithread=mymultthread)
plotErrors(errR, nominalQ=TRUE)

##dereplication of each one of the red pairs to unique sequences (collapsing of the identical sequences for each pair per sample)##
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

##Name the derep-class objects by the sample names##
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##sample composition inference after read correction##
dadaFs <- dada(derepFs, err=errF, multithread=mymultthread)
dadaRs <- dada(derepRs, err=errR, multithread=mymultthread)

##merge read pairs retaining the per sequence sample information##
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

##construct the sequence table, remove the chimeras, and create a summary##
##construct sequence table##
seqtab <- makeSequenceTable(mergers)

##chimera removal (consensus)##
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=mymultthread, verbose=TRUE)

##dimensions of seqtab.nochim##
dim(seqtab.nochim)

##record the portion of good sequences out of the total prior the chimera removal##
sum(seqtab.nochim)/sum(seqtab)  

##track reads through the pipeline##
getN <- function(x) {sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

head(track)

##save the read quality control data##
write.table(track, file = "readQC_bac_wd.txt", col.names = NA, sep = "\t", quote = FALSE)

##taxonomically classify the sequences##
taxa <- assignTaxonomy(seqtab.nochim, minBoot = 80, "silva_nr99_v138.1_train_set.fa.gz", multithread=mymultthread)

##go as low as species where possible with the species alignment##
taxa <- addSpecies(taxa, "silva_nr99_v138.1_wSpecies_train_set.fa.gz")

##proceed with analysis, load the phyloseq package##
library(phyloseq); packageVersion("phyloseq")

##load the ggplot2 package, for plotting functions##
# load the experimental design data##
samdf <- read.table(file ="design", header = T, row.names = 1, sep = "\t", check.names = F, quote = "", comment.char = "")

##construct the phyloseq object## 
ps_bac <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE) 
                   ,sample_data(samdf) 
                   ,tax_table(taxa))

##replace the taxon names (the sequences of each ASV with something easier to read and save the sequence info in the object)##
##load the Biostrings and stringr packages##
library("Biostrings")

sequences <- Biostrings::DNAStringSet(taxa_names(ps_bac))
names(sequences) <- taxa_names(ps_bac)

ps <- merge_phyloseq(ps_bac, sequences)

ps_bac_wd <- ps

library(stringr)
taxa_names(ps_bac_wd) <- paste("ASV",str_pad(1:length(taxa_names(ps_bac_wd)),4, pad = "0"),sep = "")

##Create table,and check number of features for each phyla##
table(tax_table(ps_bac_wd)[, "Kingdom"], exclude = NULL)
table(tax_table(ps_bac_wd)[, "Phylum"], exclude = NULL) 
table(tax_table(ps_bac_wd)[, "Class"], exclude = NULL) 
table(tax_table(ps_bac_wd)[, "Order"], exclude = NULL)
table(tax_table(ps_bac_wd)[, "Family"], exclude = NULL)
table(tax_table(ps_bac_wd)[, "Genus"], exclude = NULL)

##Further, features with ambiguous annotation, low abundance or non target taxa removed from phyloseq object##
##For our analysis removed Eukaryota, Archaea, NA##
ps_bac_wd <- subset_taxa(ps_bac_wd, !is.na(Kingdom) & !Kingdom %in% (c("", "uncharacterized", "Eukaryota","Archaea")))

##also remove mitochondria and chloroplasts##

ps_bac_wd <- subset_taxa(ps_bac_wd, !Order %in% c("Chloroplast"))

ps_bac_wd <- subset_taxa(ps_bac_wd, !Family %in% c("Mitochondria"))

ps_bac_wd <- prune_taxa(taxa_sums(ps_bac_wd)>0, ps_bac_wd)

##replace the NA names with the classified leftmost (higher level annotated) taxa##
ps_bac_wd <- BACTERIAFINALWOOD

for(i in 1:nrow(tax_table(BACTERIAFINALWOOD))){
  for(j in 2:ncol(tax_table(BACTERIAFINALWOOD))){
    if(is.na(tax_table(BACTERIAFINALWOOD)[i,j])){
      tax_table(BACTERIAFINALWOOD)[i,j] <- tax_table(BACTERIAFINALWOOD)[i,j-1]
    }
  }
}

##save phyloseq object and continue with the statistical analysis scripts##

saveRDS(BACTERIAFINALWOOD, file = "BACTERIAFINALWOOD.RDS")
# ***Grapevine wood microbiome analysis identifies key fungal pathogens and potential interactions with the bacterial community implicated in grapevine trunk dise ase appearance***
### By Bekris F. <sup>1</sup>, Vasileiadis S. <sup>1</sup>, Papadopoulou E. <sup>1</sup>, Samaras A. <sup>2</sup>, Testempasis S. <sup>2</sup>, Gkizi D. <sup>3</sup>, Tavlaki, G. <sup>4</sup>, Tzima A. <sup>3</sup>, Paplomatas E. <sup>3</sup>, Markakis E. <sup>4</sup>, Karaoglanidis G. <sup>2</sup>, Papadopoulou K.K. <sup>1</sup>, Karpouzas D.G. <sup>1*</sup>

### (\* corr. author)


<sup>1</sup> University of Thessaly, Department of Biochemistry and Biotechnology, Laboratory of Plant and Environmental Biotechnology, Viopolis – 41500 Larissa, Greece

<sup>2</sup> Plant Pathology Laboratory, Faculty of Agriculture, Aristotle University of Thessaloniki, Thessaloniki, Greece

<sup>3</sup> Laboratory of Plant Pathology, Agricultural University of Athens, Iera Odos 75, 11855 Athens, Greece

<sup>4</sup> Laboratory of Mycology, Department of Viticulture, Vegetable Crops, Floriculture and Plant Protection, Institute of Olive Tree, Subtropical Crops and Viticulture, Hellenic Agricultural Organization DIMITRA, 32A Kastorias street, Mesa Katsabas 71307, Heraklion, Crete, Greece



## The provided material includes the code used in the statistical analysis of the study.

For obtaining the code the users need to open a terminal and having the [GitHub tools](https://github.com/git-guides/install-git), git-clone or download the repository, and enter the base folder. E.g:

```
$ git clone https://github.com/Fotisbs/Grapevine_wood_microbiome-.git
```

In the case of the computational methods, with the "Grapevine_wood_microbiome-" folder as working directory, and assuming that the necessary software and R packages are installed, the used code can be executed as described in this Readme.md file. The necessary datasets for performing all sequencing based analysis can be downloaded implementing the code provided in the corresponding repository folders as explained below.

## Description of the order of executed scripts.

Steps 1-3 concern the data retrieval from NCBI and preprocessing, while steps 4-6 concern the actual data analysis for total fungi and bacteria. 

1) First, it is necessary to download the sequencing data.
To do so, you need to enter the "0.DownloadData" subfolder of e.g. the "Fungi" and execute the "fetch_data.sh" bash script (this assumes that you are located at the working directory "Grapevine_wood_microbiome-").
The script is based on the SRR accession numbers found in the 0.DownloadData folder.
Once the download is done, you need to combine all forward reads to a single file and all reverse reads to another file as well.
```
for i in {1..3}
do
	cd Fungi/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
	cd Bacteria/0.DownloadData/batch${i}
	sh -x fetch_data.sh
	cat *_1.fastq | gzip > forward.fastq.gz
	cat *_2.fastq | gzip > reverse.fastq.gz
	cd ../../../
done
```

2) Then you need to demultiplex the data according to our own demultiplexing method using our in-house script.
This requires Flexbar v3.0.3 to be installed as described in the manuscript.
A detailed description of our in-house multiplexing approach is provided in our [previous work] (https://github.com/SotiriosVasileiadis/mconsort_tbz_degr#16s).
You need to enter the folder Fungi/1.Demultiplex and run the following commands (change the MY_PROCS variable to whatever number of logical processors you have available and want to devote).
the following commands are going to save the demultiplexed files in the Fungi(or Bacteria)/1.Demultiplex/demux_out folder.
```
MY_WORKING_DIR_BASE=`pwd`
for i in {1..3}
do
  cd Fungi/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_out${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/reverse.fastq.gz fun${i}_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
  cd Bacteria/1.Demultiplex
  MY_PROCS=56
  bash DemuxOwnBCsys_absPATH.sh demux_ou${i} ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/forward.fastq.gz ${MY_WORKING_DIR_BASE}/Fungi/0.DownloadData/batch${i}/reverse.fastq.gz bac${i}_map_file.txt ${MY_PROCS}
  cd demux_out${i}/analysis_ready
  gunzip *.gz # unzips files skipped by the Demux script
  cd ../../../../
done

cd Fungi/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../

cd Bacteria/1.Demultiplex
mkdir -p demux_out/analysis_ready
cp demux_out[0-9]/analysis_ready/*.fastq demux_out/analysis_ready/
cd ../../
```
3) Following, the "phyloseqPrep.r" script of the Fungi(or Bacteria)/2.PhyloseqObjectPerp folder is run in order to prepare the final phyloseq object to be used in the data analysis described below. Before runnin gthe script make sure that the necessary reference databases are found in the same folder.
```
cd Fungi/2.PhyloseqObjectPrep
# fetch the databases
wget https://files.plutof.ut.ee/public/orig/1D/B9/1DB95C8AC0A80108BECAF1162D761A8D379AF43E2A4295A3EF353DD1632B645B.gz
# run the R script
Rscript phyloseqPrep.r
cd ../../
cd Bacteria/2.PhyloseqObjectPrep
# fetch the databases
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz
tar vxf *.gz
# run the R script
Rscript phyloseqPrep.r
cd ../../
```
4a) Run the overall PERMANOVA tests.
```
cd Fungi/3.DataAnalysis/PERMANOVA
Rscript PERMANOVA.R
cd ../../../
cd Bacteria/3.DataAnalysis/PERMANOVA
Rscript PERMANOVA.R
cd ../../../
```
4b) Run the NMDS and PERMANOVA tests.
```
cd Fungi/3.DataAnalysis/NMDS_PERMANOVA_PAIRWISE
Rscript NMDS_PERMANOVA.R
cd ../../../
cd Bacteria/3.DataAnalysis/NMDS_PERMANOVA_PAIRWISE
Rscript NMDS_PERMANOVA.R
cd ../../../
```
4c) Prepare the barplots.
```
cd Fungi/3.DataAnalysis/BarPlots
Rscript BarPlots.R
cd ../../../
cd Bacteria/3.DataAnalysis/BarPlots
Rscript BarPlots.R
cd ../../../
```

5) Pathogenic fungi, linked to the grapevine trunk decline (GTD) complex, were selected as described in the manuscript.

6) Network analysis was performed no attempt and identify links among the GTD complex members and between the complex and the rest fungi and total bacteria. run the following

```
cd NetWork
cp ../PathogenicFungiObjectPrep/PathogenASVsFull.RDS ./
cp ../Fungi/2.PhyloseqObjectPrep/FUNGIFINALWOOD.RDS ./
cp ../Bacteria/2.PhyloseqObjectPrep/BACTERIAFINALWOOD.RDS ./
Rscript network_vine_wood.R
cd ..
```


## Code Usage disclaimer<a name="disclaimer"></a>

The following is the disclaimer that applies to all scripts, functions, one-liners, etc. This disclaimer supersedes any disclaimer included in any script, function, one-liner, etc.

You running this script/function means you will not blame the author(s) if this breaks your stuff. This script/function is provided **AS IS** without warranty of any kind. Author(s) disclaim all implied warranties including, without limitation, any implied warranties of merchantability or of fitness for a particular purpose. The entire risk arising out of the use or performance of the sample scripts and documentation remains with you. In no event shall author(s) be held liable for any damages whatsoever (including, without limitation, damages for loss of business profits, business interruption, loss of business information, or other pecuniary loss) arising out of the use of or inability to use the script or documentation. Neither this script/function, nor any part of it other than those parts that are explicitly copied from others, may be republished without author(s) express written permission. Author(s) retain the right to alter this disclaimer at any time. This disclaimer was copied from a version of the disclaimer published by other authors in https://ucunleashed.com/code-disclaimer and may be amended as needed in the future.


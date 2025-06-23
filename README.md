# 3700 Assignment 2 : Dispersal from deep-sea vents inferred from mtDNA 

# These code blocks will help you download public dna sequences, align them, and create phylogenies to complete assignment 2. Please note that you should have the most up to date version of R installed along with R studio and that you are connected to the internet throughout. 

# To start off, if you need a new version of R, or if this is the first time you're using R, visit here https://cran.rstudio.com/ and choose the download that corresponds to your computer. AFTER installing R, you can visit here to download Rstudio desktop for your computer https://posit.co/download/rstudio-desktop/ .  Make sure to install R BEFORE installing R studio.  

# This first block of code clears your working environment in case you've been using R for something else and then confirms what version of R you are running

rm(list=ls())

# You need a new version of R - at least 4.3.3. You  check this by entering:

R.version

# The next block of code tells you where your working directory is.  The 'wd' is where to put your input files and where to check for any output files you produce

# If you want to change the 'wd' - use the optional setwd command below, and replace the ____ between the "_____" with where you'd like to save your material (and then remove the '#')

getwd()

# setwd("_____")

setwd("C:\\Users\\malex\\Sync\\R")

# This next block of code will install a package called BioManager that will then install 9 further packages if you have not already installed them. If you have already installed them, skip to the library() commands below and these will open the packages you need. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager");

BiocManager::install("bipartite");
BiocManager::install("ape");
BiocManager::install("DECIPHER");
BiocManager::install("phangorn");
BiocManager::install("ggtree");
BiocManager::install("ggplot2");
BiocManager::install("phytools");
BiocManager::install("picante")
BiocManager::install("imager")

library(ape); 
library(DECIPHER);
library(phangorn);
library(bipartite);
library(ggtree);
library(ggplot2);
library(phytools);
library(picante);
library(imager);


# This next block of code is going to plot a map of the sampling sites where your copepods were collected. 

map<-load.image("3700 genetic sampling sites.jpg")
plot(map,axes=FALSE, main = "Deep-sea vents where copepods were sampled")


# The next block of code goes to GenBank and downloads public sequences. These sequences are COI mitochondrial DNA barcodes from two species of deep-sea vent invertebrate. You can explore GenBank at https://www.ncbi.nlm.nih.gov/nucleotide/

dispersal_sequences = read.GenBank(c("OQ693582", "OQ693581", "OQ693580", "OQ693579", "OQ693578", "OQ693577", "OQ693576", "OQ693575", "OQ693574", "OQ693573", "OQ693572", "OQ693571", "OQ693570", "OQ693569", "OQ693568", "OQ693567", "OQ693566", "OQ693565", "OQ693564", "OQ693563", "OQ693497", "OQ693478", "OQ693473", "OQ693460", "OQ693458", "OQ693457", "OQ693434", "OQ693415", "OQ693414", "OQ693413", "OQ693104", "OQ693098", "OQ693097", "OQ693096", "OQ693087", "OQ693069", "OQ693054", "OQ693049", "OQ693044", "OQ693042"))


write.dna(dispersal_sequences, file = 'dispersal_sequences.fasta', format = 'fasta' )



# Next, you need to align these sequences. Aligning sequences ensures that you are comparing homologous regions with each other so that your phylognies will make appropriate branching relationships based on their similarity.

fas <- "dispersal_sequences.fasta"

dna <- readDNAStringSet(fas)

# The next command opens a web browser to visualise your un-aligned DNA sequences. See how they look like puzzle pieces you have just dumped on the table?  Aligning will resolve this puzzle.

BrowseSeqs(dna)

# The code that follows will complete an alignment with gap opening and extending costs set high. This is because your DNA is codes for a protein (ie it has a job) and so gaps would result in your DNA  not making amino acids.

DNA.no_trans.1 <- AlignSeqs(dna, gapOpening = c(-20, -10), gapExtension = c(-5, -1))

# You can visualise your alignment in a web browser 

BrowseSeqs(DNA.no_trans.1)

# Curious what your alignment did? Jump back and forth between the browser tabs of your aligned and unaligned DNA sequences. You can see homologous regions have been placed alongside each other with the use of "-" gaps to slide the sequences forwards or backwards to make sure like is being compared with like. 

# The next code block saves your aligned DNA sequences as a 'fasta' file - standard for dna sequences. If you're curious, you can use a notepad to open it and explore if you'd like. 

writeXStringSet(DNA.no_trans.1, file="dispersal_sequences_aligned.fasta")

dispersal_sequences_align <- read.dna("dispersal_sequences_aligned.fasta", format="fasta", as.matrix=TRUE)

write.dna(DNA.no_trans.1, file = 'dispersal_sequences_aligned.fasta', format = 'fasta' )

# Now that you've added sequences and aligned them it's time to make a phylogeny with your aligned sequences. There are many ways to make a phylogeny - we're going to try 2, ML and NJ. 

# This block will make a Maximum Likelihood (ML) tree called "rooted_ml_tree"

grass.phy <- read.phyDat('dispersal_sequences_aligned.fasta', format = 'fasta', type = 'DNA')

grass.align <- read.dna('dispersal_sequences_aligned.fasta', format = 'fasta')

grass.phy <- phyDat(grass.align)

grass.phy

dist <- dist.ml(grass.align)

nj.tree <- nj(dist)

fit <- pml(nj.tree, data = grass.phy)

fit

fitJC <- optim.pml(fit, rearrangement = "NNI")

plot(fitJC$tree, main = "JC, NNI rearrangement")

write.tree(fitJC$tree, file="grass_ml.tre")

phy_ml = read.tree(file = "grass_ml.tre")

rooted_ml_tree <- midpoint.root(phy_ml)

plot(rooted_ml_tree)

# This block makes an Neighbour-Joining tree called 'phy'

D <- dist.dna(grass.align, model="K80")

D

class(D)

length(D)

phy <- nj(D)

class(phy)

rooted_tree <- midpoint.root(phy)

plot(rooted_tree)

# In the next block we're going to plot your tree in a ggplot fashion rather than base R using the package ggtree. First your nj tree (tip label is set large (to 7) for what you need eventually in printing - you can reduce it if you'd like. 

njtree = ggplot(rooted_tree) + geom_tree() + theme_tree()+ geom_treescale()+geom_tiplab(size=7)

njtree

# then your ml tree

mltree = ggplot(rooted_ml_tree) + geom_tree() + theme_tree()+ geom_treescale()+geom_tiplab(size=7)

mltree


# Congratulations - you have a made a phylogeny!  Now you need to append information about these sequences.  To do this, you need the file 3700 test genbank metadata.csv. This .csv file includes the site information associated with each sample sequence. Remember to make sure that this .csv file is in whatever directory/folder you set as the working directory

genbank_seq_metadata <- read.csv(file = "3700 test genbank metadata.csv",head=TRUE, sep=",", row.names = 1)

# The gheatmap command below plots the phylogeny you've created against the site information you've just uploaded. First for the NJ tree.

gheatmap(njtree, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w NJ tree")+ theme(legend.position="none")

# Now, use the same command structure to append the information to the ML tree.

gheatmap(mltree, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w ML tree")+ theme(legend.position="none")


# Now that you've made the two phylogenies and appended the site information, you can use the pdf command below to make a pdf of your phylogenies and the map so you can print them for your video. 

pdf("dispersal sequences with metadata - 250623.pdf", width = 18, height = 12) # Open a new pdf file

gheatmap(njtree, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w NJ tree")+ theme(legend.position="none")

gheatmap(mltree, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w ML tree")+ theme(legend.position="none")

plot(map,axes=FALSE, main = "Deep-sea vents where copepods were sampled")

dev.off()

# So - hats off to you!! You've made two kinds of phylogeny from publicly available DNA sequences that were collected from two genera of deep-sea vent copepods. Now, print your pdf, examine the map from the assignment, and prepare to speak for three minutes (**without notes**) about the conclusions you might make regarding the larval dispersal of the two genera based on your phylogeny. Which taxon is likely to possess planktotrophic larvae?  Which taxon is likely to possess lecithotrophic larvae?  Why?  What consequences would mining in the Lau Basin have on species living at and around deep-sea vents? 




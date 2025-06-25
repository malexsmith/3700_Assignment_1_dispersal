*#* ZOO*3700 Assignment 1 : Dispersal among and between deep-sea vents inferred from public mtDNA sequences

*#* These code blocks will help you download public DNA sequences, align them, and create phylogenies to complete Assignment 1. Please note that you should have the most up to date version of R installed along with R studio and that you are connected to the internet throughout. 

# To start off, if this is the first time you're using R (or if you need a new version), visit here https://cran.rstudio.com/ and choose the download that corresponds to your computer. Then, AFTER installing R, visit here to download Rstudio desktop for your computer https://posit.co/download/rstudio-desktop/.  

# Make sure to install R BEFORE installing R studio.  

# If you already have R installed on your computer, check your version - you will need at least version 4.3.3 of R. You can check what version you have by entering:

R.version

# This assignment is designed to be followed along in steps from start to finish.  Each step creates an output that subsequent steps depend on.  Remember to follow the code blocks in order!

# This first block of code will clear your working environment in case you've been using R for something else in the past.

rm(list=ls())

# The next block of code will tell you where your working directory is (i.e. any figures you generate you can find in that folder).  The 'wd' is also where to put your input files.

getwd()

# If you want to change the 'wd' - use the optional setwd command below, and replace the ____ between the "_____" with where you'd like to save your material (and then remove the '#')

# setwd("_____")

setwd("C:\\Users\\malex\\Sync\\R")

# This next block of code will install a package called BioManager that we will use to install 9 further packages (if you have not already installed them). If you have already installed them, you can skip to the library() commands below which will open the packages you need to complete this assignment. 

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

# The next block are the library() commands which will 9 the packages you need to complete this assignment. 

library(ape); 
library(DECIPHER);
library(phangorn);
library(bipartite);
library(ggtree);
library(ggplot2);
library(phytools);
library(picante);
library(imager);

# This next block of code is going to plot a map of the deep-sea vent sampling sites where your deep-sea vent taxa were collected. In this case, for two species of copepod collected at two basins in the Western Pacific. 

map<-load.image("3700 genetic sampling sites.jpg")
plot(map,axes=FALSE, main = "Deep-sea vents where copepod spp. were sampled")


# The next block of code goes to GenBank and downloads the public DNA sequences for your two species of deep-sea vent copepod. These sequences are COI mitochondrial DNA barcodes. You can explore GenBank at https://www.ncbi.nlm.nih.gov/nucleotide/

dispersal_sequences = read.GenBank(c("OQ693582", "OQ693581", "OQ693580", "OQ693579", "OQ693578", "OQ693577", "OQ693576", "OQ693575", "OQ693574", "OQ693573", "OQ693572", "OQ693571", "OQ693570", "OQ693569", "OQ693568", "OQ693567", "OQ693566", "OQ693565", "OQ693564", "OQ693563", "OQ693497", "OQ693478", "OQ693473", "OQ693460", "OQ693458", "OQ693457", "OQ693434", "OQ693415", "OQ693414", "OQ693413", "OQ693104", "OQ693098", "OQ693097", "OQ693096", "OQ693087", "OQ693069", "OQ693054", "OQ693049", "OQ693044", "OQ693042"))



write.dna(dispersal_sequences, file = 'dispersal_sequences.fasta', format = 'fasta' )


# Next, you need to align these sequences. Aligning sequences ensures that you are comparing homologous regions with each other so that your phylognies will make appropriate branching relationships based on their similarity.

fas <- "dispersal_sequences.fasta"

dna <- readDNAStringSet(fas)

# The next command opens a web browser to visualise your un-aligned DNA sequences. See how they look like puzzle pieces you have just dumped on the table?  Aligning will resolve this puzzle.

BrowseSeqs(dna)

# The code that follows will complete an alignment with gap opening and extending costs set high. This is because your DNA is codes for a protein (ie it has a job) and so gaps would result in your DNA  not making amino acids.

DNA.no_trans.1 <- AlignSeqs(dna, gapOpening = c(-20, -10), gapExtension = c(-5, -1))

# Now, you can visualise your alignment as before in a new browser tab. Can you spot your two species?

BrowseSeqs(DNA.no_trans.1)

# Curious what your alignment did? Jump back and forth between the browser tabs of your aligned and unaligned DNA sequences. You can see homologous regions have been placed alongside each other with the use of "-" gaps to slide the sequences forwards or backwards to make sure like is being compared with like. 

# The next code block saves your aligned DNA sequences as a 'fasta' file - standard for dna sequences. If you're curious, use a notepad program to open the fasta file and explore. 

writeXStringSet(DNA.no_trans.1, file="dispersal_sequences_aligned.fasta")

dispersal_sequences_align <- read.dna("dispersal_sequences_aligned.fasta", format="fasta", as.matrix=TRUE)

write.dna(DNA.no_trans.1, file = 'dispersal_sequences_aligned.fasta', format = 'fasta' )

# Now that you've added sequences and aligned them it's time to make a phylogeny with your aligned sequences. There are many ways to make a phylogeny - we're going to try 2, ML and NJ. 

# This block will make a Maximum Likelihood (ML) tree called "phy_ml"

dispersal.align <- read.dna('dispersal_sequences_aligned.fasta', format = 'fasta')

dispersal.phy <- phyDat(dispersal.align)

dispersal.phy

dist <- dist.ml(dispersal.align)

nj.tree <- nj(dist)

fit <- pml(nj.tree, data = dispersal.phy)

fit

fitJC <- optim.pml(fit, rearrangement = "NNI")

plot(fitJC$tree, main = "JC, NNI rearrangement")

write.tree(fitJC$tree, file="dispersal_ml.tre")

phy_ml = read.tree(file = "dispersal_ml.tre")

# Your sequences are from two species - this next block will place the root of your phylogeny at the midpoint for "rooted_ml_tree". 

rooted_ml_tree <- midpoint.root(phy_ml)

# Plot the newly rooted phylogeny - can you see the two species? (hint - they are monophyletic) 

plot(rooted_ml_tree)

# This block makes an Neighbour-Joining tree called 'rooted_tree' based on the pairwise distances between each DNA sequence (in this case, it's Kimura's 2 Parameter distance). 

D <- dist.dna(dispersal.align, model="K80")

D

class(D)

phy <- nj(D)

class(phy)

# As before, since your sequences are from two species - this code block will root  your phylogeny at the midpoint. 

rooted_tree <- midpoint.root(phy)

plot(rooted_tree)

# In the next block we're going to plot your tree in a more visually pleasing ggplot fashion rather than base R using the package ggtree. First your nj tree (Note that the tip label is set large here (to 7) for what you need eventually in printing - you can reduce it if you'd like. 

njtree = ggplot(rooted_tree) + geom_tree(linewidth = 1) + theme_tree()+ geom_treescale(linesize = 1, fontsize =7)+geom_tiplab(size=7)

njtree

# then your ml tree

mltree = ggplot(rooted_ml_tree) + geom_tree(linewidth = 1) + theme_tree()+ geom_treescale(linesize = 1, fontsize =7)+geom_tiplab(size=7)

mltree


# Congratulations - you have a made phylogenies from publicly available deep-sea vent species!  Now you need to append the metadata about these sequences (what basin were they from). To do this, you need the file 3700 test genbank metadata.csv. This .csv file includes the site information associated with each sample sequence. Remember to make sure that this .csv file is in whatever directory/folder you set as the working directory

genbank_seq_metadata <- read.csv(file = "3700 test genbank metadata.csv",head=TRUE, sep=",", row.names = 1)

# The gheatmap command below plots the phylogeny you've created against the site information you've just uploaded. First for the NJ tree.

gheatmap(njtree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w NJ tree")+ theme(legend.position="none")

# Now, use the same command structure to append the information to the ML tree.

gheatmap(mltree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w ML tree")+ theme(legend.position="none")


# Now that you've made the two phylogenies and appended the site information, you can use the pdf command below to make a single Acrobat file of your phylogenies and the map which you should print and have on hand for your video submission of this assignment. 

pdf("ZOO3700 deep-sea vent sequences with metadata - dispersal assignment - 2506242.pdf", width = 18, height = 12) # Open a new pdf file

plot(map,axes=FALSE, main = "Deep-sea vents where copepods were sampled")

gheatmap(njtree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.05, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w NJ tree")+ theme(legend.position="none")

gheatmap(mltree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.05, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w ML tree")+ theme(legend.position="none")


dev.off()

# So - hats off to you!! You've made two kinds of phylogeny from publicly available DNA sequences that were collected from two species of deep-sea vent copepods. Now, print your pdf, examine the map from the assignment, and prepare to speak for three minutes (!!without notes!!) about the conclusions you might make regarding the larval dispersal of the two genera based on your phylogeny. 
#  
# Which taxon is likely to possess planktotrophic larvae?  Why? 
# Which taxon is likely to possess lecithotrophic larvae?  Why? 
#  
# Does the way you made your phylogeny change your prediction? 
# What consequences would mining in the Lau Basin have on species living at and around deep-sea vents? 

# good luck!


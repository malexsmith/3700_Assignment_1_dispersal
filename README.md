# 3700-dispersal
test

# code blocks to download public dna sequences, align them, and create phylogenies to complete assignment 2
# please note that you should have the most up to date version of R installed
# along with R studio
# and that you are connected to the internet throughout.  

# this first block of code clears your working environment in case you've been using 
# R for something else.  
# and then confirms what version of R you are running

rm(list=ls())

## You need a new version of R - at least 4.3.3. 
# check this by entering:
R.version

# next block of code tells you where your working directory is 
# this is where to put your input files and where to check for any output files you produce
# if you want to change it - use the optional setwd command beneath it 
## and replace the ____ between the "_____" with where you'd like to save your material

getwd()
# alter alter what you want for wd below. 
# setwd("_____")
setwd("C:\\Users\\malex\\Sync\\R")

# the next block of code will install a package called BioManager
# that will then install 8 further packages if you have not already installed them

# if you have already installed them, skip to the library() commands below and 
# and these will open the packages you need. 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bipartite")
BiocManager::install("ape")
BiocManager::install("DECIPHER")
BiocManager::install("phangorn")
BiocManager::install("ggtree")
BiocManager::install("ggplot2")
BiocManager::install("phytools")
BiocManager::install("picante")

library(ape) 
library(DECIPHER)
library(phangorn)
library(bipartite)
library(ggtree)
library(ggplot2)
library(phytools)
library(picante)

# The next block of code goes to GenBank and downloads public sequences 
# these sequences are COI mitochondrial DNA barcodes from two species of deep-sea vent invertebrate
# You can explore GenBank at https://www.ncbi.nlm.nih.gov/nucleotide/


# https://jhudatascience.org/AnVIL_Phylogenetic-Techniques/downloading-the-sequences-from-genbank.html

dispersal_sequences = read.GenBank(c("OQ693582", "OQ693581", "OQ693580", "OQ693579", "OQ693578", "OQ693577", "OQ693576", "OQ693575", "OQ693574", "OQ693573", "OQ693572", "OQ693571", "OQ693570", "OQ693569", "OQ693568", "OQ693567", "OQ693566", "OQ693565", "OQ693564", "OQ693563", "OQ693497", "OQ693478", "OQ693473", "OQ693460", "OQ693458", "OQ693457", "OQ693434", "OQ693415", "OQ693414", "OQ693413", "OQ693104", "OQ693098", "OQ693097", "OQ693096", "OQ693087", "OQ693069", "OQ693054", "OQ693049", "OQ693044", "OQ693042"))
dispersal_sequences3 = read.GenBank(c("OQ693582", "OQ693581", "OQ693580", "OQ693579", "OQ693578", "OQ693577", "OQ693576", "OQ693575", "OQ693574", "OQ693573", "OQ693572", "OQ693571", "OQ693570", "OQ693569", "OQ693568", "OQ693567", "OQ693566", "OQ693565", "OQ693564", "OQ693563", "OQ693497", "OQ693478", "OQ693473", "OQ693460", "OQ693458", "OQ693457", "OQ693434", "OQ693415", "OQ693414", "OQ693413", "OQ693104", "OQ693098", "OQ693097", "OQ693096", "OQ693087", "OQ693069", "OQ693054", "OQ693049", "OQ693044", "OQ693042"), species.names = TRUE)


write.dna(dispersal_sequences, file = 'dispersal_sequences.fasta', format = 'fasta' )


write.dna(dispersal_sequences3, file = 'dispersal_sequences3.fasta', format = 'fasta' )



# Next, you need to align these sequences
# aligning sequences ensures that you are comparing homologous regions with each other
# so that your phylognies will make appropriate branching relationships based on how similar they are


fas <- "dispersal_sequences3.fasta"

dna <- readDNAStringSet(fas)
# The next command opens a web browser to visualise your un-aligned DNA sequences 
## see how they look like puzzle pieces you have just dumped on the table?  Aligning will resolve this puzzle.
BrowseSeqs(dna)

##  align with gap opening and extending costs set high
# because your DNA is codes for a protein (ie it has a job) and so gaps implies that the DNA is not making amino acids

DNA.no_trans.1 <- AlignSeqs(dna, gapOpening = c(-20, -10), gapExtension = c(-5, -1))

## now you can visualise your alignment in a web browser 
BrowseSeqs(DNA.no_trans.1)

# jump back and forth between your aligned and unaligned DNA sequences
## you can see homologous regions have been placed alongside each other with the use of "-" gaps to slide 
# the sequences forwards or backwards to make sure like is being compared with like. 

# next code block saves your aligned DNA sequences as a 'fasta' file - standard for dna sequences.  
# you can use a notepad to open it and explore if you'd like. 

writeXStringSet(DNA.no_trans.1, file="dispersal_sequences_aligned.fasta")
dispersal_sequences_align <- read.dna("dispersal_sequences_aligned.fasta", format="fasta", as.matrix=TRUE)
write.dna( DNA.no_trans.1, file = 'dispersal_sequences_aligned.fasta', format = 'fasta' )


# Now that you've added sequences and aligned them
## it's time to make a phylogeny with your aligned sequences
## There are many ways to make a phylogeny - we're going to try 2, ML and NJ. 


## this block will make a Maximum Likelihood (ML) tree called "rooted_ml_tree"
## clean this block
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


fitJC <- optim.pml(fit, rearrangement = "stochastic", control = pml.control(trace = 0))
plot(fitJC$tree, main = "JC, stochastic rearrangement")

write.tree(fitJC$tree, file="grass_ml.tre")

phy_ml = read.tree(file = "grass_ml.tre")

rooted_ml_tree <- midpoint.root(phy_ml)
plot(rooted_ml_tree)

# http://adegenet.r-forge.r-project.org/files/PRstats/practical-introphylo.1.0.pdf


# this block makes an Neighbour-Joining tree called 'phy'
D <- dist.dna(grass.align, model="K80")
D
class(D)
## [1] "dist"
length(D)

phy <- nj(D)
class(phy)
plot(phy)

rooted_tree <- midpoint.root(phy)

## In the next block we're going to plot your tree in a ggplot fashion rather than base R using ggtree

# first your nj tree (tip label is set large (to 7) for what you need eventually in printing - you can reduce it if you'd like. 
test = ggplot(rooted_tree) + geom_tree() + theme_tree()+ geom_treescale()+geom_tiplab(size=7)
test

# then your ml tree
test2 = ggplot(rooted_ml_tree) + geom_tree() + theme_tree()+ geom_treescale()+geom_tiplab(size=7)
test2

pdf("dispersal sequences test - 25066.pdf", width = 18, height = 12) # Open a new pdf file
test
dev.off()

## this .csv file includes the site information associated with each sample

genbank_seq_metadata <- read.csv(file = "3700 test genbank metadata.csv",head=TRUE, sep=",", row.names = 1)

## - the gheatmap command below plots the phylogeny you've create
## against the site information you've just uploaded.
# first for the NJ tree

gheatmap(test, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal")+ theme(legend.position="none")

## then for the ML tree

gheatmap(test2, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal")+ theme(legend.position="none")


## the pdf command below creates a pdf of the plots above. 

pdf("dispersal sequences with metadata - 250622.pdf", width = 18, height = 12) # Open a new pdf file
gheatmap(test, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal")+ theme(legend.position="none")
gheatmap(test2, genbank_seq_metadata , low = "white",high = "#1099dd",offset=0.03, width=0.15, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal")+ theme(legend.position="none")
dev.off()

## so - hats off to you - you've made two kinds of  phylogeny from publicly available DNA sequences 
## that were collected from 

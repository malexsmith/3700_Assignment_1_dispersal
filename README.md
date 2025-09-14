# ZOO*3700 Assignment 1 : Dispersal among and between deep-sea vents inferred from public mtDNA sequences

# Welcome to your first invertebrateR assignment! These code blocks will help you download public DNA sequences, align them, and create phylogenies to complete Assignment 1. Please note that you should have the most up to date version of R installed along with R studio and that you are connected to the internet throughout. 

# This assignment has been written in R using Windows and I recommend that you complete the assisgnment on a Windows machine.  

# Even if you have used R and RStudio in the past, I suggest that you uninstall that verison and reinstall following the instructions here.  This will help you avoid MANY unnecesssary headaches!

# To start off, you need to download R.  So visit here https://cran.rstudio.com/ and choose the download that corresponds to your computer. Make sure to install R BEFORE installing R studio.

# Then, AFTER installing R, visit here https://posit.co/download/rstudio-desktop/ to download Rstudio desktop for your computer.   

# Make sure to install R BEFORE installing R studio.  

# Once you have both R and R studio installed, download the three files (.R, .csv and .jpg) from the GitHub repository and put them in your working directory (don't know here that is?  we'll get to that in a second!)

# Open up the .R extension file from R.Studio.  You can either work from the console window the (">") and follow along here by copying the commands from GitHub (the smaller text are the commands to copy, the larger text is me talking to you) to the R console window, or by selecting the command blocks in the .R window and then pressing "Run". This option will be faster. 

# As you work through this assignment do not run the entire code at once.  This assignment is designed to be followed along in steps (code blocks) from start to finish.  Some are very simple, some are more complicated, but generally, each step will create an output that subsequent steps will depend on. So follow the code blocks in order and run them one at a time!

# Here's the first command block - let's confirm what version of R you downloaded and installed (If you are persevering with the version of R you have already installed, confirm that you have at least version 4.3.3.

R.version

# The second block of code will clear your working environment in case you've been using R for something else in the past.

rm(list=ls())

# The next block of code will tell you where your working directory is (i.e. any figures you generate you can find in that folder).  The 'wd' is also where to put your input files.

getwd()

# This next block of code will install a package called BioManager that we will need to install 9 further packages 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager");

# The next block uses BioManager package to install the 8 separate packages you need to run the assignment (if you have not already installed them). If you have already installed them, you can skip to the library() commands below which will open the packages you need to complete this assignment. 

# Note - that if you have previously installed any of these R packages, you can likely skip ahead to the library() commands below.  If you proceed with the install commands you might be asked whether to update all, some or none ('a' 's' 'n') each time to move forward with the code.

BiocManager::install("ape");
BiocManager::install("DECIPHER");
BiocManager::install("phangorn");
BiocManager::install("ggtree");
# BiocManager::install("ggplot2");
BiocManager::install("phytools");
BiocManager::install("picante");
BiocManager::install("imager");

# Note 240914 - the most recent ggplot2 has a bug in it.  Please make sure you are using version 3.5.2 The line below will install that for you. 

packageurl <- "https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.5.2.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

# The next block are the library() commands which will open the 8 packages needed to work through this assignment. 

library(ape); 
library(DECIPHER);
library(phangorn);
library(ggtree);
library(ggplot2);
library(phytools);
library(picante);
library(imager);

# This next block of code is going to plot a map of the deep-sea vent sampling sites where your deep-sea vent taxa were collected. In this case, for two species of copepod collected at two basins in the Western Pacific. 

map<-load.image("3700 genetic sampling sites.jpg");

plot(map,axes=FALSE, main = "Deep-sea vents where copepod spp. were sampled")


# The next block of code goes to GenBank and downloads the public DNA sequences for your two species of deep-sea vent copepod. These sequences are COI mitochondrial DNA barcodes. GenBank is a global repository for all DNA sequences and if that's interesting to you, you can explore this fantastic resource further at https://www.ncbi.nlm.nih.gov/nucleotide/

dispersal_sequences = read.GenBank(c("OQ693582", "OQ693581", "OQ693580", "OQ693579", "OQ693578", "OQ693577", "OQ693576", "OQ693575", "OQ693574", "OQ693573", "OQ693572", "OQ693571", "OQ693570", "OQ693569", "OQ693568", "OQ693567", "OQ693566", "OQ693565", "OQ693564", "OQ693563", "OQ693497", "OQ693478", "OQ693473", "OQ693460", "OQ693458", "OQ693457", "OQ693434", "OQ693415", "OQ693414", "OQ693413", "OQ693104", "OQ693098", "OQ693097", "OQ693096", "OQ693087", "OQ693069", "OQ693054", "OQ693049", "OQ693044", "OQ693042"))



write.dna(dispersal_sequences, file = 'dispersal_sequences.fasta', format = 'fasta' )


# In the next code block, you will align these sequences. Aligning DNA sequences ensures that you are comparing homologous regions with each other so that your phylognies will make appropriate branching relationships based on their similarity.

fas <- "dispersal_sequences.fasta"

dna <- readDNAStringSet(fas)

# The next code block will open a web browser to visualise your un-aligned DNA sequences. See how they look like puzzle pieces you have just dumped on the table?  Aligning will resolve this puzzle.

BrowseSeqs(dna)

# The next code block will align your DNA sequences with gap opening and extending costs set high. This is because your DNA is codes for a protein (ie it has a job) and so and actual gaps between nucleotides would result in your DNA not making amino acids.

DNA.no_trans.1 <- AlignSeqs(dna, gapOpening = c(-20, -10), gapExtension = c(-5, -1))

# Now, you can visualise your alignment as before in a new browser tab. Can you spot your two species in this visualisation?

BrowseSeqs(DNA.no_trans.1)

# Curious what your alignment did? Jump back and forth between the browser tabs of your aligned and unaligned DNA sequences. You can see homologous regions have been placed alongside each other with the use of "-" gaps to slide the sequences forwards or backwards to make sure like is being compared with like. 

# The next code block saves your aligned DNA sequences as a 'fasta' file - standard format for dna sequences. If you're curious, use a free program called BioEdit to open the fasta file and explore. https://thalljiscience.github.io/. It's an older program but does a lot and has been a workhorse for years for those that analyse DNA sequences. 

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

# Your sequences are from two species and this next block will place the root of your phylogeny at the midpoint for "rooted_ml_tree". 

rooted_ml_tree <- midpoint.root(phy_ml)

# The next code block will plot the newly rooted phylogeny - can you see the two species? (hint - they are monophyletic) 

plot(rooted_ml_tree)

# This block makes an Neighbour-Joining tree called 'rooted_tree' based on the pairwise distances between each DNA sequence (in this case, it's Kimura's 2 Parameter distance, or K2P). 

D <- dist.dna(dispersal.align, model="K80")

D

class(D)

phy <- nj(D)

class(phy)

# As before, since your sequences are from two species - this code block will root your phylogeny at the midpoint. 

rooted_tree <- midpoint.root(phy)

plot(rooted_tree)

# In the next block we're going to plot your tree in a more visually pleasing manner rather than base R using the packages ggtree and ggplot. First your NJ tree (Note that the tip label is set large here (to 7) for what you need eventually in printing to pdf - you can reduce it if you'd like using the fontsize parameter below - try 3 - but return it to seven before you print). 

njtree = ggplot(rooted_tree) + geom_tree(linewidth = 1) + theme_tree()+ geom_treescale(linesize = 1, fontsize =7)+geom_tiplab(size=7)

njtree

# The next code block is do the same for your ML tree

mltree = ggplot(rooted_ml_tree) + geom_tree(linewidth = 1) + theme_tree()+ geom_treescale(linesize = 1, fontsize =7)+geom_tiplab(size=7)

mltree


# Congratulations - you have a made phylogenies from publicly available deep-sea vent species!  Now you need to append the metadata about these sequences (what basin on your map were they from). 

# To do this, you need the file 3700 test genbank metadata.csv. This .csv file includes the site information associated with each sample sequence. Remember to make sure that this .csv file is in whatever directory/folder you set as the working directory.  The next code block uploads that .csv into your R environment. 

genbank_seq_metadata <- read.csv(file = "3700 test genbank metadata.csv",head=TRUE, sep=",", row.names = 1)

# The next code block uses the gheatmap command below to plot the phylogeny you've created against the site information you've just uploaded (this is possible because the tip labels of the phylogeny and the header of the column of the .csv file are the same values). First for the NJ tree. (Make sure to add you name to the title! Find ggtitle and insert your name inbetween the quotes). 

gheatmap(njtree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w NJ tree")+ theme(legend.position="none")

# Don't worry if the tree and matrix overlap in the plot window - the sizing has been made for the final pdf output.  

# Now, use the same command structure to append the information to the ML tree. (Again, don't worry about overlap and make sure to add you name to the title! Find ggtitle and insert your name inbetween the quotes). 

gheatmap(mltree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.15, font.size=3, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w ML tree")+ theme(legend.position="none")

# Now that you've made the two phylogenies and appended the site information, the next code block uses the pdf command below to make a single Acrobat file of your map and phylogenies. A printout of this pdf is what you should have on hand as a visual aid for your video submission of this assignment. 

pdf("ZOO3700 deep-sea vent sequences with metadata - dispersal assignment - 2506242.pdf", width = 18, height = 12) # Open a new pdf file

plot(map,axes=FALSE, main = "Deep-sea vents where copepods were sampled")

gheatmap(njtree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.05, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w NJ tree")+ theme(legend.position="none")

gheatmap(mltree, genbank_seq_metadata , low = "white",high = "#1099dd",color="grey", offset=0.03, width=0.05, font.size=6, 
         colnames_angle=90, hjust=1)+vexpand(.1, -1)+ ggtitle("Deep Sea Vent Dispersal w ML tree")+ theme(legend.position="none")

dev.off()

# So - hats off to you!! You've made two kinds of phylogeny from publicly available DNA sequences that were collected from multiple populations of two species of deep-sea vent copepods. 

# Now, print your pdf (hard copy), examine the map and phylogenies you created in this assignment. 

# The final part of your assignment is to record yourself using the print out as a visual aid as you speak for three minutes (!!without notes!!) about the conclusions you made regarding the larval dispersal of the two genera based on your phylogeny. 

# Which taxon is likely to possess planktotrophic larvae? Which taxon is likely to possess lecithotrophic larvae? Why?

# Does the way you made your phylogeny change your prediction? Why or why not? 

# Based on your phylogeny, which of the two basins contains populations that would you estimate to be more fragmented/isolated? Why? Is a metapopulation approach useful here?  

# If mining companies were to target the region for development, would this have the same effect on all species living at and around deep-sea vents in each basin? Why or why not? 

# good luck!


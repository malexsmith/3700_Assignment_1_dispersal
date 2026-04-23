# ZOO*3700 Assignment 1 : Dispersal among and between deep-sea hydrothermal vents inferred from public mtDNA sequences

# Welcome to your first Invert-R assignment! These code blocks will help you map collection localities, download, and align, public DNA sequences, and create phylogenies to complete Assignment 1. Please note that you should have the most up to date version of R installed along with R studio and that you are connected to the internet throughout. 

# This assignment has been written in R using Windows and I recommend that you complete the assignment on a Windows machine.  

# Even if you have used R and RStudio in the past, I suggest that you uninstall that version and re-install following the instructions here.  This will help you avoid MANY unnecessary headaches!

# To start off, you need to download R.  So visit here https://cran.rstudio.com/ and choose the download that corresponds to your computer. Make sure to install R BEFORE installing R studio.

# Then, AFTER installing R, visit here https://posit.co/download/rstudio-desktop/ to download Rstudio desktop for your computer.   

# Make sure to install R BEFORE installing R studio.  

# Once you have both R and R studio installed, download the three files (.R and 2 .csv files) from the GitHub repository and put them in your working directory (don't know here that is?  we'll get to that in a second!)

# Open up the .R extension file from R.Studio.  You can either work from the console window the (">") and follow along here by copying the commands from GitHub (the smaller text are the commands to copy, the larger text is me talking to you) to the R console window, or by selecting the command blocks in the .R window and then pressing "Run". This option will be faster. 

# As you work through this assignment do not run the entire code at once.  This assignment is designed to be followed along in steps (code blocks) from start to finish.  Some are very simple, some are more complicated, but generally, each step will create an output that subsequent steps will depend on. So follow the code blocks in order and run them one at a time!

# Here's the first command block - let's confirm what version of R you downloaded and installed (If you are persevering with the version of R you have already installed, confirm that you have at least version 4.3.3.

R.version

# The second block of code will clear your working environment in case you've been using R for something else in the past.

rm(list=ls())

# The next block of code will tell you where your working directory is (i.e. any figures you generate you can find in that folder).  The 'wd' is also where to put your input files.

getwd()
setwd("C://Users//malex//Documents//3700_Assignment_1_dispersal")

# This next block of code will install a package called BioManager that we will need to install 9 further packages 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager");

# The next block uses BioManager package to install the 7 separate packages you need to run the assignment (if you have not already installed them). If you have already installed them, you can skip to the library() commands below which will open the packages you need to complete this assignment. 

# Note - that if you have previously installed any of these R packages, you can likely skip ahead to the library() commands below.  If you proceed with the install commands you might be asked whether to update all, some or none ('a' 's' 'n') each time to move forward with the code.

BiocManager::install("ape");
BiocManager::install("Biostrings");
BiocManager::install("DECIPHER");
BiocManager::install("phangorn");

BiocManager::install("ggtree");
BiocManager::install("phytools");
BiocManager::install("ggplot2");
BiocManager::install("sf");
BiocManager::install("rnaturalearth");
BiocManager::install("ggspatial");
BiocManager::install("vegan");
BiocManager::install("cowplot");

# The next block are the library() commands which will open the 8 packages needed to work through this assignment. 

library(ape); 
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(ggtree);
library(ggplot2);
library(phytools);
library(sf)
library(rnaturalearth)
library(ggspatial)
library(vegan)
library(cowplot)

# This next block of code is going to plot a map of the deep-sea vent sampling sites where your deep-sea vent taxa were collected. In this case, for two species of copepod collected at two basins in the Western Pacific. 

map_sites = read.csv("3700 invertR1 data.csv")

# Convert lat/long to sf format
sites_sf <- st_as_sf(map_sites, coords = c("lon", "lat"), crs = 4326)

# Load Map Data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the Map
vent_map = ggplot(data = world) +
  geom_sf(fill = "antiquewhite", color = "gray50") +
  # Plot points
  geom_sf(data = sites_sf, aes(color = basin), size = 3) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  coord_sf(xlim = c(100, 190), ylim = c(-40, 10), expand = FALSE) +
  theme_minimal() +
  labs(title = "Deep-sea hydrothermal vents where copepods were sampled",
       x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "aliceblue"))
vent_map


# The next block of code goes to GenBank and downloads the public DNA sequences for your two species of deep-sea vent copepod. These sequences are COI mitochondrial DNA barcodes. GenBank is a global repository for all DNA sequences and if that's interesting to you, you can explore this fantastic resource further at https://www.ncbi.nlm.nih.gov/nucleotide/

dispersal_sequences1 = map_sites$accession

dispersal_sequences = read.GenBank(dispersal_sequences1)

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

# In the next block we're going to plot your tree in a more visually pleasing manner rather than base R using the packages ggtree and ggplot. First your NJ tree (Note that the tip label is set large here (to 7) for what you need eventually in printing to pdf - you can reduce it if you'd like using the fontsize parameter below - try 3 - but return it to seven before you print). 

mltree = ggplot(rooted_ml_tree) + geom_tree(linewidth = 1) + theme_tree()+ geom_treescale(fontsize =7)+geom_tiplab(size=7)

mltree

# Congratulations - you have constructed a ML phylogeny from publicly available deep-sea vent species!  
# Now you need to append the metadata about these sequences (what basin on your map were they from). 

# To do this, you need the file 3700 test genbank metadata.csv. This .csv file includes the site information associated with each sample sequence. Remember to make sure that this .csv file is in whatever directory/folder you set as the working directory.  The next code block uploads that .csv into your R environment. 

matrix = with(map_sites, table(map_sites$accession, map_sites$basin))


# Maintaining this structure while transforming to data.frame (using nested unclass)
matrix <- as.data.frame.matrix(unclass(matrix))
matrix

p <- ggtree(rooted_ml_tree)

# EXPLANATION 

colored_tree = p %<+% map_sites + 
  geom_tiplab(aes(color = sp), size = 2) + 
  scale_color_viridis_d(option = "viridis")+
  geom_treescale(x=0, y=25, fontsize=4, linesize=1)+# Color by column
  theme(legend.position = "null")
colored_tree

p_colored <- p %<+% map_sites +
  geom_tippoint(aes(color = sp), size = 2) +
  geom_treescale(x=0, y=25, fontsize=4, linesize=1)+# Color by column
  scale_color_viridis_d(option = "viridis") + # Use viridis
  theme(legend.position = "none")
p_colored

# Don't worry if the tree and matrix overlap in the plot window - the sizing has been made for the final pdf output.  

## 
colored_tree_heatmap = gheatmap(p_colored, matrix , low = "white",high = "#1099dd",color="grey",offset=0.0015, width=0.1, font.size=2.5, 
                                colnames_angle=45, hjust=1)+vexpand(.1, -1)+ theme(legend.position="none")
colored_tree_heatmap


#Now it's time to consider how genetic variation is related to geographic distance
# Isolation by distance (or IBD> 

# these are intra-specific calculations - so first you must break your dataset into two


# Let's first calculate species A genetic distances

target_id_a1 <- map_sites[map_sites$sp == "A", ]
target_id_a = target_id_a1$accession
subset_dna_a <- dispersal.align[target_id_a, ]
dist_a <- dist.ml(subset_dna_a)

# Now we will use this species A subset to calculate the pairwise geographic distance between samples
# Use the package sf to create a mappable simple feature object
sites_a <- st_as_sf(target_id_a1, coords = c("lon", "lat"), crs = 4326)
# Now calculate the pairwise distance matrix (result in meters by default for WGS84)
distance_matrix_a <- st_distance(sites_a, sites_a)
# You can set units explicitly if needed, e.g., for kilometers:
distance_matrix_km_a <- units::set_units(st_distance(sites_a, sites_a), "km")
distance_matrix_km_a <- as.dist(distance_matrix_km_a)

# Because these distances are pairwise - the measures are not independent of each other
# therefore rather than a regression, perform a Mantel test in vegan to test the relationship between genetic and geographic distance

mantel_res_a <- mantel(dist_a, distance_matrix_km_a, method = "spearman", permutations = 999)
print(mantel_res_a)

# Combine distances into a single data frame so you can plot them
plot_data_a <- data.frame(
  Geography = as.vector(distance_matrix_km_a),
  Diversity = as.vector(dist_a)
)

options(scipen = 999) # Set a high penalty to avoid scientific notation

mantel_plot_a = ggplot(plot_data_a, aes(x = Geography, y = Diversity)) +
  geom_point() +
  # geom_point(position = position_jitter(width = 20)) +
  # above applies a 20km jitter to points 
  geom_smooth(method = "glm") +
  theme_minimal() +
  labs(title = "Distance Decay in Genetic Diversity species A",
       x = "Geographic Distance (km)",
       y = "Genetic Diversity")
mantel_plot_a

# 3. species A genetic distances
# Names in FASTA object are usually in names(sequences)
target_id_b1 <- map_sites[map_sites$sp == "B", ]
target_id_b = target_id_b1$accession
subset_dna_b <- dispersal.align[target_id_b, ]
dist_b <- dist.ml(subset_dna_b)

sites_b <- st_as_sf(target_id_b1, coords = c("lon", "lat"), crs = 4326)
# Calculate the pairwise distance matrix (result in meters by default for WGS84)
distance_matrix_b <- st_distance(sites_b, sites_b)

# You can set units explicitly if needed, e.g., for kilometers:
distance_matrix_km_b <- units::set_units(st_distance(sites_b, sites_b), "km")
distance_matrix_km_b <- as.dist(distance_matrix_km_b)

# mantel test in vegan
mantel_res_b <- mantel(dist_b, distance_matrix_km_b, method = "spearman", permutations = 999)
print(mantel_res_b)

# Combine distances into a data frame
plot_data_b <- data.frame(
  Geography = as.vector(distance_matrix_km_b),
  Diversity = as.vector(dist_b)
)

# Plot
options(scipen = 999) # Set a high penalty to avoid scientific notation

mantel_plot_b = ggplot(plot_data_b, aes(x = Geography, y = Diversity)) +
  geom_point() +
  # geom_point(position = position_jitter(width = 20)) +
  # above applies a 20km jitter to points 
  geom_smooth(method = "glm") +
  theme_minimal() +
  labs(title = "Distance Decay in Genetic Diversity species B",
       x = "Geographic Distance (km)",
       y = "Genetic Diversity")
mantel_plot_b


# plot together
library(dplyr)
long_df <- bind_rows(list("Species B" = plot_data_b, "Species A" = plot_data_a), .id = "source")

mantel_plot_both = ggplot(long_df, aes(x = Geography, y = Diversity, color = source)) +
  geom_point() +
  # geom_point(position = position_jitter(width = 20)) +
  # above applies a 20km jitter to points 
  geom_smooth(method = "glm") +
  scale_color_viridis_d(option = "viridis")+
  theme_minimal() +
  labs(title = "Distance Decay with Genetic Diversity for two species of hydrothermal copepods",
       x = "Geographic Distance (km)",
       y = "Genetic Distance")
mantel_plot_both


## explanation
mantel_plot_both_no_colour = ggplot(long_df, aes(x = Geography, y = Diversity, color = source)) +
  geom_point() +
  # geom_point(position = position_jitter(width = 20)) +
  # above applies a 20km jitter to points 
  geom_smooth(method = "glm") +
 scale_color_viridis_d(option = "turbo")+
  theme_minimal() +
  labs(title = "Distance Decay with Genetic Diversity for two species of hydrothermal copepods",
       x = "Geographic Distance (km)",
       y = "Genetic Distance")+
  theme(legend.position = "none") 
mantel_plot_both_no_colour


# Now that you've made the two phylogenies and appended the site information, the next code block uses the pdf command below to make a single Acrobat file of your map and phylogenies. A printout of this pdf is what you should have on hand as a visual aid for your video submission of this assignment. 

pdf("ZOO3700 deep-sea vent sequences with metadata - dispersal assignment - 260423.pdf", width = 10, height = 6) # Open a new pdf file

vent_map
colored_tree_heatmap
mantel_plot_both_no_colour

dev.off()

# So - hats off to you!! You've made a map and a phylogeny and plotted isolation by distance from publicly available DNA sequences that were collected from two species of deep-sea vent copepods. 

# Now, print your pdf (hard copy), examine the map, the phylogeny and the plot you created in this assignment. 

# The final part of your assignment is to record yourself using the print out as a visual aid as you speak for three minutes (!!without notes!!) about the conclusions you made regarding the larval dispersal of the two genera based on your phylogeny. 

# Which taxon is likely to possess planktotrophic larvae? Which taxon is likely to possess lecithotrophic larvae? Why?

# ################# Which species is which in the IBD and phylogeny?

# Based on your phylogeny, which of the basins would you estimate is more fragmented/isolated? 

# If mining companies were to target the region for development, would this have the same effect on all species living at and around deep-sea vents? Why or why not? 

# GOOD LUCK!


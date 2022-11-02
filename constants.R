library(treeio) # Many functions to read in different tree formats, includes raxml with bootstrap.
library(seqinr) # Deal with fasta files.
library(ggplot2)
library(dplyr) # Format dataframes.
library(lubridate) # Deal with date format.
library(sf) # Deal with maps in R.
library(ape) # Read tree files.
library(tidyverse)
library(gridExtra)
library(ggtree)
library(ggrepel)
library(lemon) # Format nice looking data frames and tables.
library(readODS)

# Load fasta file and metadata files.
fasta_path = "XXXX"
fasta_file <- read.fasta(fasta_path, as.string = TRUE, forceDNAtolower = FALSE, set.attributes = FALSE)
metadata_path = "YYYY"
metadata_file <- read.table(metadata_path,  sep = '\t', header = TRUE, stringsAsFactors = FALSE)
lab_path = "ZZZZ"
lab_file <- read.table(lab_path,  sep = '\t', header = TRUE, stringsAsFactors = FALSE, quote = "£")
# Load laboratory tags file.
labtags_file <- read_ods("path/to/labtagsfile.ods", strings_as_factors = FALSE)
# Load map of Ireland file. The .shp file needs to be on same directory as .shx, .dbf, .cpg and .prj files.
# If .shp map doesn't exist or is missing some files, open with openjump and save as shape file to generate all required files.
map_path = "path/to/mapfile.shp"
ROI_map <- st_read(map_path)
# Make a list of mutations to report.
target_mutations <- list()
target_mutations[[1]] = c("mutA", "mutB")
target_mutations[[2]] = c("mutC")
target_mutations[[3]] = c("mutA", "mutB", "mutC")
# Get mafft and raxml commands on variables.
# mafft documentation: https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
mafft_command = paste0("mafft --auto --reorder ", fasta_path, " > mafft_aligned.fasta")
# raxml cocumentation: https://isu-molphyl.github.io/EEOB563/computer_labs/lab4/raxml-ng.html
# raxml github installation : https://github.com/amkozlov/raxml-ng
raxml_command = "/usr/bin/raxmlHPC-PTHREADS-AVX /usr/share/man/man1/raxmlHPC-PCTHREADS-AVX -T 12 -f a -x 123 -p 123 -N 100 -m GTRCAT -k -O -s mafft_aligned.fasta -n raxml_tree -w `pwd`"



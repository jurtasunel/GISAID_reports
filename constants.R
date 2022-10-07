library(seqinr) # Deal with fasta files.
library(ggplot2)
library(dplyr) # Format dataframes.
library(lubridate) # Deal with date format.
library(sf) # Deal with maps in R.
library(ape) # Read tree files.
library(tidyverse)
library(dplyr)

fasta_path = "/home/gabriel/Desktop/Jose/Projects/GISAID_reports/Data/August/gisaid_hcov-19_2022_09_29_15.fasta"
fasta_file <- read.fasta(fasta_path, as.string = TRUE, forceDNAtolower = FALSE, set.attributes = FALSE)
metadata_path = "/home/gabriel/Desktop/Jose/Projects/GISAID_reports/Data/August/gisaid_hcov-19_2022_09_29_15.tsv"
metadata_file <- read.table(metadata_path,  sep = '\t', header = TRUE, stringsAsFactors = FALSE)
map_path = "/home/gabriel/Desktop/Jose/Maps/26 Counties/Census2011_Admin_Counties_generalised20m.shp"
# If .shp map doesn't exist or is missing some files, open with openjump and safe as shape file to generate all required files.
ROI_map <- st_read(map_path)



source("/home/gabriel/Desktop/Jose/Projects/GISAID_reports/Scripts/constants.R") # Source file has fasta_file and metadata_file path.

### Quality filter:
# Make vectors to store the coverage of each sequence, and coverage and ID of low coverage sequences.
nt_perc <- c()
low_quality_seqs_ID <- c()
low_quality_seqs_Nperc <- c()
# Loop through the fasta file.
for (i in 1:length(fasta_file)){
  # Rename all the headers to be only the accession ID.
  names(fasta_file)[i] <- unlist(strsplit(names(fasta_file)[i], "|", fixed = TRUE))[2]
  # Get the percentage of N of each sequence and the percente of nt.
  N_perc <- (str_count(fasta_file[[i]], pattern = "N") * 100) / nchar(fasta_file[[i]])
  nt_perc <- c(nt_perc, (100 - N_perc))
  # If N percentage is higher than 20, get the value and the Id to separate vectors.
  if (N_perc > 20){
    low_quality_seqs_ID <- c(low_quality_seqs_ID, names(fasta_file)[i])
    low_quality_seqs_Nperc <- c(low_quality_seqs_Nperc, N_perc)
  }
}
# Make a dataframe to plot the coverage of all IDs
df_plot_quality <- data.frame(names(fasta_file), nt_perc)
colnames(df_plot_quality) <- c("Accession_ID", "Coverage")
# Plot barplot
quality_plot <- ggplot(df_plot_quality, aes(x = Accession_ID, y = Coverage)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_hline(yintercept = 80, color = "red") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.y = element_text(),
        axis.text.y = element_text(size = 9),
        axis.ticks.y = element_line(),
        axis.title.x = element_text(),
        axis.text.x = element_text(angle = 90, size = 0.8)) +
  ggtitle("Quality filter: Coverage > 80%") +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  labs(caption = paste0("Sequences excluded from analysis: ", length(low_quality_seqs_Nperc)))

quality_plot
ggsave("Coverage_barplot", quality_plot, device = "png", dpi = 700, bg = "white")

# Remove sequences of low quality from fasta file and metadata file
fasta_file <- fasta_file[! names(fasta_file) %in% low_quality_seqs_ID]
metadata_file <- metadata_file[! metadata_file$Accession.ID %in% low_quality_seqs_ID, ]
lab_file <- lab_file[! lab_file$gisaid_epi_isl %in% low_quality_seqs_ID, ]

#Change names of Lineages
metadata_file$Lineage <- gsub(" (marker override based on Emerging Variants AA substitutions)", "", metadata_file$Lineage, fixed = TRUE)
# Change all names of labs based on the lab tags file.
for (i in 1:nrow(lab_file)){
  if (lab_file$submitting_lab[i] %in% labtags_file$Laboratory_name){
    lab_file$submitting_lab[i] <- labtags_file$lab_tag[labtags_file$Laboratory_name == lab_file$submitting_lab[i]]
  }
}
lab_file$submitting_lab <- gsub("Microbiology,  Beaumont Hospital", "Beaumont Hospital", lab_file$submitting_lab, fixed = TRUE)
### Lineage Pie Chart:
# Get frequency of each lineage in a table and remove levels of lineage with as.character.
lineage_freq <- as.data.frame(table(metadata_file$Lineage))
lineage_freq$Var1 <- as.character(lineage_freq$Var1)
# Change column names and reorder table so high frequent lineages are on first rows.
colnames(lineage_freq) <- c("Lineage", "Frequency")
lineage_freq <- arrange(lineage_freq,desc(by_group = Frequency))

# Create empty vectors to store the top occurring lineages and their frequencies.
lineage <- c()
freq <- c()
# Get a count for the frequency of the less common lineages.
others_count <- 0
other_lineages <- c()
other_freq <- c()
# Loop through the number of rows.
for (i in 1:nrow(lineage_freq)){
  # The first five lineages get appended directly.
  if (i < 11){
    lineage <- c(lineage, lineage_freq$Lineage[i])
    freq <- c(freq, lineage_freq$Frequency[i])
  } else{
    # From the forth lineage onwards, increment the count.
    others_count <- others_count + lineage_freq$Frequency[i]
    other_lineages <- c(other_lineages, lineage_freq$Lineage[i])
    other_freq <- c(other_freq, lineage_freq$Frequency[i])
  }
}

# Append the vectors with the final count and the tag "others".
lineage <- c(lineage, "Others*")
freq <- c(freq, others_count)

# Make data frame for plotting with updated values.
df_plot <- data.frame(lineage, freq)

# Make the pie chart.
lineages_piechart <- ggplot(df_plot, aes(x = "", y = freq, fill = lineage)) + 
  geom_col() +
  geom_label(aes(label = freq),
             colour = "black",
             size = 2.8,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) + # Position of lineage frequency within pie.
  coord_polar(theta = "y") +
  guides(fill = guide_legend(title = "Lineage")) +
  theme_void() + # Remove grid and background lines.
  scale_fill_brewer(palette = "Spectral",
                    breaks = as.character(df_plot$lineage)) + # Rearrange legend on decreasing freq.
  ggtitle("Top 10 most occurring lineages") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  labs(caption = (paste0("*Lineages with less than ", max(other_freq), " sequences are grouped together"))) +
  theme(plot.margin = unit(c(10,5,10,5), "mm"))
lineages_piechart
ggsave("Top10_lineages", lineages_piechart, device = "png", dpi = 700, bg = "white")

### Origin lab piechart
# Get frequency of each lab in a table and remove levels of labs with as.character.
lab_freq <- as.data.frame(table(lab_file$submitting_lab))
lab_freq$Var1 <- as.character(lab_freq$Var1)
# Change column names and reorder table so high frequent labs are on first rows.
colnames(lab_freq) <- c("Laboratory", "Frequency")
lab_freq <- arrange(lab_freq,desc(by_group = Frequency))
# Create empty vectors to store the top occurring labs and their frequencies.
lab <- c()
freq <- c()
# Get a count for the frequency of the less common labs
others_count <- 0
others_name <- c()
others_freqs <- c()
# Loop through the number of rows.
for (i in 1:nrow(lab_freq)){
  # The first five lineages get appended directly.
  if (i < 10){
    lab <- c(lab, lab_freq$Laboratory[i])
    freq <- c(freq, lab_freq$Frequency[i])
  } else{
    # From the forth lineage onwards, increment the count and store the other lab names
    others_count <- others_count + lab_freq$Frequency[i]
    others_name <- c(others_name, lab_freq$Laboratory[i])
    others_freqs <- c(others_freqs, lab_freq$Frequency[i])
  }
}
# Append the vectors with the final count and the tag "others".
lab <- c(lab, "Others*")
freq <- c(freq, others_count)
# Make data frame for plotting with updated values.
df_plot <- data.frame(lab, freq)
df_plot
# Drop last row if others= 0
if (df_plot[df_plot$lab == "Others*",]$freq == 0){
  df_plot <- df_plot[-nrow(df_plot),]
}

# Make the pie chart.
lab_piechart <- ggplot(df_plot, aes(x = "", y = freq, fill = lab)) + 
  geom_col() +
  geom_label(aes(label = freq),
             colour = "black",
             size = 2.8,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) + # Position of lineage frequency within pie.
  coord_polar(theta = "y") +
  guides(fill = guide_legend(title = "Laboratory")) +
  theme_void() + # Remove grid and background lines.
  scale_fill_brewer(palette = "Spectral",
                    breaks = as.character(df_plot$lab)) + # Rearrange legend on decreasing freq.
  ggtitle("Submitting Laboratories") +
  #theme(plot.margin = unit(c(20,10,20,10), "mm")) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
  theme(plot.margin = unit(c(10,5,10,5), "mm"))

lab_piechart
ggsave("Submitting_laboratories", lab_piechart, device = "png", dpi = 700, bg = "white")

### Weeks of collection Bar Chart.
# Get the collection dates from the metadata file.
coll_dates <- as.Date(metadata_file$Collection.date)
# Create vectors to store each date and the week of date.
coll_date <- c()
coll_week <- c()
for (i in 1:length(coll_dates)){
  if (i == 1){
    coll_date <- coll_dates[i]
    coll_week <- week(coll_dates)[i]
  } else{
  coll_date <- c(coll_date, coll_dates[i])
  coll_week <- c(coll_week, week(coll_dates)[i])
  }
}

# Bind the two vector in data frame.
sorted_dates <- data.frame(sort(coll_date), sort(coll_week), stringsAsFactors = FALSE)

# Get the frequency of each week.
weeks_freq <- as.data.frame(table(sorted_dates$sort.coll_week.))

# Compare the frequency table with the sorted dates to get the first day of each week.
week_first_day <- c()
for (i in 1:nrow(weeks_freq)){
  if (i == 1){
    week_first_day <- c(sorted_dates[sorted_dates$sort.coll_week. == weeks_freq$Var1[i],][1,1])
  } else{
    week_first_day <- c(week_first_day, sorted_dates[sorted_dates$sort.coll_week. == weeks_freq$Var1[i],][1,1])
  }
}

# Add a column to the frequency df with the first date of each week.
weeks_freq <- cbind(weeks_freq, week_first_day)
colnames(weeks_freq) <- c("Week", "N_Sequences", "Week_first_day")

weeks_of_collection_barplot <- ggplot(weeks_freq, aes(x = Week, y = N_Sequences, fill = N_Sequences)) +
  geom_bar(stat = "identity", color = "black") +
  geom_label(aes(label = N_Sequences),
             colour = "black",
             size = 3,
             position = position_stack(vjust = 1),
             show.legend = FALSE) +
  theme_minimal() +
  scale_fill_viridis_b() +
  #xlab("Week") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 5, size = 9, face = "bold")) +
  ggtitle("Number of sequences collected per week") +
  scale_x_discrete(labels = week_first_day) +
  theme(plot.title = element_text(size = 18, face = "bold")) +
  theme(plot.margin = unit(c(10,5,10,5), "mm"))
weeks_of_collection_barplot
ggsave("Weeks_of_collection", weeks_of_collection_barplot, device = "png", dpi = 700, bg = "white")

### County map
# Get the frequency of each county.
county_freq <- as.data.frame(table(metadata_file$Location))
county_freq$Var1 <- as.character(county_freq$Var1)
colnames(county_freq) <- c("County", "Frequency")
for (i in 1:nrow(county_freq)){
  county_freq$County[i] <- unlist(strsplit(county_freq$County[i], "/ "))[3]
}

# Get the names of counties on order from the map file.
map_counties <- c()
for (i in 1:nrow(ROI_map)){
  current_county <- as.character(ROI_map["COUNTYNAME"]$COUNTYNAME[i])
  map_counties <- c(map_counties, unlist(strsplit(current_county, " "))[1])
}

# Make a dataframe with the counties and a count number starting at 0, and the county ID of the map.
county_count <- numeric(length(map_counties))
df_county_count <- data.frame(map_counties, county_count, ROI_map["COUNTY"]$COUNTY)

# Fill the county count data frame with the frequency of each county.
for (i in 1:nrow(county_freq)){
  current_county <- county_freq$County[i]
  if (current_county %in% df_county_count$map_counties){
    df_county_count[df_county_count$map_counties == current_county,]$county_count <- county_freq$Frequency[i]
  }
}

# Make a modified map with the frequency of each county.
modified_map <- ROI_map %>%
  left_join(df_county_count, by = c("COUNTY" = "ROI_map..COUNTY...COUNTY"))

roi_map_plot <- ggplot(modified_map) +
  geom_sf(aes(fill = county_count), color = "black", size = 0.1) + # fill counties with sequences frequency 
  geom_sf_text(aes(label = map_counties), colour = "white", fontface = "bold", size = 1.78) + # add name of each county
  theme_void() + # remove grid, background and axis.
  scale_fill_viridis_c(name = "") +
  ggtitle(paste0("Sequences per county. Total sequences: ", nrow(metadata_file))) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
roi_map_plot
ggsave("Sequences_per_county", roi_map_plot, device = "png", dpi = 700, bg = "white")

### Sequences phiylogeny
# Run Mafft for multiple seq alignment and raxml for phylogeny.
system(mafft_command)
system(raxml_command)

# Load the file as a tree and rename the tib labels.
tree <- read.tree(file = "RAxML_bipartitionsBranchLabels.raxml_tree")
for (i in 1:length(tree$tip.label)){
  tree$tip.label[i] <- unlist(strsplit(tree$tip.label[i], "|", fixed = TRUE))[2]
}

#tree_text <- paste(readLines("RAxML_bipartitionsBranchLabels.raxml_tree"), collapse="\n")
#tree <- read.nhx(textConnection(treetext))

# Make a data frame to store the submitting lab of each ID.
tree_df <- data.frame(lab_file$gisaid_epi_isl, lab_file$submitting_lab)
# Get the lineage of each ID and add it as an additional column to the data frame.
tree_lineages <- c()
tree_counties <- c()
for (i in tree_df$lab_file.gisaid_epi_isl){
  tree_lineages <- c(tree_lineages, metadata_file[metadata_file$Accession.ID == i, "Lineage"])
  tree_counties <- c(tree_counties, unlist(strsplit(metadata_file[metadata_file$Accession.ID == i, "Location"], "/ "))[3])
  unlist(strsplit(county_freq$County[i], "/ "))[3]
}
tree_df <- cbind(tree_df, tree_lineages, tree_counties)
colnames(tree_df) <- c("Accession_ID", "Laboratory", "Lineage", "County")
# Make the tree and link it with the tree dataframe to colour tips
plot_tree <- ggtree(tree) %<+% tree_df +
  geom_tiplab(aes(fill = factor(County)),
              color = "black",
              geom = "label",
              size = 1.1,
              label.padding = unit(0.02, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label borders
  guides(fill = guide_legend(title = "County")) +
  geom_tippoint(aes(color = factor(Lineage)),
                size = 0.7,
                alpha = 0.85) +
  guides(color = guide_legend(title = "Lineage"))
plot_tree
ggsave("Tree_counties", plot_tree, device = "png", height = 420, width = 420, units = "mm", dpi = 700, )

#
#raxml <- read.raxml("RAxML_bipartitionsBranchLabels.raxml_tree")
#for (i in raxml[,"bootstrap"]){
#  print(i)
#}

#bootstrap_tree <- ggtree(raxml) +
#  geom_label(aes(label=bootstrap, fill=bootstrap)) +
#  scale_fill_viridis_c() +
#  theme(legend.position='top',
#        legend.justification='left',
#        legend.direction='horizontal')
#bootstrap_tree

### Variants mutations
# 1-TARGET VARIANTS
# Create a list to store the result where each element will be the row number of the metadata file presenting the targeted mutations.
targets_result <- list()
# Loop through the indices of the target mutations
for (i in 1:length(target_mutations)){
  # Get all the row numbers of the metadata file.
  candidate_rows <- 1:length(metadata_file$AA.Substitutions)
  # Loop through the elements of each entry of the list.
  for (j in target_mutations[[i]]){
    # Get all the rows of the metadata file whose mutation section contains the target mutations of each entry. 
    current_rows <- which(grepl(j, metadata_file$AA.Substitutions))
    candidate_rows <- intersect(current_rows, candidate_rows)
  }
  targets_result[[i]] = candidate_rows
}
# Create a variable to store the final table with three columns being accession ID, lineage and list of mutations.
target_mut_results <- c()
# Loop through the number of entries of the target mutations.
for (i in 1:length(target_mutations)){
  # Get the Accession ID ad Lineage corresponding to each row number.
  col1_2 <- metadata_file[targets_result[[i]], c("Accession.ID", "Lineage")]
  # Get the mutations of that entry pasted together as many times as the length of the results of that entry and bind it with the first columns.
  col3 <- rep(paste(target_mutations[[i]], collapse=", "), length(targets_result[[i]]))
  result_columns <- cbind(col1_2, col3)
  # Bind the three columns to the results table.
  target_mut_results <-rbind(target_mut_results, result_columns)
}

# Loop through the accession IDs
for (i in 1:nrow(target_mut_results)){
  current_ID <- target_mut_results$Accession.ID[i]
  # If the Id is repeated.
  if (length(target_mut_results[target_mut_results$Accession.ID == current_ID,]) > 1){
    # Store the entries with repeated ID.
    repeated <- target_mut_results[target_mut_results$Accession.ID == current_ID,]
    # Remove the rows whose mutations are shorter, so the remaining one is the on with all mutations.
    condition_to_remove <- which(nchar(repeated$col3) != max(nchar(repeated$col3)))
    rownames_to_remove <- rownames(repeated[c(condition_to_remove), ])
    target_mut_results <- target_mut_results[!(rownames(target_mut_results) %in% rownames_to_remove),]
  }
}
# Remove unnecessary strings and change column names
colnames(target_mut_results) <- c("Accession.ID", "Lineage", "Mutations")

# Export as png
png("Target_mutations.png", height = 50*nrow(target_mut_results), width = 200*ncol(target_mut_results))
grid.table(target_mut_results)
dev.off()

### Lineage percentages
# Get lineages and collection dates in a data frame.
lineages <- metadata_file$Lineage
lin_dates <- as.Date(metadata_file$Collection.date)
lineage_dates <- data.frame(lineages, lin_dates)
lineage_dates <- na.omit((lineage_dates))
# Fill a vector with the week corresponding to each date.
lin_weeks <- c()
for (i in 1:nrow(lineage_dates)){
  lin_weeks <- c(lin_weeks, week(lineage_dates$lin_dates[i]))
}
lineage_dates <- data.frame(lineage_dates, lin_weeks)
colnames(lineage_dates) <- c("lineage", "date", "week")
# Get the unique lineages and unique weeks
unique_lineages <- sort(unique(lineage_dates$lineage))
unique_weeks <- sort(unique(lineage_dates$week))
# Make a dataframe with weeks and lineages as columns and rows names.
lindates_matrix <- data.frame(matrix(nrow = length(unique_lineages) + 1, ncol = length(unique_weeks)))
colnames(lindates_matrix) <- c(unique_weeks)
rownames(lindates_matrix) <- c(unique_lineages, "Total")
# Loop through the column names
for (i in colnames(lindates_matrix)){
  # The total of samples for each week is extracted by subseting the lineage dates with the current weeks and counting the row number.
  column_total = nrow(lineage_dates[lineage_dates$week == i,])
  lindates_matrix["Total", i] <- column_total
  # Loop through the row names
  for (j in rownames(lindates_matrix)){
    # The value for each cell is extracted by subseting the lineage dates with the current weeks and lineage and counting the row number.
    value <- nrow(lineage_dates[lineage_dates$lineage == j & lineage_dates$week == i,])
    # Make the percentage by multiplying 100 and dividing by the total
    percent_val <- (value * 100) / column_total
    lindates_matrix[j, i] <- round(percent_val, 2)
  }
  # Make the total number of each week the row total instead of the calculated percentage.
  lindates_matrix["Total",i] <- column_total
}
# Format the name of the columns.
for (i in 1:ncol(lindates_matrix)){
  colnames(lindates_matrix)[i]<- paste0("Week ", colnames(lindates_matrix)[i])
}
# Get the total occurence of each lineages form the lineage dates dataframe. 
total_lineages <- c()
for (i in rownames(lindates_matrix)){
  total_lineages <- c(total_lineages, nrow(lineage_dates[lineage_dates$lineage == i,]))
}
# Add the total lineages as a colum to the lindates matrix.
total_lineages <- c(total_lineages)
lindates_matrix$Total <- total_lineages

png("Lineage_percentages", height = 50*nrow(lindates_matrix), width = 200*ncol(lindates_matrix))
grid.table(lindates_matrix)
dev.off()

# Create PDF and save all the plots.
#pdf(paste0("SARS-CoV-2_NVRLreport_", Sys.Date(), ".pdf", height=11, width=10))
#print(quality_plot)
#print(lineages_piechart)
#print(lab_piechart)
#print(weeks_of_collection_barplot)
#print(roi_map_plot)
#print(phylogenetic_tree)
#othertree
#grid::grid.newpage()
#grid.table(target_mut_results)
#grid::grid.newpage()
#grid.table(lindates_matrix)
#dev.off()









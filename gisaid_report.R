source("constants.R") # Source file has fasta_file and metadata_file path.

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
# Loop through the number of rows.
for (i in 1:nrow(lineage_freq)){
  # The first five lineages get appended directly.
  if (i < 6){
    lineage <- c(lineage, lineage_freq$Lineage[i])
    freq <- c(freq, lineage_freq$Frequency[i])
  } else{
    # From the forth lineage onwards, increment the count.
    others_count <- others_count + lineage_freq$Frequency[i]
  }
}
# Append the vectors with the final count and the tag "others".
lineage <- c(lineage, "Others")
freq <- c(freq, others_count)

# Make data frame for plotting with updated values.
df_plot <- data.frame(lineage, freq)

# Make the pie chart.
piechart <- ggplot(df_plot, aes(x = "", y = freq, fill = lineage)) + 
  geom_col() +
  geom_label(aes(label = freq),
             colour = "black",
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) + # Position of lineage frequency within pie.
  coord_polar(theta = "y") +
  guides(fill = guide_legend(title = "Lineage")) +
  theme_void() + # Remove grid and background lines.
  scale_fill_brewer(palette = "Spectral",
                    breaks = as.character(df_plot$lineage)) + # Rearrange legend on decreasing freq.
  ggtitle("Most occurring lineages") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
  
piechart

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
sorted_dates <- data.frame(sort(coll_date), sort(coll_week))
# Get the frequency of each week.
weeks_freq <- as.data.frame(table(sorted_dates$sort.coll_week.))
# Compare the frequency table with the sorted dats to get the first day of each week.
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

brplot <- ggplot(weeks_freq, aes(x = Week_first_day, y = N_Sequences, fill = N_Sequences)) +
  geom_bar(stat = "identity", color = "black") +
  geom_label(aes(label = N_Sequences),
             colour = "black",
             position = position_stack(vjust = 1),
             show.legend = FALSE) +
  theme_minimal() +
  scale_fill_viridis_b() +
  #xlab("Week") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 1.5, size = 9)) +
  ggtitle("Number of sequences collected per week") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
  #scale_fill_brewer(palette = "Spectral")

brplot

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

map_plot <- ggplot(modified_map) +
  geom_sf(aes(fill = county_count), color = "black", size = 0.1) + # fill counties with sequences frequency 
  geom_sf_text(aes(label = map_counties), colour = "white", fontface = "bold", size = 1.78) + # add name of each county
  theme_void() + # remove grid, background and axis.
  scale_fill_viridis_c(name = "") +
  ggtitle(paste0("Number of sequences per couny. Total sequences: ", nrow(metadata_file)))
map_plot

#print(piechart)
### Sequences philogeny
# Run Mafft for multiple seq alignment.
# Documentation: https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
#system(paste0('mafft --auto --reorder ', fasta_name, ' > "mafft_aligned.fasta"'))
# Run raxml-ng for phylogeny
# Documentation: https://biohpc.cornell.edu/doc/alignment_exercise3.html
#system('raxml-ng --msa mafft_aligned.fasta --model GTR+G')
tree_path = "/home/gabriel/Desktop/Jose/Projects/GISAID_reports/Data/August/A1.raxml.bestTree"
tree_file <- read.tree(file = tree_path)

### Variants mutations
# Get all the mutations and remove parenthesis.
mutations <- str_remove_all(metadata_file$AA.Substitutions, "[()]")

# Make a list and fill it with all the unique mutations.
mutation_list <- c()
# Unlist every list of mutations for each seq.
for (i in 1:length(mutations)){
  current_mut <- unlist(strsplit(mutations[i], ","))
  # Add every mutation to the mutation list if it is not present yet.
  for (j in 1:length(current_mut)){
    if (!(current_mut[j] %in% mutation_list)){
      mutation_list <- c(mutation_list, current_mut[j])
    }
  }
}

# Create data frame with unique lineages as column names and unique mutations as column names.
mutations_df <- data.frame(matrix(ncol = length(unique(metadata_file$Lineage)), nrow = length(mutation_list)))
colnames(mutations_df) <- unique(metadata_file$Lineage)
rownames(mutations_df) <- mutation_list
# Fill all na cells with 0.
mutations_df[is.na(mutations_df)] <- 0

# Loop through the metadata file and fill in the mutations df.
for (i in metadata_file$Lineage){
  for (j in unlist(strsplit(mutations, ","))){
    mutations_df[j,i] <- mutations_df[j,i] + 1
    
  }
}

for (j in unlist(strsplit(mutations, ","))){
  print(j)
}






# Get the top 5 most occuring lineages.
top5_lin <- lineage_freq$Lineage[1:5]











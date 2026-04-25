# AUTHOR / EMAIL: Bethany Mariel Suliguin (bs7902@princeton.edu) | Modified code from Michelle White and Stavi Tennenbaum
# CREATION DATE: 03/23/26
# MODIFIED DATE: 04/15/26
# PURPOSE: up- and down-valley analysis (PCA, FROH)

############################################################################
# read in installed packages
library(readxl) # for reading Excel files
library(dplyr) # for data manipulation
library(ggplot2) # for plotting
library(wesanderson) # to get some fun colors
library(readxl) # for reading Excel files
library(stringr) # used to work with strings
library(cowplot) # for combining plots
library(sf) # for map shapefiles and map overlays
library(ggspatial) # for making maps with ggplot
library(maptiles) # for maptiles 
############################################################################

# set working directory
setwd("/Users/bethanysuliguin/Documents/EEB331")

# set up plot aesthetics
eeb_theme <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Set all text to Arial
      text = element_text(family = "Arial"),
      # Remove grid lines and background
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      # Plot title
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, margin = margin(b = 10)),
      # Axes
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 9),
      # Legend
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_blank(),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9),
    ) 
}

# Read in marmot metadata
metadata <- read_excel("radseq_pedigree_metadata.xlsx")

# add column that shows up/down classification (taken from col_area names) 
metadata <- metadata |> 
  mutate(col_class = case_when(
    grepl("river", col_area) ~ "down",
    grepl("mm", col_area) ~ "up",
    grepl("picnic", col_area) ~ "up",
    grepl("cliff", col_area) ~ "up",
    col_area %in% c("bench", "horsemound", "rvannex", "gothictown", "avalanche") ~ "down",
    col_area %in% c("boulder", "northpk") ~ "up"
  ))

# read in minimally filtered PCA files
pve = read.table("pve.txt", header=F)
head(pve) # PC1: 6.8%, PC2: 5.2%

pcs = read.table("pcs.txt", header=T)

# find how many marmots in pcs
dim(pcs) # 221 marmots

############################################################################

# PCA analysis

# remove invalid IDs and unwanted data; make IDs match UID format for future merging
pcs <- pcs |> 
  select(-FID) |> 
  rename(uid = 1) |> 
  slice(-(1:3)) |> 
  mutate(uid = gsub("^UID_|_sorted$", "", uid)) # find the beginning and end of each string -> replace with whitespace in the column 'uid'

# combine mm_ colonies to make "marmot meadow"; river_ colonies to make "river"; picnic_ to make "picnic"
metadata <- metadata %>%
  mutate(col_area = gsub("_(lower|upper|middle)$", "", col_area)) %>%
  mutate(col_area = gsub("_(aspen|maintalus)$", "", col_area)) %>%
  mutate(col_area = gsub("_(rivermound|sagemound)$", "", col_area))

# merge new PCA data with new metadata
merged_data <- merge(pcs, metadata, by.x = "uid", by.y = "uid")

# count how many marmots in either up or down valley
table(merged_data$col_class) # 72 down, 146 up

# PCA by Valley Location
col_type <- ggplot(merged_data, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill=col_class), color = "black", size = 3, pch = 21) + 
  labs(title = "PCA by Valley Location", # rename and add PC %s
       x = "PC 1 (6.8%)",
       y = "PC 2 (5.2%)",
       fill = "Valley Location") +
  eeb_theme() +
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))

# ANOVAs: PCs with up/down valley
col_type_aov_pc1 <- aov(merged_data$PC1 ~ as.factor(merged_data$col_class)) 
summary(col_type_aov_pc1)  # F = 137.5, p = <2e-16 ***
col_type_aov_pc2 <- aov(merged_data$PC2 ~ as.factor(merged_data$col_class))  
summary(col_type_aov_pc2)  # F = 65.15, p = 4.8e-14 ***

############################################################################

# FROH Analysis

# constants
genome_size <- 2319465413 # in nucleotide base pairs, I cal'd this for you!

# read in the ROH output from bcftools
lines <- readLines("Marmots_maf3_geno90_mind20_ROH.txt") 
lines <- lines[!grepl("^#", lines)]
rg_lines <- lines[grepl("^RG", lines)]

# extract ROH segments based on file structure
roh_segments <- read.table(
  text = rg_lines,
  header = FALSE,
  stringsAsFactors = FALSE
)

# add back in column names based on read group (RG) headers
colnames(roh_segments) <- c(
  "type",
  "sample",
  "chr",
  "start",
  "end",
  "length_bp",
  "n_markers",
  "quality"
)

# calculate some basic metrics
roh_summary <- roh_segments %>%
  mutate(length_bp = as.numeric(end) - as.numeric(start)) %>%
  group_by(sample) %>%
  summarise(
    NSEG = n(),
    total_bp = sum(length_bp),
    mean_bp = mean(length_bp),
    mean_length = total_bp / NSEG
  )

# calculate FROH per individual
roh_summary <- roh_summary %>%
  mutate(FROH = total_bp / genome_size)

# clean roh_summary columns for future merging
roh_clean <- roh_summary |> 
  slice(-221) |> 
  slice(-(1:2)) |> 
  mutate(uid = gsub("UID_|_sorted.*$", "", sample)) # find the beginning and end of each string -> replace with whitespace in the column 'uid'

# merge roh_clean and metadata
merged_roh <- merge(roh_clean, merged_data, by.x = "uid", by.y = "uid")

# plot FROH by up/down valley
valley_FROH <- ggplot(merged_roh, aes(x = col_class, y = FROH)) +
  geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
  geom_boxplot(fill = NA, color = 'black') +
  labs(x = "Colony Classification", y = "FROH", title = "Distribution of FROH by Valley Location") +
  eeb_theme() +
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))

# t-test/statistics for FROH by colony classification (up/down)
t.test(FROH ~ col_class, data = merged_roh)
# t = 0.35656, df = 128.36, p-value = 0.722
# mean in down-valley = 0.2446503; mean in up-valley = 0.2409272

# combine PCA and FROH into one row for a figure
combined_valley_analysis <- plot_grid(col_type, valley_FROH, ncol = 2, labels = c("a", "b"))

# save multiplot
ggsave("Up and Down Valley Multiplot.png", combined_valley_analysis, width = 10, height = 5, units = "in", dpi=1600, bg="white")

############################################################################

# make a map to show where colonies are

# load in map
shp_map <- st_read("/Users/bethanysuliguin/Documents/EEB331/marmot_polygons_wgs84.shp")
shp_map <- st_transform(shp_map, 3857)

# remove unrepresented colonies
shp_map <- shp_map |> 
  filter(Site != "Cliff") |> 
  filter(Site != "Avery") |> 
  filter(Site != "Bellview") |> 
  filter(Site != "Stonefield")

# load in satellite image
esri_tiles <- get_tiles(shp_map, provider = "Esri.WorldImagery", zoom = 16)

# make the plot
colony_map <- ggplot() +
  layer_spatial(esri_tiles) +
  geom_sf(data = shp_map, fill = NA, color = "red", linewidth = 0.5) +
  coord_sf(expand = FALSE) +
  geom_sf_text(
    data = shp_map, 
    aes(label = Site), 
    color = 'white',
    size = 1.5,
    nudge_x = -200,
    nudge_y = 100) +
  eeb_theme() +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

ggsave("Colony Plot.png", colony_map, width = 7, height = 7, units = "in", dpi=1600, bg="white")

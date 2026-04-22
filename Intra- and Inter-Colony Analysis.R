##########Part 1: Group Based analyses of FROH (intra- and inter-colony)
### Made by: Marina Cookson (mc3861@princeton.edu)
### Date finished and checked: April 21, 2026

###################################################################
################### Loading Packages being used ###################
###################################################################

library(readr) # for reading files
library(dplyr) # for data manipulation
library(ggplot2) # for plotting
library(wesanderson) # to get some fun colors
library(cowplot) # for plotting multiple plots in one
install.packages("readxl")
library(readxl)

# set up plot aesthetics we can use for all figures
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
###################################################################
################### Loading Dataset ###################
###################################################################

setwd("/Users/marinacookson/Desktop/EEB331")

###Read in data 
MYDATA <- read_excel("radseq_pedigree_metadata.xlsx")
View(MYDATA) # take a look

###################################################################
################### Running PCA and graphing ###################
###################################################################

pve = read.table("pve.txt", header=F)
head(pve) # record % variance explained on PC1 and PC2, control F to add updated values to all PCA plots!
pcs = read.table("pcs.txt", header=T)


#### data cleaning #### 
# we can see here that column 1 is redundant and the first 3 rows don't have valid IDs, so let's drop them
# let's also edit the column of unique animal IDs to match the "uid" format of our metadata for merging
pcs <- pcs %>%
  select(-FID) %>%
  rename(uid = 1) %>%
  slice(-(1:3)) %>%
  mutate(uid = gsub("^UID_|_sorted$", "", uid)) #find the beginning and end of each string -> replace with whitespace in the column 'uid'

View(pcs)

# we'll also need to clean up our metadata a bit before we proceed
# we can ignore some of the columns for now, like 'furmark', but no need to drop them
# picnic_lower and picnic_upper are the same colony, so edit this be replacing all white space after picnic_
# follow the same logic to clean up the names for the mm_ and river_ colonies
MYDATA <- MYDATA %>%
  mutate(col_area = gsub("_(lower|upper|middle)$", "", col_area)) %>%
  mutate(col_area = gsub("_(aspen|maintalus)$", "", col_area)) %>%
  mutate(col_area = gsub("_(rivermound|sagemound)$", "", col_area))

View(MYDATA)

# now MERGE new PCA data with new metadata so we can easily make updated plots!
MYDATA2 <- merge(pcs, MYDATA, by.x = "uid", by.y = "uid")
View(MYDATA2)

#### merging river and river southmound under RIVER
MYDATA2$col_area[MYDATA2$col_area %in% c("river", "river_southmound")] <- "river"

#### Plotting PCA with elypses
colony <- ggplot(MYDATA2, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = col_area),
             color = "black",
             size = 2.5,
             pch = 21,
             alpha = 0.8) +
  stat_ellipse(aes(color = col_area),
               linewidth = 0.8,
               show.legend = FALSE) +
  labs(title = "Colony Area",
       x = "PC 1 (6.8%)",
       y = "PC 2 (5.2%)",
       fill = "Colony") +
  scale_fill_discrete(labels = c(
    "avalanche" = "Avalanche",
    "bench" = "Bench",
    "boulder" = "Boulder",
    "cliff" = "Cliff",
    "gothictown" = "Gothic Town",
    "horsemound" = "Horsemound",
    "mm" = "Marmot Meadow",
    "northpk" = "North Picnic",
    "picnic" = "Picnic",
    "river" = "River",
    "rvannex" = "River Annex" 
  )) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  eeb_theme()

colony

### Checking numbers of individuals in each colony 
MYDATA2 %>%
count(col_area)

#### Graphing the number of individuals in each colony 
count_plot <- ggplot(MYDATA2, aes(x = col_area, fill = col_area)) +
  geom_bar() +
  labs(x = "Colony", y = "Count", title = "Number of Individuals Per Colony")+
  scale_fill_manual(values = c(
    "avalanche" = "#F8766D",
    "bench" = "#D89000",
    "boulder" = "#A3A500",
    "cliff" = "#39B600",
    "gothictown" = "#00BF7D",
    "horsemound" = "#00BFC4",
    "mm" = "#00B0F6",
    "northpk" = "#619CFF",
    "picnic" = "#9590FF",
    "river" = "#E76BF3",
    "rvannex" = "#FF62BC"
  )) +
  scale_x_discrete(labels = c(
    "avalanche" = "Avalanche",
    "bench" = "Bench",
    "boulder" = "Boulder",
    "cliff" = "Cliff",
    "gothictown" = "Gothic Town",
    "horsemound" = "Horsemound",
    "mm" = "Marmot Meadow",
    "northpk" = "North Picnic",
    "picnic" = "Picnic",
    "river" = "River",
    "rvannex" = "River Annex"
  ))+
  eeb_theme() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

count_plot

##########################
#### FST #################
#########################

setwd("/Users/marinacookson/Desktop/EEB331")

# Install packages
install.packages("adegenet")
install.packages("hierfstat")
install.packages("ggplot2")
install.packages("reshape2")

# Load libraries
library(adegenet)
library(hierfstat)
library(ggplot2)
library(reshape2)

### LOAD GENOTYPE DATA
geno <- read.table("marmots_2k.raw", header = TRUE)

### Remove metadata columns (first 6 columns)
geno_clean <- geno[, -c(1:6)]

### Extract IDs
ids <- geno$IID

### Clean IDs
ids_clean <- gsub("^UID_", "", ids)
ids_clean <- gsub("_sorted$", "", ids_clean)

### Match metadata
idx <- match(ids_clean, MYDATA$uid)
MYDATA2 <- MYDATA[idx, ]

### CLEAN COLONY NAMES 
# robust replacement for any "river_south..." variation
MYDATA2$col_area <- gsub("river.*south.*", "river", MYDATA2$col_area)

# check it worked
unique(MYDATA2$col_area)

### REMOVE NA
valid <- !is.na(MYDATA2$col_area)

MYDATA2 <- MYDATA2[valid, ]
geno_clean <- geno_clean[valid, ]

### Remove loci with NA
geno_clean2 <- geno_clean[, colSums(is.na(geno_clean)) == 0]

### Create population vector
pop <- MYDATA2$col_area

### Build hierfstat input
hf <- data.frame(pop = as.factor(pop), geno_clean2)

### Calculate pairwise FST
fst_matrix <- pairwise.WCfst(hf)

### Clean negative values
fst_matrix[fst_matrix < 0] <- 0

### View results
fst_matrix

### Heatmap (base R)
heatmap(as.matrix(fst_matrix))

### Prepare for ggplot
fst_df <- melt(as.matrix(fst_matrix))
colnames(fst_df) <- c("Colony1", "Colony2", "FST")

# remove diagonal
fst_df <- fst_df[fst_df$Colony1 != fst_df$Colony2, ]

### Recode names
fst_df$Colony1 <- recode(fst_df$Colony1,
                         "mm" = "Marmot Meadow",
                         "northpk" = "North Picnic",
                         "picnic" = "Picnic",
                         "river" = "River",
                         "rvannex" = "River Annex",
                         "gothictown" = "Gothic Town",
                         "horsemound" = "Horsemound",
                         "bench" = "Bench",
                         "boulder" = "Boulder",
                         "cliff" = "Cliff",
                         "avalanche" = "Avalanche")

fst_df$Colony2 <- recode(fst_df$Colony2,
                         "mm" = "Marmot Meadow",
                         "northpk" = "North Picnic",
                         "picnic" = "Picnic",
                         "river" = "River",
                         "rvannex" = "River Annex",
                         "gothictown" = "Gothic Town",
                         "horsemound" = "Horsemound",
                         "bench" = "Bench",
                         "boulder" = "Boulder",
                         "cliff" = "Cliff",
                         "avalanche" = "Avalanche")

### Plot heatmap
fst_plot <- ggplot(fst_df, aes(x = Colony1, y = Colony2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(FST, 3)), size = 1.5)+
  scale_fill_gradient(low = "steelblue1", high = "steelblue4") +
  theme_minimal() +
  eeb_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.line = element_blank(), 
    panel.grid = element_blank()
  ) +
  labs(
    title = "Pairwise FST Between Colonies",
    x = "Colony",
    y = "Colony",
    fill = "FST")

fst_plot

##################################
#### MANTEL TEST #################
##################################

install.packages("vegan")
library(vegan)
install.packages("sf")
library(sf)

##### create FST matrix character 
fst_mat <- as.matrix(fst_matrix)
rownames(fst_mat) ###checking info
colnames(fst_mat) ###checking info

###Get coordinate information for geographic matrix
shp <- st_read("marmot_polygons_wgs84.shp") ## make sure all 4 polygon files are in home directory file
names(shp)

###extract coordinates 
unique(shp$Site)  
unique(shp$Site2) ### this one has the names i want so im using this
rownames(fst_matrix) ### but the names arent exactly what i want 

#### fix names
shp$colony <- shp$Site2
shp$colony[shp$colony == "marmotmeadow"] <- "mm"
shp$colony[shp$colony == "northpicnic"] <- "northpk"

### check if everything matches
setdiff(rownames(fst_matrix), shp$colony) ## if character(0) is the output it is correct

### now extract coordinates
centroids <- st_centroid(shp)

centroids <- centroids[match(rownames(fst_matrix), shp$colony), ]

### make geographic distance matrix 
geo_mat <- as.matrix(st_distance(centroids))

### visualize matrix 
library(reshape2)
library(ggplot2)
library(units)

geo_mat <- as.matrix(geo_mat)
geo_mat <- units::drop_units(geo_mat)

geo_mat <- geo_mat / 1000

# Make sure matrix has names
rownames(geo_mat) <- rownames(fst_matrix)
colnames(geo_mat) <- colnames(fst_matrix)

# convert to dataframe
geo_df <- reshape2::melt(geo_mat)
colnames(geo_df) <- c("Colony1", "Colony2", "Distance")

# change names
name_map <- c(
  "avalanche" = "Avalanche",
  "bench" = "Bench",
  "boulder" = "Boulder",
  "cliff" = "Cliff",
  "gothictown" = "Gothic Town",
  "horsemound" = "Horsemound",
  "mm" = "Marmot Meadow",
  "northpk" = "North Picnic",
  "picnic" = "Picnic",
  "river" = "River",
  "rvannex" = "River Annex")

geo_df$Colony1 <- name_map[geo_df$Colony1]
geo_df$Colony2 <- name_map[geo_df$Colony2]

geo_df$Colony1 <- factor(geo_df$Colony1, levels = name_map[rownames(geo_mat)])
geo_df$Colony2 <- factor(geo_df$Colony2, levels = name_map[colnames(geo_mat)])

# plot
geo_plot <- ggplot(geo_df, aes(x = Colony1, y = Colony2, fill = Distance)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Distance, 1)), size = 2) +
  scale_fill_gradient(low = "orchid1", high = "orchid4") +
  theme_minimal() +
  eeb_theme() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.line = element_blank(), 
    panel.grid = element_blank()
  ) +
  labs(
    title = "Geographic Distance Between Colonies",
    x = "Colony",
    y = "Colony",
    fill = "Distance (km)")

geo_plot


######## RUNNING MANTEL TEST

# check row names before test = must be the same
rownames(fst_matrix)
rownames(geo_mat)

# make sure there are no NA values = output should be zero
sum(is.na(fst_matrix)) ## was not zero 
sum(is.na(geo_mat))

#replace NA with 0 
fst_matrix[is.na(fst_matrix)] <- 0

#check 
sum(is.na(fst_matrix))

# run test
mantel_result <- mantel(
  as.dist(fst_matrix),
  as.dist(geo_mat),
  method = "pearson",
  permutations = 999)     #999 is standard practice for randomness, basically means how many randomness simulations are done

mantel_result

####RESULTS 
    #Mantel statistic based on Pearson's product-moment correlation 
      #Call:
         #mantel(xdis = as.dist(fst_matrix), ydis = as.dist(geo_mat), method = "pearson",      permutations = 999) 

      #Mantel statistic r: 0.4806 
      #Significance: 0.002 

      #Upper quantiles of permutations (null model):
          # 90%   95% 97.5%   99% 
       # 0.158 0.221 0.290 0.329 
          #Permutation: free
          #Number of permutations: 999

### RUNNING IBD PLOT

#Transform FST to genetic distance
gen_dist <- fst_matrix / (1 - fst_matrix)

# remove problematic values (Inf / NaN if FST = 1 or NA)
gen_dist[!is.finite(gen_dist)] <- NA

#Match row/column order
# ensure same order
gen_dist <- gen_dist[rownames(geo_mat), colnames(geo_mat)]


#Convert to vectors (upper triangle only)
gen_vec <- as.vector(gen_dist[upper.tri(gen_dist)])
geo_vec <- as.vector(geo_mat[upper.tri(geo_mat)])


### Remove NA pairs 
valid <- !is.na(gen_vec) & !is.na(geo_vec)

gen_vec <- gen_vec[valid]
geo_vec <- geo_vec[valid]


### Basic IBD plot
plot(geo_vec, gen_vec,
     xlab = "Geographic Distance",
     ylab = "Genetic Distance (FST / (1 - FST))",
     main = "Isolation by Distance",
     pch = 19)

# add regression line
abline(lm(gen_vec ~ geo_vec), col = "red", lwd = 2)

#replotting 
library(ggplot2)

ibd_df <- data.frame(
  Geographic = geo_vec,
  Genetic = gen_vec
)

ibd_plot <- ggplot(ibd_df, aes(x = Geographic, y = Genetic)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Isolation by Distance",
    x = "Geographic Distance",
    y = "Genetic Distance (FST / (1 - FST))"
  ) + 
  eeb_theme() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank()
  )


ibd_plot

######################################
########## AMOVA #####################
######################################

install.packages("poppr")
library(poppr)
library(adegenet)

# convert names
colnames(geno_clean) <- paste0("L", 1:ncol(geno_clean))

#convert genotypes
convert_geno <- function(x){
  ifelse(is.na(x), NA,
         ifelse(x == 0, "A/A",
                ifelse(x == 1, "A/B",
                       ifelse(x == 2, "B/B", NA))))}

geno_allele <- as.data.frame(lapply(geno_clean, convert_geno))

# convert to proper adagenet syntax
genind_obj <- df2genind(
  geno_allele,
  sep = "/",
  ploidy = 2,
  pop = as.factor(pop))

# run amova 
amova_result <- poppr.amova(genind_obj, ~pop)

amova_result

#### RESULTS 
  #$call
  #ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

  #$results
    #Df   Sum Sq   Mean Sq
    #Between pop                 10 10529.81 1052.9812
    #Between samples Within pop 207 39709.53  191.8335
    #Within samples             218 49117.57  225.3100
    #Total                      435 99356.92  228.4067

  #$componentsofcovariance
    #Sigma          %
    #Variations  Between pop                 24.80055  10.627032
    #Variations  Between samples Within pop -16.73825  -7.172338
    #Variations  Within samples             225.30998  96.545306
    #Total variations                       233.37227 100.000000

  #$statphi
    #Phi
    #Phi-samples-total  0.03454694
    #Phi-samples-pop   -0.08025176
    #Phi-pop-total      0.10627032

#AMOVA Significance test
randtest(amova_result, nrepet = 999)

#RESULTS
  #class: krandtest lightkrandtest 
  #Monte-Carlo tests
  #Call: randtest.amova(xtest = amova_result, nrepet = 999)

  #Number of tests:   3 

  #Adjustment method for multiple comparisons:   none 
  #Permutation number:   999 
  #Test       Obs   Std.Obs   Alter Pvalue
  #1  Variations within samples 225.30998 -1.353488    less  0.094
  #2 Variations between samples -16.73825 -5.138102 greater  1.000
  #3     Variations between pop  24.80055 58.901513 greater  0.001

######### MULTIPANEL CONSOLIDATION 

#FIGURE 1: COUNT + PCA + IBD
install.packages("patchwork")
library (patchwork)

fig1 <- (count_plot | colony | ibd_plot) +
  plot_annotation(tag_levels = "A")

fig1

#FIGURE 2: HEADMAPS
fig2 <- fst_plot | geo_plot +
  plot_annotation(tag_levels = "A")

fig2







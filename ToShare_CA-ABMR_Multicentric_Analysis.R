############### NANOSTRING NORMALIZATION AND ANALYSIS PIPELINE ################

# Place this R script in the analysis folder.
# This folder should contain:
#   a) an annotation file named "clinicallabels.csv" reporting sample names as rownames 
#      and other additional columns for the study group and the batch label
#   b) a file called "BHOT_entrez_mapping.csv" including the BHOT panel entrezIDs list
#   c) a file called "custom_pathways.rds"
#   d) a file called "custom_celltypes.rds"
#   d) a "Normalization_Folder" subfolder containing all RCC files

######## LOAD OR INSTALL NECESSARY PACKAGES ########

# Install if needed the NanoNormIter pipeline
#devtools::install_github("bhattacharya-a-bt/NanoNormIter")

#install.packages("pacman")
library("pacman")
pacman::p_load(MASS,
               tibble,
               gage,
               DESeq2,
               tidyverse,
               Biobase,
               NanoNormIter,
               rstudioapi,
               nanostringr,
               stringr,
               plotly,
               ggfortify,
               dplyr,
               gplots,
               RColorBrewer,
               ggplot2,
               ggrepel,
               umap,
               reticulate) 

######## DATA IMPORT AND PREPROCESSING ########

# Set Working Directory to folder of the script
setwd(dirname(getActiveDocumentContext()$path))

# Import Annotations File
clin_lab <- read.csv("clinicallabels.csv", row.names = 1) # Annotations File

## Reading RCC files
rcc_files <- read_rcc("Normalization_Folder/")
raw <- rcc_files$raw %>% # Import raw counts
  select(-c("Code.Class", "Accession")) %>%
  column_to_rownames(var = "Name")
# Column Name has to match SampleID as in the annotation file loaded later

p_data <- NanoStringQC(rcc_files$raw, rcc_files$exp) %>% # Import Patient Data
  column_to_rownames(var = "File.Name")
# Column Name has to match SampleID as in the annotation file loaded later

f_data <- rcc_files$raw %>% # Import probeset information
  select(c("Name", "Code.Class", "Accession")) %>%
  column_to_rownames(var = "Name")
colnames(f_data)[1] <- "Class"

# Filter annotations file and add information from the annotations file in p_data

clin_lab <- clin_lab[row.names(clin_lab) %in% row.names(p_data),] # Filter only samples that are also present in the dataset
p_data$Group <- clin_lab$Group # Adds study group labels to p_data for HK association test
p_data$batch = clin_lab$Batch # Adds batch labels to p_data for HK association test

#### Verify if HK genes are associated with phenotypes or biological conditions ####

# Create list of HK genes and identify samples with no expression of HK genes
c_idx <- rownames(f_data[f_data$Class == "Housekeeping", ])
p_data$HK_Gene_Miss <- colSums(raw[c_idx, ] == 0)
hk_raw <- raw[c_idx, ] # HK genes raw matrix

# Calculate association between HK and the investigated variable
pval <- vector(length = nrow(hk_raw)) # Create p-values vector for associations
for (i in seq_len(nrow(hk_raw))){
  reg <- glm.nb(as.numeric(hk_raw[i, ]) ~ as.factor(p_data$Group))
  pval[i] <- coef(summary(reg))[2, 4]
}

# Add p-values to the hk table
hk_raw <- add_column(hk_raw,
                     pval,
                     .before = 1)

# Create df of the HK to remove from the HK list prior to normalization
biased_hk <- hk_raw[hk_raw$pval < 0.05, ]

# Eventually label associated HK genes as endogenous
c_idx <- c_idx[!(c_idx %in% row.names(biased_hk))] # Remove from c_idx biased HK genes
f_data[row.names(f_data) %in% row.names(biased_hk) & f_data$Class == "Housekeeping",]$Class <- "Endogenous"

#### Flag samples with low values close to the limit of detection ####

neg_raw <- raw[f_data$Class == "Negative", ] # Creates Negative Controls data frame
lod <- colMeans(neg_raw) - apply(neg_raw, 2, sd) # Calculate LOD
num_endogenous_blod <- colSums(raw[f_data$Class == "Endogenous", ] < lod) # Count endogenous genes per sample below LOD
num_hk_blod <- colSums(raw[f_data$Class == "Housekeeping", ] < lod) # Count housekeeping genes per sample below LOD

# Remove unwanted samples by selecting number of housekeeping genes below lod to dump the sample
raw <- raw[, ! (colnames(raw) %in% names(num_hk_blod[num_hk_blod > 7]))] # Remove from raw counts
p_data <- p_data[!(rownames(p_data) %in% names(num_hk_blod[num_hk_blod > 7])), ] # Remove from p_data
clin_lab <- clin_lab[!(row.names(clin_lab) %in% names(num_hk_blod[num_hk_blod > 7])), ] # Remove from annotations

#### Remove Positive and Negative controls ####
raw <- raw[rownames(f_data[f_data$Class != "Positive" & f_data$Class != "Negative",]), ]
f_data <- f_data[f_data$Class != "Positive" & f_data$Class != "Negative",]

#### Eventually remove samples following post-normalization/visualization ####
# Try to start conservative and drop after visualizing.

######## NORMALIZATION ########

# Iteratively try tuning k, removing flagged samples, and changing the control genes
# Use visualization to achieve optimal results
# Look at expression heatmaps and see if things look odd on the genes correlation level
# Work with LOD flags and HK genes to normalize as well as possible.

k <- 1 # Try with k values ranging from 1 to 5
norm_dat <- RUV_total(raw, p_data, f_data, k = k) # In the iterative process, run this line substituting raw to renormalize

# Visualize relative log expression and color code on known sources of technical variation
# to make decisions about samples or k
col_batch <- vector() # Create color coding vector
col_batch[p_data$batch == "EMC"] <- "red" # Insert Factor name/s
col_batch[p_data$batch == "MGH"] <- "blue"

plotRLE(norm_dat$set,
        col = col_batch,
        style = "full",
        main = "",
        xlab = "Sample",
        xaxt = "n")

# Extract normalized log counts to plot
log_norm_counts <- log(normCounts(norm_dat$set) + 1) # Rerun this line at every renormalization for downstream analysis

# Generate Heatmap of normalized data
pal <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)
heatmapgrouping <- vector(length =length(p_data$Group))
heatmapgrouping[p_data$Group == "positive"] <- "#b502cf" 
heatmapgrouping[p_data$Group == "negative"] <- "#02bf87"
heatmap.2(log_norm_counts,
          scale = "row",
          trace = "none",
          col = pal,
          ColSideColors = heatmapgrouping,
          main = "Heatmap of Normalized mRNA Counts",
          ylab = "Gene",
          xlab = "Sample") # Heatmaps of log normalized data
par(lend = 1)           # square line ends for the color legend
legend(x = -4.33, # location of the legend on the heatmap plot
       y = 3.44,
       xpd = T,
       legend = c("Negative", "Positive"), # category labels
       col = c( "#b502cf" ,"#02bf87"),  # color key
       lty= 1,             # line style
       lwd = 10,            # line width
)

# PCA of non-normalized data for comparison

pca_test_nn <- prcomp(t(log(counts(norm_dat$set) + 1)), scale. = F)
ggplotly(autoplot(pca_test_nn,
                  data = clin_lab,
                  colour = "Batch",
                  label = F))

# PCA of normalized data

pca_test <- prcomp(t(log_norm_counts), scale. = T)
#pca_test[["scale"]] <- FALSE # Optional
ggplotly(autoplot(pca_test,
                  data = clin_lab,
                  colour = "Batch",
                  label = F))
ggplotly(autoplot(pca_test,
                  data = clin_lab,
                  colour = "Group",
                  label = F))

# Extract Set Counts and Variance Stabilizing Data for iterative normalization and DEG

vsd = norm_dat$vsd
set = norm_dat$set

##### Negative Control Selection and Outlier Replacement #####

# Steps 1-3 Should be followed only at the first iteration. Successively, work only on 4 to select control genes.

## Step 1: Create DESeqDataSet object without including RUV W_1 correction variable
dds <- DESeqDataSetFromMatrix(countData = counts(set),
                              colData = pData(set),
                              design = ~ Group)

## Step 2: Run DESeq without continuous variable
dds <- DESeq(dds)

## Step 3: (Optional) Identify and replace outliers. Otherwise skip to 4
dds <- replaceOutliers(dds)
raw <- assay(dds)

## Step 4: Identify empirical control genes. Tweak the lfc and padj thresholds to taste
nondegs <- as.data.frame(results(dds,contrast = c('Group','positive','negative'), altHypothesis = "lessAbs",
                                 lfcThreshold=0.6)) # Performs 1st line DEG analysis with differential expression as null hypothesis to identify genes with equal expression
empirical <- rownames(nondegs[nondegs$padj < 0.05 & !(is.na(nondegs$padj)),])
# Changes the gene class of the newly selected HK genes to HK
c_idx <- unique(c(c_idx, empirical))
f_data[row.names(f_data) %in% empirical & f_data$Class == "Endogenous",]$Class <- "Housekeeping"

## Step 5: Re-normalize and re-visualize results. 
# Try to improve normalization acting on control genes. Once satisfied, proceed to Step 6a. Proceed to step 6b if outliers
# replacement has been performed but no control gene has been changed 

## Step 6a: Rerun normalization without outliers and run this   
# to create new DESeqDataSet object with corrected counts and full design correcting with RUV
dds <- DESeqDataSetFromMatrix(countData = counts(set),
                              colData = pData(set),
                              design = ~ W_1 + Group)

# or

# Step 6b: Run directly DESeq from outliers replaced dataset
dds <- DESeqDataSetFromMatrix(countData = assay(dds),
                              colData = pData(set),
                              design = ~ W_1 + Group)

######## DIFFERENTIAL GENE EXPRESSION ANALYSIS #########

dds <- DESeq(dds) # Run DESeq on normalized data with full design
plotPCA(rlog(dds), intgroup = "Group") # Re-check outliers
resDEG <- as.data.frame(results(dds,contrast = c('Group','positive','negative'))) # DEG Analysis

##### PCA Top 15 Genes #####

# Identify 15 most signifcant DEGs 
topDEGs <- as.data.frame(resDEG) %>%
  arrange(padj)
topDEGs <- topDEGs[1:15,]
exp_topDEGs <- t(log_norm_counts[rownames(log_norm_counts) %in% rownames(topDEGs),]) # Extract Most differentilly expressd genes
PCAtest_top15 <- prcomp(exp_topDEGs, scale.= T) # Perform PCA. Try to change "scale." to see the effect on the plot
ggplotly(autoplot(PCAtest_top15, # Generate the plot
                  data = clin_lab,
                  colour = "Group",
                  label = F,
                  main = "Principal Component Analysis",
                  hjust = 0.5))
ggplotly(autoplot(PCAtest_top15, # Generate the plot
                  data = clin_lab,
                  colour = "Batch",
                  label = F))

##### Pathway analysis #####

### Kegg Pathways

bhot_entrez <- read.csv("BHOT_entrez_mapping.csv",
                        sep = ";") # Import B-HOT Panel Entrez IDs
data("kegg.gs") # Import KEGG pathway definitions
pathdata <- log_norm_counts # Prepare input file by filtering and renaming log counts by Entrez ID
pathdata <- pathdata[rownames(pathdata) %in% bhot_entrez$Original.BHOT.names,] # Retain only matching probes
bhot_entrez <- bhot_entrez[bhot_entrez$Original.BHOT.names %in% rownames(pathdata),] # Retain only matching probes
pathdata <- pathdata[bhot_entrez$Original.BHOT.names,] # Reorder to match
rownames(pathdata) <- bhot_entrez$Entrez.Gene # Reassign 
hn= match(row.names(clin_lab[clin_lab$Group == "negative",]), colnames(pathdata)) # Column Index of reference class
dcis = match(row.names(clin_lab[clin_lab$Group == "positive",]), colnames(pathdata)) # Column Index of test class

pathdata.kegg.p <- gage(pathdata, # Perform Pathway analysis
                        gsets = kegg.gs,
                        ref = hn,
                        samp = dcis,
                        compare = "unpaired")
inc_path <- pathdata.kegg.p$greater # Identify Pathways with increased expression
dec_path <- pathdata.kegg.p$less # Identify Pathways with decreased expression

### Harvard's custom Pathways analysis

# Import custom pathways definition
pathlist <- readRDS("custom_pathways.rds")
pathdata.custom.p <- gage(pathdata, # Perform Pathway analysis
                          gsets = pathlist,
                          ref = hn,
                          samp = dcis,
                          set.size = c(5,500),
                          compare = "unpaired",
                          full.table = T)

custom_inc_path <- pathdata.custom.p$greater # Identify Pathways with increased expression
custom_dec_path <- pathdata.custom.p$less # Identify Pathways with decreased expression

### Cell type analysis

# Import custom cell type definitions
celllist <- readRDS("custom_pathways.rds")
celldata.custom.p <- gage(pathdata, # Perform Pathway analysis
                          gsets = celllist,
                          ref = hn,
                          samp = dcis,
                          set.size = c(3,500),
                          compare = "unpaired")
custom_inc_cell <- celldata.custom.p$greater # Identify Pathways with increased expression
custom_dec_cell <- celldata.custom.p$less # Identify Pathways with decreased expression
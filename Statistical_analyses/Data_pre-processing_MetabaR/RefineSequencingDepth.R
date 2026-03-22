# Script Title: Sequencing Depth Assessment and Sample Filtering for Fish Gut Metabarcoding
#      
# Author: Quentin PETITJEAN
# Date Created: 03/12/2025
# Last Modified: 03/12/2025
# ==============================================================================
# Requirements: 
# - R version 4.2.3
# - Packages: 
#   - BiocManager v1.30.26: For installing and managing Bioconductor packages.
#   - phyloseq v1.52.0: For organizing microbiome data and computing diversity metrics.
#   - metabaR 1.0.0: For handling metabarcoding data, rarefaction, and coverage-based QC.
#   - ggplot2 v4.0.0: For generating diagnostic plots.
#   - remote v2.5.0: for installing R package from github
# ==============================================================================
# Script Overview:
# This script evaluates sequencing depth and removes unreliable fish gut samples from an
# already cleaned metabarcoding dataset using coverage-based criteria:
# 1. Install and load the required R packages (CRAN and Bioconductor).
# 2. Specify the working directory ('savingDir') where input data and objects are stored.
# 3. Convert the metabaR object to a phyloseq object and compute, for each sample:
#    sequencing depth, number of MOTUs, and Shannon/Simpson diversity indices.
# 6. Generate diagnostic plots of sequencing depth versus MOTU richness and alpha diversity
#    (Shannon, Simpson) to visualise depth–diversity relationships.
# 7. Summarise the number and proportion of samples falling below a series of candidate
#    sequencing-depth thresholds (e.g. 500, 1000, 2000, 5000, 10000, 15000 reads).
# 8. Perform Hill-number-based rarefaction with metabaR to obtain coverage curves and
#    plot rarefaction patterns across all samples.
# 9. For each PCR/sample, determine the minimum read count required to reach 95% coverage,
#    and identify samples that never achieve this coverage.
# 10. Derive a common rarefaction depth from samples reaching 95% coverage and flag
#     samples that either fail to reach 95% coverage or fall below this depth.
# 11. Subset the metabarlist to retain only reliable samples ('FishRefined') and compute
#     how many samples and MOTUs are removed by this filtering step.
#
# Note that this script assumes that the metabaR object has already undergone basic
# cleaning; it focuses specifically on sequencing-depth/coverage diagnostics and
# coverage-based sample filtering.
# ==============================================================================

##############################################
#       	Install needed packages            #
##############################################
remotes::
if(!require(remotes)){
  install.packages("remotes")
}
if(!require(BiocManager)){
  BiocManager::install("BiocManager")
}
if(!require(phyloseq)){
  BiocManager::install("phyloseq")
}
if(!require(metabaR)){
  remotes::install_github("metabaRfactory/metabaR")
}
if(!require(ggplot2)){
  install.packages("ggplot2")
}
if (!require(plyr)) {
  install.packages("plyr")
}

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported and saved
savingDir <- "C:/Users/qpetitjean/Desktop/WORK/POSTDOC_INP_GOLFECH_2023/R scripts/GutMicrobiomeCaging/ZenodoRepo"

#######################
# Import the data     #
#######################

# import the full dataset from the experiment (treatment, sex and other fish variable)
#FullDat <-
#  read.csv2(file.path(savingDir, "Data/DesignData", "DesignData.csv"), dec = ".", sep = ";")

# import the cleaned dataset (metabaR object) and select only gut samples
Microbiome <- readRDS(file.path(savingDir, "Data/CleanedData", "fguts_Bact_agg_Caging_MergedRepSum.RDS"))
Fish <- Microbiome

###########################
# Check sequencing depth  #
###########################

# transform data into phyloseq object
read.phylo <-
  phyloseq::otu_table(Fish[["reads"]], taxa_are_rows = F) #makes otu table for phyloseq
sample.phylo <- phyloseq::sample_data(Fish[["samples"]])
ps  <- phyloseq::phyloseq(read.phylo, sample.phylo)

# Compute sequencing depth
seq_depth <- phyloseq::sample_sums(ps)

# Compute the number of MOTUs per sample
n_otus <- apply(phyloseq::otu_table(ps), 1, function(x) sum(x > 0))

# Compute Shannon diversity indices
alphaDiv <- phyloseq::estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Group into a data frame
df <- data.frame(
  SampleID = names(seq_depth),
  SequencingDepth = seq_depth,
  NumOTUs = n_otus,
  Shannon = alphaDiv$Shannon,
  Simpson = alphaDiv$Simpson
)

# Plot 1: MOTUs vs Sequencing Depth
ggplot2::ggplot(df, ggplot2::aes(x = SequencingDepth, y = NumOTUs)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "loess") +
  ggplot2::labs(title = "MOTUs vs Sequencing Depth", x = "Sequencing Depth (reads)", y = "Number of MOTUs")

# Plot 2: Shannon vs Sequencing Depth
ggplot2::ggplot(df, ggplot2::aes(x = SequencingDepth, y = Shannon)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "loess") +
  ggplot2::labs(title = "Shannon Diversity vs Sequencing Depth", x = "Sequencing Depth (reads)", y = "Shannon Diversity Index")

# Plot 3: Simpson vs Sequencing Depth
ggplot2::ggplot(df, ggplot2::aes(x = SequencingDepth, y = Simpson)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "loess") +
  ggplot2::labs(title = "Shannon Diversity vs Sequencing Depth", x = "Sequencing Depth (reads)", y = "Shannon Diversity Index")

# check the number of sample to remove according to several rarefying treshold
thresholds <- c(500, 1000, 2000, 5000, 10000, 15000)

# Initialize vectors
below_n <- numeric(length(thresholds))
below_percent <- numeric(length(thresholds))

# Loop to fill values
for (i in seq_along(thresholds)) {
  count <- sum(df$SequencingDepth < thresholds[i])
  below_n[i] <- count
  below_percent[i] <- round(100 * count / nrow(df), 1)
}

threshold_summary <- data.frame(
  Threshold = paste0("<", thresholds, " reads"),
  Num_Samples = below_n,
  Percent_Removed = below_percent
)

# rarefaction curve
Fish.raref <- metabaR::hill_rarefaction(Fish, nboot = 20, nsteps = 10)
material <- paste(Fish$samples$experiment, Fish$samples$matrix)
material <- setNames(material,rownames(Fish$samples))
options(max.print=99999)
head(Fish.raref$hill_table)

#plot
p <- metabaR::gghill_rarefaction(Fish.raref, group=material)
p

# which sequencing depth is needed to reach 95% coverage
hill_table <- Fish.raref$hill_table
pcr_ids <- unique(hill_table$pcr_id)

# Initialize result container
reads_for_95cov <- data.frame(pcr_id = integer(),
                              min_reads_95 = numeric(),
                              stringsAsFactors = FALSE)

# Loop over each pcr_id
for (id in pcr_ids) {
  sub <- hill_table[hill_table$pcr_id == id, ]
  sub_95 <- sub[sub$coverage >= 0.95, ]
  
  if (nrow(sub_95) > 0) {
    min_row <- sub_95[which.min(sub_95$reads), ]
    reads_for_95cov <- rbind(reads_for_95cov,
                             data.frame(pcr_id = id,
                                        min_reads_95 = min_row$reads))
  } else {
    # No 95% coverage reached
    reads_for_95cov <- rbind(reads_for_95cov,
                             data.frame(pcr_id = id,
                                        min_reads_95 = NA))
  }
}

median(reads_for_95cov$min_reads_95, na.rm = T)

# identify samples that do not reach 95% coverage
unreliable <- reads_for_95cov[which(is.na(reads_for_95cov$min_reads_95)), "pcr_id"]
hill_table[which(hill_table$pcr_id %in% unreliable),]

filteredDf <- df[which(!df$SampleID %in% unreliable),]
min(filteredDf$SequencingDepth)
hist(filteredDf$SequencingDepth, breaks=500, xlim=c(0,5000))


#######################################################################
# Identify samples with sequencing depth below rarefy_depth treshold #
######################################################################

# Merge the observed depth and the min_reads needed for 95% coverage
merged <- merge(df, reads_for_95cov, by.x = "SampleID", by.y = "pcr_id", all.x = TRUE)

# Identify valid samples (reaching 95% coverage)
valid <- merged[!is.na(merged$min_reads_95), ]

# Determine the rarefaction depth across those valid samples (for 95% coverage = 1000 reads)
rarefy_depth <- median(valid$min_reads_95)

# Print the number and preview
to_remove <- merged[is.na(merged$min_reads_95) | merged$SequencingDepth < rarefy_depth, ]
samples_to_remove <- to_remove$SampleID

#########################################
# Remove samples that are not reliable #
########################################
# Remove unreliable samples
FishRefined <- 
  metabaR::subset_metabarlist(Microbiome, 
                              "samples", 
                              indices = !rownames(Microbiome$samples) %in% samples_to_remove)

# how many samples removed - 88
metabaR::summary_metabarlist(labMicrobiome)[[1]]["samples",1] - 
  metabaR::summary_metabarlist(FishRefined)[[1]]["samples",1]

# how many Motus removed - 0
metabaR::summary_metabarlist(labMicrobiome)[[1]]["motus",1] - 
  metabaR::summary_metabarlist(FishRefined)[[1]]["motus",1]

# save the dataset for further analyses
saveRDS(FishRefined, file.path(savingDir, "Data/CleanedData", "fguts_Bact_agg_Caging_MergedRepSumRefined.RDS"))

# check the final content of the dataset for gut samples (reported in material and methods in the 2.4.2	Data Pre-Processing section)
FishFinal <- metabaR::subset_metabarlist(FishRefined,
                                       table = "pcrs",
                                       indices = FishRefined$pcrs$matrix == "fish_gut")
metabaR::summary_metabarlist(FishFinal)

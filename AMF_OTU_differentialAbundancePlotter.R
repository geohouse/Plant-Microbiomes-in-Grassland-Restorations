# This script re-creates the basic graph shown in Figure 4 (differential abundance of each OTU
# in either remnant prairies or nearby post-agricultural fields). The graph plotted by 
# this script is a horizontal version of the bargraph presented vertically in Figure 4
# (labels [and alternating light/dark bars] from left-to-right on the horizontal graph 
# from this script are in the same order as they are from top-to-bottom on the vertical
# graph in Figure 4). The phylogeny was manually added in Figure 4 following the accepted
# AM fungal genus (or family) phylogenetic relationships (Redecker et al 2013 Mycorrhiza
# DOI 10.1007/s00572-013-0486-y).

# Input data is from three sources:
#   1) OTU table from which the differential abundances are calculated
#   2) Metadata about whether each sample was collected from a remnant or a disturbed site
#   3) Taxonomic attribution of each OTU

library('phyloseq')
library('DESeq2')

#-----------------------------
# Import the OTU table
OTUTable <- read.table("~/Box Sync/R_code/Gradient_dataAndPlotter_forLizs_review/AMF_OTU_tableForPlotting.tsv", header=TRUE, sep="\t", stringsAsFactors = FALSE)

# metadata for each of the samples, particularly whether they are from remnant or disturbed sites.
metadataForSamples <- read.table("~/Box Sync/R_code/Gradient_dataAndPlotter_forLizs_review/metadataForSamples.tsv", sep = "\t", header = TRUE)

# taxonomic attribution for each OTU
taxonomicAttributions <- read.table("~/Box Sync/R_code/Gradient_dataAndPlotter_forLizs_review/taxonomicAttributionOfAMF_OTUs.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Need to add a leading 'X' to the taxonomic attributions before using them to match
# against the column names in the OTUTable
OTUNamesMatchOTUTable <- paste("X", taxonomicAttributions$OTU_num, sep = "")

# Order the OTU name indices in OTUTable the same as they appear in the taxonomicAttributions
OTUTableIndicesInPhylogOrderForLabels <- match(OTUNamesMatchOTUTable, colnames(OTUTable))

#OTUTableIndicesInPhylogOrderForLabels <- na.omit(OTUTableIndicesInPhylogOrderForLabels)

# Add pseudo count - this is REQUIRED  to be able to run the rlog transformation below and the DESeq analysis.
OTUTable <- OTUTable + 1

# Convert OTU table to phyloseq OTU table
phyloseq_OTUTable <- otu_table(OTUTable, taxa_are_rows = FALSE)

# This seems to work, but View() fails on it....
phyloseq_sampleData <- sample_data(metadataForSamples)

phyloseq_combined <- phyloseq(phyloseq_OTUTable, phyloseq_sampleData)

OTU_sampleDataForDESeq <- phyloseq_to_deseq2(phyloseq_combined, ~ Disturbed_Remnant)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in remnant; - val = higher in disturbed sites), along with p-values for the log2 fold change.
DESeqOutput <- DESeq(OTU_sampleDataForDESeq, fitType = 'local')

DESeqOutput_results_DistHistCompare <- results(DESeqOutput)
DESeqOutput_resultsPValSorted_DistHistCompare <- DESeqOutput_results_DistHistCompare[order(DESeqOutput_results_DistHistCompare$padj),]
DESeqOutput_resultsEffectSizeSorted_DistHistCompare <- DESeqOutput_results_DistHistCompare[order(abs(DESeqOutput_results_DistHistCompare$log2FoldChange), decreasing = TRUE),]
DESeqOutput_resultsLabelSorted_DistHistCompare <- DESeqOutput_results_DistHistCompare[OTUTableIndicesInPhylogOrderForLabels,]

# Make the plot.  
par(mar=c(4.1,6.1,3.1,1.1))
# Ordered by phylogeny
# Leave width at 0.835 for the segments to fit well.
phyloAdjPValColor <- ifelse(DESeqOutput_resultsLabelSorted_DistHistCompare$padj < 0.05, 'green', 'white')
barplot(DESeqOutput_resultsLabelSorted_DistHistCompare$log2FoldChange, col = phyloAdjPValColor,
        ylab = 'Log 2 fold diff Remn (+) /Dist (-)',
        ylim = c(-8,8), width = 0.835, cex.lab = 1.5, cex.axis = 1.5, yaxt = "n")
# Manually add the axis to get the ticks at the correct spacing.
axis(2, at = seq(-8,8,2), las = 1)
    
# This is to add the genus/demarking color line and dashed lines 
lineThickness <- 3
abLineColor <- "gray"
abLineThickness <- 0.5
segmentBarWidth <- 1
# Keep at 0.3 and save as 16"W by 4" H
barOffset <- 0.3

contrastBarColor1 <- "black"
contrastBarColor2 <- "gray"

yLocation = -7.9

# For Div. (9)
segments(x0 = barOffset, 
         y0 = yLocation, 
         x1 = barOffset + 9 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 9 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Aca. (15)
segments(x0 = barOffset + 9 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 24 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 24 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Scu. (2)
segments(x0 = barOffset + 24 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 26 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 26 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Cet. (1)
segments(x0 = barOffset + 26 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 27 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 27 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Gig. (2)
segments(x0 = barOffset + 27 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 29 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 29 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Cla. (13)
segments(x0 = barOffset + 29 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 42 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 42 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Glomeraceae. (52)
segments(x0 = barOffset + 42 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 94 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 94 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Rhi. (52)
segments(x0 = barOffset + 94 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 146 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 146 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Glomus. (17)
segments(x0 = barOffset +146 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 163 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 163 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Fun. (3)
segments(x0 = barOffset + 163 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 166 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 166 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Sep. (9)
segments(x0 = barOffset + 166 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 175 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 175 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Amb. (2)
segments(x0 = barOffset + 175 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 177 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 177 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Arc. (1)
segments(x0 = barOffset + 177 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 178 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor1, lwd = lineThickness, lend = "butt")
abline(v = barOffset + 178 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)

# For Par. (3)
segments(x0 = barOffset + 178 * segmentBarWidth, 
         y0 = yLocation, 
         x1 = barOffset + 181 * segmentBarWidth, 
         y1 = yLocation, 
         col = contrastBarColor2, lwd = lineThickness, lend = "butt")


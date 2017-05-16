### Adapted from
## RNA-seq analysis with DESeq2
## Stephen Turner, @genetics_blog
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

## These are the required packages. Uncomment these lines if they are not already installed -------
#install.packages(DESeq2)
#install.packages(ggplot2)
#install.packages(gplots)
#install.packages(RColorBrewer)
#install.packages(gridExtra)
#install.packages(genefilter)
#install.packages(calibrate)

## Dataset information ----------------------------------------------------------------------------
print("Setting the dataset information ...")
# What is the name of the project that should be included in each output file
project_name = "PRJNA344306"
# The directory where the output files will be placed
output_dir <- paste("/media/ab/RData/PRJNA344306/R_Output_", project_name) 
# The highest directory level common to all input datasets
file_dir <- "/media/ab/RData/PRJNA344306"
# The name of each input file, including the subdirectory location relative to the file directory. 
# The input files should have at least 2 columns: The gene IDs and the counts.
# The control datasets must be listed first, then the experimental datasets. 
dataset_list <- c("SRR4292750/SRR4292750.featureCounts.txt", 
                  "SRR4292751/SRR4292751.featureCounts.txt", 
                  "SRR4292752/SRR4292752.featureCounts.txt", 
                  "SRR4292753/SRR4292753.featureCounts.txt"
                  )

# The number of control and experimental datasets
num_control_datasets = 2
num_experimental_datasets = 2
# The column the count data located is in, from original file. 
count_column <- 7

## Setup the output directory ---------------------------------------------------------------------
print("Setting the output directory ...")
# Create the output directory
new_dir <- ifelse(test = !dir.exists(output_dir), yes = dir.create(output_dir, showWarnings = FALSE), no = FALSE)
# If the directory already exists, ask if it is ok to overwrite any existing files
if (new_dir == FALSE) {
  answer <- readline(prompt="The output folder already exists. Do you wish to continue and overwrite any existing contents? (Y or N)")
  if ((answer != "y") & (answer != "Y") & (answer != "yes") & (answer != "Yes") & (answer != "YES")) {
    stop()    
  } 
}
# Set the working directory
setwd(output_dir)

## Import & pre-process the count data ------------------------------------------------------------
print("Importing the count datasets ...")
first_dataset <- TRUE
# Load the datasets and combine into a single dataframe
for (val in dataset_list) {
  if (first_dataset == TRUE) {
    print("Adding the first dataset ...")
    countdata <- read.table(paste(file_dir, val, sep = "/"), header=TRUE, row.names=1)[count_column-1]
    first_dataset <- FALSE
  } else {
    print("Adding another dataset ...")
    tmp_countdata <- read.table(paste(file_dir, val, sep = "/"), header=TRUE, row.names=1)[count_column-1]
    countdata <- merge(countdata, tmp_countdata, by='row.names')
    countdata <- data.frame(countdata[,-1], row.names=countdata[,1])
  }
}

# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <- sapply(colnames(countdata), function(colname) 
  tail(unlist(strsplit(x = colname, split = "\\.", fixed=FALSE)), 1))

# Now that the data is all numerical, convert the dataframe to a matrix
countdata <- as.matrix(countdata)

# Assign the conditions (first four are controls, second four contain the expansion)
condition <- factor(c(rep("ctl", num_control_datasets), rep("exp", num_experimental_datasets)))

## Start the analysis with DESeq2 -----------------------------------------------------------------
print("Loading the necessary libraries ...")
# Load the necessary libraries
library(DESeq2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(gridExtra)


## Create the main DESeq2 dataset -----------------------------------------------------------------
print("Running the main DESeq2 function ...")
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# Run the main DESeq function
dds <- DESeq(dds)

# Colors for plots
# (mycols <- 1:length(unique(condition)))
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]

## Plot the dispersions  --------------------------------------------------------------------------
print("Creating the dispersion plot ...")
# Show the level of dispersion for each gene
png(paste(project_name, "_qc-dispersions.png"), 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)

# Turn the log transformed data into a data frame, and replace the columns with the correct columns
df <- data.frame((assay(rld)))
colnames(df) <- colnames(assay(rld))

## Create histograms of log transformed expression data -------------------------------------------
print("Creating the histograms of log2-transformed expression data ...")
# Initilize the png for the histograms
png(paste(project_name, "_log-transformed.png"), 750, 750, pointsize=40)
# Function that creates a histogram for every dataset
plot_data_column = function (df, i) 
  qplot(df[i], geom = "histogram",
        binwidth = 0.5,  
        main = paste("Histogram for: ", i),  
        ylab = "Count",
        xlab = "Log transformation of values",
        fill="blue") + ggplot2::scale_x_continuous() + theme(legend.position="none")
# Combine the histograms into 1 file
myplots <- lapply(X = colnames(df), FUN = plot_data_column, df = df)
# Show the histograms
n <- length(myplots) 
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(myplots, ncol=nCol))
dev.off()

## Heatmaps of the distances between samples ------------------------------------------------------
print("Creating the heatmap of sample distances ...")
sampleDists <- as.matrix(dist(t(assay(rld))))
png(paste(project_name, "_qc-heatmap-samples.png"), w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

## Principal components analysis to group datasets ------------------------------------------------
print("Creating the PCA plot of samples ...")
# Could do with built-in DESeq2 function:
# DESeq2::plotPCA(rld, intgroup="condition")
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  #require(calibrate)
  #require(MASS)
  #require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("blue", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png(paste(project_name, "_qc-pca.png"), 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()

## Save the differential expression results to a csv file -----------------------------------------
print("Saving the differential experssion results to a CSV file ...")
res <- results(dds, alpha=0.05, independentFiltering = TRUE)
table(res$padj<0.05)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
# Write results
write.csv(resdata, file=paste(project_name, "_diffexpr-results.csv"))

## Histogram of p-values --------------------------------------------------------------------------
print("Creating the histogram of p-values ...")
png(paste(project_name, "_p-values.png"), 500, 500, pointsize=80)
qplot(x = res$pvalue[!is.na(res$pvalue)], 
      geom = "histogram", 
      bins = 100,
      xlim = c(0,1),
      main = paste("Histogram of p-values for ",  project_name),  
      ylab = "Count",
      xlab = "p values" 
      ) + theme(legend.position="none")
dev.off()

## Make plot of theta rejections for independent fiiltering ---------------------------------------
print("Creating the plot of genes rejected by independent filtering ...")
number_of_rejections <- data.frame(attr(res, "metadata")["filterNumRej"])
colnames(number_of_rejections) <- c("x", "y")
png(paste(project_name, "_theta-rejection.png"), 750, 750, pointsize=40)
qplot(x = number_of_rejections$x, y=number_of_rejections$y,
      geom = "point", 
      main = paste("Rejection at theta: ", attr(res, "metadata")["filterTheta"]), 
      xlab="quantiles of baseMean", 
      ylab="number of rejections",
      xlim = c(0,1))
dev.off()

## MA plot ----------------------------------------------------------------------------------------
print("Creating the MA plot ...")
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  all_datapoints <- subset(res, baseMean>0 & !is.na(baseMean))
  all_datapoints$sig <- ifelse(all_datapoints$padj<thresh, 1, 0)
  qplot(data=all_datapoints, x=baseMean, y=log2FoldChange, 
        geom = "point", log = "x", shape=20, 
        xlab = "Normalized counts per gene",
        ylab = "log2 fold change in expression between conditions") +
    geom_point(aes(color = cut(all_datapoints$sig, c(-Inf, 0.5, Inf)))) +
    scale_color_manual(name = "all_datapoints$sig", 
                       values = c("(-Inf,0.5]" = mycols[1], "(0.5, Inf]" = mycols[2]),
                       labels = c("NS", "Sig")) +
    geom_hline(yintercept = 1, linetype="dashed", size=1.5) +
    geom_hline(yintercept = -1, linetype="dashed", size=1.5) +
    theme_set(theme_grey(base_size = 30)) + 
    theme(legend.position="none") +
    geom_text(aes(label=ifelse(abs(all_datapoints$log2FoldChange)>=1 & labelsig,all_datapoints$Gene,''), 
                  hjust=0,vjust=0), color = "blue", show.legend = NA) +
    scale_shape_identity()
}
png(paste(project_name, "_diffexpr-maplot.png"), 1000, 1000, pointsize=40)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled  -------------------------------------------------
print("Creating the volcano plot ...")
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png(paste(project_name, "_diffexpr-volcanoplot.png"), 1200, 1000, pointsize=50)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()
print("DESeq2 analysis is complete")
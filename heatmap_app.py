############################################################
## Quick Heatmap Generator for RNA-seq expression data
## Author: Prekovic Lab
## Usage:
##   source("custom_heatmap.R")
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(readxl)
  library(pheatmap)
  library(viridis)
  library(matrixStats)
})

## --- 1) File input ----------------------------------------------------------
message("📂 Reading input files...")

# ---- User file paths (adjust or use file.choose())
expr_file <- file.choose()      # e.g. counts_matrix_gene_with_symbols.txt
annot_file <- file.choose()     # e.g. annotations_samples.xlsx

# ---- Read files
count_data <- read.table(expr_file, header=TRUE, sep="\t", check.names=FALSE)
rownames(count_data) <- count_data$Geneid
counts <- count_data[, 8:ncol(count_data)]

annotation <- as.data.frame(read_excel(annot_file))
annotation$SampleName <- annotation$...1
annotation <- annotation[,-1]
rownames(annotation) <- annotation$SampleName

# ---- Match order of samples
counts <- counts[, colnames(counts) %in% annotation$SampleName, drop=FALSE]
annotation <- annotation[colnames(counts), , drop=FALSE]

# ---- Normalize (VST)
sym <- ifelse(nchar(count_data$Gene_Symbol)>0, count_data$Gene_Symbol, count_data$Geneid)
sym <- make.unique(sym)
rownames(counts) <- sym[rownames(counts)]
dds <- DESeqDataSetFromMatrix(as.matrix(counts), annotation, design=~1)
vsd <- vst(dds, blind=TRUE)
expr <- assay(vsd)
coldata <- as.data.frame(colData(dds))
rownames(coldata) <- annotation$SampleName

## --- 2) User input for heatmap ---------------------------------------------
cat("\n✅ Data loaded. Available groups:\n")
print(unique(annotation$group))
cat("\nEnter groups to include, separated by commas (e.g. VEH,TACT,TACT_BAF): ")
grp_input <- scan(what="character", sep=",")
grp_input <- trimws(grp_input)
samples <- annotation$SampleName[annotation$group %in% grp_input]

cat("\nEnter genes of interest separated by commas (e.g. FOXP3,CTLA4,PDCD1): ")
gene_input <- scan(what="character", sep=",")
gene_input <- trimws(gene_input)

## --- 3) Prepare matrix ------------------------------------------------------
pal <- viridis::magma(100)
sel_genes <- intersect(gene_input, rownames(expr))
if (length(sel_genes) == 0) stop("No matching gene symbols found.")

mat <- expr[sel_genes, samples, drop=FALSE]
ann <- data.frame(group=annotation[samples, "group", drop=TRUE])
rownames(ann) <- samples

## --- 4) Plot heatmap --------------------------------------------------------
message("\n🎨 Generating heatmap...")

pdf("custom_heatmap.pdf", width=8, height=6, useDingbats=FALSE)
pheatmap(mat,
         scale="row",
         color=pal,
         show_rownames=TRUE,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         annotation_col=ann,
         border_color=NA,
         fontsize=9,
         main=paste("Custom Heatmap –", paste(sel_genes, collapse=", ")))
dev.off()

message("\n✅ Heatmap saved as 'custom_heatmap.pdf' in current directory.")

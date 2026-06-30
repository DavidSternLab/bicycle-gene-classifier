#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rtracklayer)
  library(optparse)
})

# Command line argument parsing
option_list = list(
  make_option(c("-g", "--gff"), type="character", default=NULL,
              help="Input GFF3 file with CDS annotations", metavar="character"),
  make_option(c("-m", "--model"), type="character",
              default=Sys.getenv("BICYCLE_MODEL", unset = ""),
              help="Path to trained GLM model file [default from BICYCLE_MODEL env var]", metavar="character"),
  make_option(c("-c", "--cutoff"), type="numeric", default=0.72,
              help="Classification cutoff threshold [default= %default]", metavar="number"),
  make_option(c("-o", "--output"), type="character", default="bicycle_output",
              help="Output prefix [default= %default]", metavar="character"),
  make_option(c("-d", "--outdir"), type="character", default="bicycle_results",
              help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check required arguments
if (is.null(opt$gff) || is.null(opt$model) || opt$model == ""){
  print_help(opt_parser)
  stop("Both --gff and --model arguments are required.\nFor --model, either use the flag or set BICYCLE_MODEL environment variable.", call.=FALSE)
}

# Check if files exist
if (!file.exists(opt$gff)) {
  stop(paste("GFF file not found:", opt$gff), call.=FALSE)
}
if (!file.exists(opt$model)) {
  stop(paste("Model file not found:", opt$model,
             "\nUse `bicycle_classifier --download-model` to fetch the default Hcor model,",
             "\nor set BICYCLE_MODEL to point at a local copy."), call.=FALSE)
}

# Numeric / output validation
if (is.na(opt$cutoff) || opt$cutoff < 0 || opt$cutoff > 1) {
  stop(paste("--cutoff must be a number in [0, 1]; got:", opt$cutoff), call.=FALSE)
}
if (nchar(opt$output) == 0) {
  stop("--output prefix cannot be empty", call.=FALSE)
}
# Create outdir early so we fail fast if it's unwritable.
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(opt$outdir)) {
  stop(paste("Could not create output directory:", opt$outdir), call.=FALSE)
}
if (file.access(opt$outdir, mode = 2) != 0) {
  stop(paste("Output directory not writable:", opt$outdir), call.=FALSE)
}

cat("Starting bicycle gene classification...\n")
cat("GFF file:", opt$gff, "\n")
cat("Model file:", opt$model, "\n")
cat("Cutoff:", opt$cutoff, "\n")
cat("Output prefix:", opt$output, "\n")
cat("Output directory:", opt$outdir, "\n\n")

# Load the trained model
load(opt$model)
if (!exists("glm.full")) {
  stop("Model file does not contain 'glm.full' object", call.=FALSE)
}

# Main prediction function
predict_bicycle = function(gff.cds, cutoff, output_prefix, output_folder_name){
  
  cat("Processing annotation file...\n")
  
  #format annotation file
  gff.cds$Parent = unlist(gff.cds$Parent)
  parents = unique(gff.cds$Parent)
  gff.cds$Parent =factor(gff.cds$Parent,levels = parents)
  gff.cds = data.frame(gff.cds)
    
  #general genes statistics
  exons.summary = gff.cds %>%
    mutate(size = end - start + 1) %>% group_by(Parent) %>%
    dplyr::summarise(num_total_exons = n(), total_exon_length = sum(size), .groups = 'drop')
  
  #length of gene from start of first exon to end of last exon
  gene.length.summary = gff.cds %>%
    group_by(Parent) %>%
    dplyr::summarise(gene_start = dplyr::first(start), gene_end = dplyr::last(end), 
                    gene_length = abs(gene_end-gene_start)+1, .groups = 'drop') %>%
    dplyr::select(Parent, gene_length)
    
  pos_strand = gff.cds %>% filter(strand == "+")
  neg_strand = gff.cds %>% filter(strand == "-")
  
  pos.first.exons.summary = pos_strand %>% group_by(Parent) %>% filter(row_number() == 1) %>% 
    mutate(size = end - start + 1) %>%
    dplyr::summarise(first_exon_length = size, .groups = 'drop')
  pos.last.exons.summary = pos_strand %>% group_by(Parent) %>% 
    filter(row_number() != 1 & row_number() == n()) %>% 
    mutate(size = end - start + 1) %>%
    dplyr::summarise(last_exon_length = size, .groups = 'drop')
    
  neg.first.exons.summary = neg_strand %>% group_by(Parent) %>% 
    filter(row_number() == n()) %>% 
    mutate(size = end - start + 1) %>%
    dplyr::summarise(first_exon_length = size, .groups = 'drop')
  neg.last.exons.summary = neg_strand %>% group_by(Parent) %>% 
    filter(row_number() == 1 & row_number() != n()) %>% 
    mutate(size = end - start + 1) %>%
    dplyr::summarise(last_exon_length = size, .groups = 'drop')
  
  first.exons.summary = rbind(pos.first.exons.summary, neg.first.exons.summary)
  last.exons.summary = rbind(pos.last.exons.summary, neg.last.exons.summary)
    
  #internal exons statistics
  internal.exons = gff.cds %>% group_by(Parent) %>% filter(row_number() != 1 & row_number() != n())
  
  internal.exons.summary = internal.exons %>%
    mutate(size = end - start + 1) %>%
    group_by(Parent) %>%
    dplyr::summarise(num_internal_exons = n(), exon_mean_length = mean(size), exon_var = var(size), 
              mode0=sum(phase==0), mode1=sum(phase==1), 
              mode2=sum(phase==2), .groups = 'drop')
  
  #put all above tables together 
  gene.summary = merge(merge(merge(merge(exons.summary, gene.length.summary, by = "Parent", all = TRUE), 
                                   first.exons.summary, by = "Parent", all = TRUE), last.exons.summary, by = "Parent", all = TRUE),internal.exons.summary, by = "Parent", all = TRUE)
  
  cat("Before filtering:", nrow(gene.summary), "genes\n")
  
  #remove all genes with 3 exons or less 
  gene.summary = gene.summary[complete.cases(gene.summary), ]
  
  cat("After filtering (complete cases only):", nrow(gene.summary), "genes\n")
  
  if(nrow(gene.summary) == 0) {
    stop("No genes with complete feature data found. Check that your annotation has:\n", 
         "  - Genes with >3 exons\n",
         "  - Phase information\n", 
         "  - Proper Parent/transcript_id attributes", call.=FALSE)
  }
  
  cat("Extracted features for", nrow(gene.summary), "genes\n")
  cat("Running classifier...\n")
  
  #run classifier to predict bicycle genes
  full.predict = predict(glm.full, gene.summary, type="response")
  full.predict.merged = cbind.data.frame(gene.summary$Parent, full.predict)
  colnames(full.predict.merged) = c("Parent", "response")

  #output classifier response for all transcripts
  dir.create(output_folder_name, showWarnings = FALSE)
  write.table(full.predict.merged, paste0(output_folder_name, "/", output_prefix, "_classifier_all_transcripts_response.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')
    
  #filter for transcripts above the cutoff threshold
  bicycle.glm = full.predict.merged[full.predict.merged$response > cutoff,]
  bicycle.gene.names = unique(gsub("\\.[^.]*$","", bicycle.glm$Parent))
  
  cat("Found", length(bicycle.gene.names), "potential bicycle genes above cutoff\n")
  
  #output gene names for genes above the cutoff threshold
  write.table(bicycle.gene.names, paste0(output_folder_name, "/", output_prefix, "_classifier_bicycle_gene_names.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
  
  #plot histogram distribution of classifier response for all transcripts
  pdf(paste0(output_folder_name, "/", output_prefix, "_classifier_response_histogram.pdf"), width=4, height=4)
  print(ggplot() + geom_histogram(aes(full.predict.merged$response), bins = 100) + scale_y_log10() +
  geom_vline(aes(xintercept=cutoff, color="red")) +
  theme_classic() + theme(legend.position = "none", 
        axis.text=element_text(size=16),
        axis.title=element_text(size=18)) + 
  xlab("response") + ylab("transcript count"))
  dev.off()
  
  cat("Results written to:", output_folder_name, "\n")
  cat("Classification complete!\n")
}

# Read input files and run analysis
cat("Reading annotation file...\n")

# Try to read as GFF/GTF and handle different formats
gff.cds = tryCatch({
  readGFF(opt$gff)
}, error = function(e) {
  stop(paste("Error reading annotation file:", e$message), call.=FALSE)
})

cat("Total features read:", nrow(gff.cds), "\n")

if(nrow(gff.cds) == 0) {
  stop("Annotation file appears to be empty or unreadable", call.=FALSE)
}

# Show file structure for debugging
cat("Columns in file:", paste(colnames(gff.cds), collapse=", "), "\n")
cat("First few feature types:", paste(head(unique(gff.cds$type), 10), collapse=", "), "\n")

# Check what feature types are available
feature_types = unique(gff.cds$type)
cat("All available feature types:", paste(feature_types, collapse=", "), "\n")

# Filter for CDS features
if("CDS" %in% feature_types) {
  gff.cds = gff.cds[gff.cds$type == "CDS", ]
  cat("Using CDS features\n")
} else if("exon" %in% feature_types) {
  gff.cds = gff.cds[gff.cds$type == "exon", ]
  cat("Using exon features (no CDS found)\n")
} else {
  cat("Available feature types in detail:\n")
  print(table(gff.cds$type))
  stop("No CDS or exon features found in annotation file", call.=FALSE)
}

cat("Found", nrow(gff.cds), "CDS/exon features\n")

# Handle different attribute naming (GTF vs GFF3).
# readGFF() returns Parent as a CharacterList (one list-element per row, each a
# vector of parent IDs). Flatten it to a plain character vector up front so the
# subsequent `||` short-circuit doesn't get a length>1 vector (which is a
# warning in R 4.2 and an error in R 4.3+).
flat_parent <- if (!is.null(gff.cds$Parent)) as.character(unlist(gff.cds$Parent)) else NULL
parent_missing <- is.null(flat_parent) || length(flat_parent) == 0 || all(is.na(flat_parent))
if (parent_missing) {
  if(!is.null(gff.cds$transcript_id)) {
    cat("GTF format detected, using transcript_id as Parent\n")
    gff.cds$Parent = gff.cds$transcript_id
  } else if(!is.null(gff.cds$gene_id)) {
    cat("Using gene_id as Parent\n")
    gff.cds$Parent = gff.cds$gene_id
  } else {
    cat("Available attribute columns:", paste(colnames(gff.cds), collapse=", "), "\n")
    stop("Cannot find Parent, transcript_id, or gene_id attributes", call.=FALSE)
  }
} else {
  # Replace the CharacterList with its flat form so downstream code (and the
  # `length(unique(...))` count below) get the right per-CDS Parent value.
  gff.cds$Parent <- flat_parent
}

# Check that we have phase information (required for the classifier)
if(is.null(gff.cds$phase) || all(is.na(gff.cds$phase))) {
  cat("Warning: No phase information found. Setting all phases to 0.\n")
  gff.cds$phase = 0
}

cat("Processed", length(unique(gff.cds$Parent)), "unique transcripts/genes\n")

# Run prediction function
predict_bicycle(gff.cds, opt$cutoff, opt$output, opt$outdir)

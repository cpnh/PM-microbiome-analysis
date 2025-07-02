#!/usr/bin/env Rscript

#######################################
#                                     #
#       Microbiome Rarefaction       #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                 #
# GitHub: promicrobial (https://github.com/promicrobial)                     #
# Date: 04-02-25                                                             #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: This script performs rarefaction of microbiome data,          #
#              separating and processing Bacterial and Archaeal sequences    #
#              independently. It includes data import, quality control,      #
#              rarefaction, and saving processed data.                       #
#                                                                            #
# Dependencies:                                                              #
#   - phyloseq: For microbiome data handling                                 #
#   - tidyverse: For data manipulation                                       #
#   - microbiome: For phyloseq object summaries                              #
#   - optparse: For command line argument parsing                            #
#   - ape: For phylogenetic tree handling                                    #
#                                                                            #
# Input:                                                                     #
#   - OGU feature table (TSV format e.g. ogu.tsv)                            #
#   - Sample metadata (TSV format)                                           #
#   - Taxonomy file (TSV format)                                             #
#   - Phylogenetic tree file (Optional, Newick format)                       #
#                                                                            #
# Output:                                                                    #
#   - Rarefied phyloseq objects (RDS format)                                 #
#   - Rarefied feature tables (CSV format)                                   #
#   - Summary statistics (TXT format)                                        #
#                                                                            #
# Usage: Rscript read-rarefaction.r  \                                       #
#     -f ogu_woltka_feature_table.tsv \                                      #
#     -m sample_metadata.tsv \                                               #
#     -t lineages.txt \                                                      #
#     -p tree.nwk \                                                          #
#     -o rarefied_output \                                                   #
#     -v                                                                     #
#                                                                            #
# Last updated: 04-02-25                                                     #
##############################################################################

# Options

options(warn=FALSE, max.print=50, width=150, scipen = 999)

# Load required packages
required_packages <- c("phyloseq", "tidyverse", "microbiome", "optparse", "ape")

# Install and load packages
load_packages <- function(packages) {
    suppressPackageStartupMessages({
        new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
        if(length(new_packages)) {
            message("Installing packages: ", paste(new_packages, collapse=", "))
            install.packages(new_packages)
        }
        invisible(lapply(packages, library, character.only = TRUE))
    })
}

load_packages(required_packages)

# Parse command line arguments
option_list <- list(
    make_option(c("-f", "--feature_table"), 
                type="character", 
                help="Path to OGU feature table (TSV format)"),
    make_option(c("-m", "--metadata"), 
                type="character", 
                help="Path to sample metadata (TSV format)"),
    make_option(c("-t", "--taxonomy"), 
                type="character", 
                help="Path to taxonomy file (TSV format)"),
    make_option(c("-p", "--phylotree"), 
                type="character", 
                default=NULL,
                help="Path to phylogenetic tree file (Newick format) [optional]"),
    make_option(c("-o", "--output"), 
                type="character", 
                default="rarefied_output",
                help="Output directory [default=%default]"),
    make_option(c("-v", "--verbose"), 
                action="store_true", 
                default=FALSE,
                help="Print detailed progress information")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to print verbose output
verbose_print <- function(x) {
    if(opt$verbose) cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", x, "\n"))
}

# Function to import and process data
import_data <- function(feature_table_path, metadata_path, taxonomy_path, tree_path = NULL) {
    verbose_print("Importing feature table...")
    ogu_file <- read.table(feature_table_path, header = TRUE, row.names = 1)
    colnames(ogu_file) <- str_remove_all(colnames(ogu_file), "_filtered")
    ogumat <- as.matrix(ogu_file)
    ogu <- otu_table(t(ogumat), taxa_are_rows = FALSE)
    
    verbose_print("Importing metadata...")
    metadata <- read.delim(metadata_path, header = TRUE, row.names = 1)
    meta <- sample_data(metadata)
    
    verbose_print("Importing taxonomy...")
    taxonomy <- read.table(taxonomy_path, header = FALSE, sep = "\t", row.names = 1)
    taxonomy <- taxonomy %>% 
        separate(V2, sep = "; ", #include single space
                c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
        mutate_at(vars(everything()), str_remove,".__")

    tax <- tax_table(as.matrix(taxonomy))
    
    # Verify sample names match
    if(!all(sort(sample_names(ogu)) == sort(sample_names(meta)))) {
        message("Warning: Sample names do not match between feature table and metadata")
    }
        
    # Add phylogenetic tree if provided
    if(!is.null(tree_path)) {
        verbose_print("Importing phylogenetic tree...")
        phyloTree <- read.tree(tree_path)
        seqdata <- phyloseq(ogu, tax, meta, phyloTree)
    } else {
        # Create phyloseq object
        seqdata <- phyloseq(ogu, tax, meta)

    }
    
    return(seqdata)
}

# Function to generate rarefaction curves
create_rarefaction_curves <- function(ogu_matrix, output_dir) {
    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    raremax <- min(rowSums(ogu_matrix))

    verbose_print(paste("Subsample size for rarefying microbiome community:", raremax))

    verbose_print("Generating rarefaction curves...")

    curves <- rarecurve(ogu_matrix, step = 10000, sample = raremax, col = "blue", cex=0.6, label = FALSE, tidy = TRUE)
    
    # Save results
    verbose_print("Saving results...")

    plot <- curves %>%
                ggplot(aes(x = Sample, y = Species)) +
                    geom_line(aes(colour = Site)) +
                    scale_x_continuous(n.breaks = 10) +
                    theme(legend.position = "none") +
                    labs(x = "Number of Reads", y = "Number of Observed Species")

    ggsave(file.path(output_dir, "rarefaction-curves.png"), plot = plot, width = 10, height = 7.5, dpi = 300)

        zoom <- curves %>%
                ggplot(aes(x = Sample, y = Species)) +
                    geom_line(aes(colour = Site)) +
                    coord_cartesian(xlim = c(0,1000000)) +
                    scale_x_continuous(n.breaks = 10) +
                    theme(legend.position = "none") +
                    labs(x = "Number of Reads", y = "Number of Observed Species")

    ggsave(file.path(output_dir, "rarefaction-curves-zoomed.png"), plot = zoom, width = 10, height = 7.5, dpi = 300)

    
    message("Rarefaction curves saved in: ", output_dir)
}

# Function to process and rarefy data
process_and_rarefy <- function(seqdata, output_dir) {
    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Separate Bacteria and Archaea
    verbose_print("Separating Bacteria and Archaea...")
    seqdata_bac <- subset_taxa(seqdata, Kingdom == "Bacteria")
    seqdata_arc <- subset_taxa(seqdata, Kingdom == "Archaea")

    verbose_print("Creating rarefaction curves with subsetted data...")
    
    create_rarefaction_curves(as.data.frame(otu_table(seqdata_bac)), output_dir)

    # Get minimum read depths
    min_depth_bac <- min(rowSums(otu_table(seqdata_bac)))
    min_depth_arc <- min(rowSums(otu_table(seqdata_arc)))
    
    verbose_print(sprintf("Rarefying Bacteria to depth: %d", min_depth_bac))
    verbose_print(sprintf("Rarefying Archaea to depth: %d", min_depth_arc))
    
    # Rarefy
    set.seed(42) # for reproducibility
    seqdata_bac_rare <- rarefy_even_depth(seqdata_bac, 
                                         sample.size = min_depth_bac,
                                         rngseed = 42,
                                         replace = FALSE)
    seqdata_arc_rare <- rarefy_even_depth(seqdata_arc,
                                         sample.size = min_depth_arc,
                                         rngseed = 42,
                                         replace = FALSE)
    
    # Add metadata variables
    sample_data(seqdata_bac_rare)$timepoint <- 
        sample_data(seqdata_bac_rare)$any_other_relevant_information
    sample_data(seqdata_bac_rare)$trt_timept <- 
        paste(sample_data(seqdata_bac_rare)$trt,
              sample_data(seqdata_bac_rare)$timepoint,
              sep = "_")
    
    sample_data(seqdata_arc_rare)$timepoint <- 
        sample_data(seqdata_arc_rare)$any_other_relevant_information
    sample_data(seqdata_arc_rare)$trt_timept <- 
        paste(sample_data(seqdata_arc_rare)$trt,
              sample_data(seqdata_arc_rare)$timepoint,
              sep = "_")
    
    # Save results
    verbose_print("Saving results...")
    saveRDS(seqdata_bac_rare, file.path(output_dir, "rarefied_bacteria.rds"))
    saveRDS(seqdata_arc_rare, file.path(output_dir, "rarefied_archaea.rds"))
    
    # Save feature tables
    write.csv(otu_table(seqdata_bac_rare),
              file.path(output_dir, 
                       sprintf("rarefied%d_bacteria_feature_table.csv", 
                             min_depth_bac)))
    write.csv(otu_table(seqdata_arc_rare),
              file.path(output_dir, 
                       sprintf("rarefied%d_archaea_feature_table.csv", 
                             min_depth_arc)))
    
    # Generate and save summary statistics
    sink(file.path(output_dir, "rarefaction_summary.txt"))
    cat("Rarefaction Summary\n")
    cat("==================\n\n")
    cat("Bacteria:\n")
    print(microbiome::summarize_phyloseq(seqdata_bac_rare))
    cat("\nArchaea:\n")
    print(microbiome::summarize_phyloseq(seqdata_arc_rare))
    sink()
    
    return(list(bacteria = seqdata_bac_rare, 
                archaea = seqdata_arc_rare))
}

# Set plot theme
fontsize <- 14
theme_set(
  theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size = fontsize, face='bold', color='black'), 
      axis.text.y = element_text(size = fontsize, face='bold', color='black'),
      axis.title.x = element_text(size = fontsize, face='bold', color='black'),
      plot.title = element_text(hjust=0.5,face = "bold", size= 26),
      legend.title=element_blank(),
      legend.position = 'bottom',
      legend.text=element_text(size = fontsize, face='bold', color='black'),
      axis.title.y = element_text(size = fontsize, face='bold', color='black'),
      strip.text = element_text(size = fontsize, face = 'bold', color = 'black')
      ))

main <- function() {
    tryCatch({
        # Validate inputs
        if(is.null(opt$feature_table) || is.null(opt$metadata) || is.null(opt$taxonomy)) {
            stop("Feature table, metadata, and taxonomy files are required")
        }
        
        verbose_print("Starting microbiome data processing...")
        
        # Import data
        seqdata <- import_data(opt$feature_table, 
                              opt$metadata, 
                              opt$taxonomy, 
                              opt$phylotree)
        
        # Process and rarefy data
        rarefied_data <- process_and_rarefy(seqdata, opt$output)
        
        verbose_print("Processing complete! Results saved in output directory.")
        
    }, error = function(e) {
        cat("ERROR: ", conditionMessage(e), "\n")
        quit(status = 1)
    })
}

# Execute main function
if(!interactive()) {
    main()
}
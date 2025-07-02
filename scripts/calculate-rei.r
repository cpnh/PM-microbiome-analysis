#!/usr/bin/env Rscript

#######################################
#                                     #
#    Rarefaction Efficiency Index     #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                   #
# GitHub: promicrobial (https://github.com/promicrobial)                       #
# Date: 04-02-25                                                             #
# License: MIT                                                                #
# Version: 1.0                                                               #
#                                                                            #
# Description: This script calculates the Rarefaction Efficiency Index (REI)  #
#              for microbiome data at different sequencing depths. REI helps  #
#              determine optimal rarefaction depth by measuring information    #
#              retention at different depths.                                 #
#                                                                            #
# Dependencies:                                                              #
#   - phyloseq: For handling microbiome data                                #
#   - optparse: For parsing command line arguments                          #
#   - data.table: For efficient data manipulation                           #
#                                                                            #
# Citations:                                                                 #
#   - REI method: Hong et al. (https://github.com/jcyhong/rarefaction)      #
#                                                                            #
# Usage: Rscript calculate_rei.R [options]                                   #
#   Required:                                                               #
#     -i, --input     Input OTU/feature table (TSV format)                 #
#     -m, --metadata  Sample metadata file (TSV format)                    #
#   Optional:                                                              #
#     -d, --depths    Comma-separated list of rarefaction depths           #
#     -o, --output    Output file prefix (default: rei_results)            #
#     -g, --group     Metadata column to use for grouping                  #
#     -v, --verbose   Print detailed progress information                  #
#                                                                          #
# Output:                                                                  #
#   - CSV file with REI values for each specified depth                   #
#   - Summary plot of REI values across depths                            #
#                                                                         #
##############################################################################

# Load required packages
suppressPackageStartupMessages({
    if (!require("optparse")) install.packages("optparse")
    if (!require("phyloseq")) install.packages("phyloseq")
    if (!require("ggplot2")) install.packages("ggplot2")
    if (!require("data.table")) install.packages("data.table")
})

library(optparse)
library(phyloseq)
library(ggplot2)
library(data.table)

# Parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"), 
                type="character", 
                help="Input OTU/feature table (TSV format)"),
    make_option(c("-m", "--metadata"), 
                type="character", 
                help="Sample metadata file (TSV format)"),
    make_option(c("-d", "--depths"), 
                type="character", 
                default="1000000,2000000,3000000", 
                help="Comma-separated list of rarefaction depths [default=%default]"),
    make_option(c("-o", "--output"), 
                type="character", 
                default="rei_results", 
                help="Output file prefix [default=%default]"),
    make_option(c("-g", "--group"), 
                type="character", 
                default="timepoint", 
                help="Metadata column to use for grouping [default=%default]"),
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

# REI calculation function
REI <- function(otu_matrix, rarefied_depth, group, taxa_are_rows=FALSE) {
    if (taxa_are_rows == TRUE) {
        otu_matrix <- t(otu_matrix)
    }
    library_size <- rowSums(otu_matrix)
    
    # Drop all observations with library sizes less than rarefied depth
    otu_matrix <- otu_matrix[library_size >= rarefied_depth, ]
    group <- group[library_size >= rarefied_depth]
    library_size <- library_size[library_size >= rarefied_depth]
    
    # Check if we have enough samples after filtering
    if(length(unique(group)) < 2) {
        warning("Not enough groups after filtering at depth ", rarefied_depth)
        return(NA)
    }
    
    group_lvl <- unique(group)
    n1 <- sum(group == group_lvl[1])
    n2 <- sum(group == group_lvl[2])
    
    # Calculate variances
    var_prop1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                           library_size[group == group_lvl[1]], 
                       2, var)
    var_prop_rrf1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                              library_size[group == group_lvl[1]], 
                          2, function(col) {
                              var(col) + 
                                  1 / rarefied_depth * mean(col * (1 - col) * 
                                                              (library_size[group == group_lvl[1]] - rarefied_depth) / 
                                                              (library_size[group == group_lvl[1]] - 1))
                          })
    
    var_prop2 <- apply(otu_matrix[group == group_lvl[2], ] / 
                           library_size[group == group_lvl[2]], 
                       2, var)
    var_prop_rrf2 <- apply(otu_matrix[group == group_lvl[2], ] / 
                              library_size[group == group_lvl[2]], 
                          2, function(col) {
                              var(col) + 
                                  1 / rarefied_depth * mean(col * (1 - col) * 
                                                              (library_size[group == group_lvl[2]] - rarefied_depth) / 
                                                              (library_size[group == group_lvl[2]] - 1))
                          })
    
    (var_prop1 / n1 + var_prop2 / n2) / 
        (var_prop_rrf1 / n1 + var_prop_rrf2 / n2)
}

main <- function() {
    tryCatch({
        verbose_print("Starting REI calculation...")
        
        # Input validation
        if(is.null(opt$input) || is.null(opt$metadata)) {
            stop("Input and metadata files are required")
        }
        
        # Read input files
        verbose_print("Reading input files...")
        ogu_file <- read.delim(opt$input, header = TRUE, skip = 1)
        metadata <- read.delim(opt$metadata, header = TRUE)
        
        # Prepare OTU table
        rownames(ogu_file) <- ogu_file[,1]
        ogu_file <- ogu_file[, 2:ncol(ogu_file)]
        ogu <- otu_table(ogu_file, taxa_are_rows = TRUE)
        
        # Get grouping variable
        if(!opt$group %in% colnames(metadata)) {
            stop("Grouping column '", opt$group, "' not found in metadata")
        }
        group_var <- metadata[[opt$group]]
        
        # Parse depths
        depths <- as.numeric(strsplit(opt$depths, ",")[[1]])
        
        # Calculate REI for each depth
        verbose_print("Calculating REI for specified depths...")
        results <- data.frame(
            depth = depths,
            rei = sapply(depths, function(d) {
                verbose_print(paste("Processing depth:", d))
                rei <- REI(otu_matrix = ogu, 
                          rarefied_depth = d, 
                          group = group_var, 
                          taxa_are_rows = TRUE)
                mean(rei, na.rm = TRUE)
            })
        )
        
        # Save results
        write.csv(results, 
                 file = paste0(opt$output, "_rei_values.csv"), 
                 row.names = FALSE)
        
        # Create visualization
        p <- ggplot(results, aes(x = depth, y = rei)) +
            geom_line() +
            geom_point() +
            theme_minimal() +
            labs(x = "Rarefaction Depth",
                 y = "Rarefaction Efficiency Index",
                 title = "REI vs Rarefaction Depth") +
            theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(paste0(opt$output, "_rei_plot.pdf"), p)
        
        verbose_print("Analysis complete! Results saved.")
        
    }, error = function(e) {
        cat("ERROR: ", conditionMessage(e), "\n")
        quit(status = 1)
    })
}

# Execute main function
if(!interactive()) {
    main()
}
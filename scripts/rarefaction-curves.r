#!/usr/bin/env Rscript

#######################################
#                                     #
#       Rarefaction Curves            #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                 #
# GitHub: promicrobial (https://github.com/promicrobial)                     #
# Date: 04-02-25                                                             #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: This script performs generates rarefaction curves from        #
#              feature tables e.g. OTU/OGU tables                            #
#                                                                            #
# Dependencies:                                                              #
#   - vegan: rarefaction curve generation                                    #
#   - tidyverse: For data manipulation                                       #
#   - optparse: For command line argument parsing                            #
#                                                                            #
# Input:                                                                     #
#   - OGU feature table (TSV format e.g. ogu.tsv)                            #
#                                                                            #
# Output:                                                                    #
#   - Rarefaction curves (PNG file)                                          #
#                                                                            #
# Usage: Rscript rarefaction-curves.r \                                      #
#     -f ogu_woltka_feature_table.tsv \                                      #
#     -v                                                                     #
#                                                                            #
# Last updated: 04-02-25                                                     #
##############################################################################

# Options

options(warn=FALSE, max.print=50, width=150)

# Load required packages
required_packages <- c("tidyverse", "optparse", "vegan")

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
    make_option(c("-o", "--output"), 
                type="character", 
                default="read-rarefaction",
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

# Function to process and rarefy data
create_rarefaction_curves <- function(ogu_matrix, output_dir) {
    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    raremax <- min(rowSums(ogu_matrix))

    verbose_print(paste("Subsample size for rarefying microbiome community:", raremax))

    verbose_print("Generating rarefaction curves...")

    curves <- rarecurve(ogu_matrix, step = 10000, sample = raremax, col = "blue", cex=0.6, label = FALSE)
    
    # Save results
    verbose_print("Saving results...")

    png(paste(output_dir, "rarefaction-curves", format(Sys.time(), "%b_%d_%Y_%I:%M"), ".png", sep = "_"), width = 1200, height = 800)

    plot(curves)

    dev.off()

}

main <- function() {
    tryCatch({
        # Validate inputs
        if(is.null(opt$feature_table)) {
            stop("Feature table is required")
        }
        
        verbose_print("Starting rarefaction curve generation...")
        
        # Import data
        seqdata <- import_data(opt$feature_table)     
        
        # Draw curves
        create_rarefaction_curves(seqdata, opt$output)
        
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
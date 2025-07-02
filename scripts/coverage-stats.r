#!/usr/bin/env Rscript

#######################################
#                                     #
#        Coverage Statistics          #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                 #
# GitHub: promicrobial (https://github.com/promicrobial)                     #
# Date: 2024-02-13                                                           #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: Processes and analyzes coverage output files from bedtools    #
#              genomecov. Combines multiple coverage files, calculates       #
#              summary statistics, and generates visualizations.             #
#                                                                            #
# Dependencies:                                                              #
#   - tidyverse: For data manipulation and visualization                     #
#   - optparse: For command line argument parsing                            #
#                                                                            #
# Input:                                                                     #
#   - Directory containing bedtools genomecov output files (*.txt)           #
#   Format: chr depth count bases fraction                                   #
#                                                                            #
# Output:                                                                    #
#   - Combined coverage statistics (CSV)                                     #
#   - Coverage summary per sample (CSV)                                      #
#   - Coverage distribution plots (PDF)                                      #
#                                                                            #
# Usage:                                                                     #
#   Rscript coverage-stats.r   \                                             #
#     -i /path/to/coverage/files \                                           #
#     -o coverage_analysis \                                                 #
#     -p "*.txt" \                                                           #
#     -v                                                                     #
##############################################################################

# Load required packages
required_packages <- c("tidyverse", "optparse")

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
    make_option(c("-i", "--input_dir"), 
                type="character",
                help="Directory containing bedtools genomecov output files"),
    make_option(c("-o", "--output_dir"), 
                type="character",
                default="coverage_analysis",
                help="Output directory for results [default=%default]"),
    make_option(c("-p", "--pattern"), 
                type="character",
                default="*.txt",
                help="Coverage file pattern [default=%default]"),
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

# Function to combine coverage files
combine_coverage_files <- function(input_dir, pattern) {
    verbose_print("Reading coverage files...")
    
    # List all coverage files
    coverage_files <- list.files(path = input_dir, 
                               pattern = pattern, 
                               full.names = TRUE)
    
    if(length(coverage_files) == 0) {
        stop("No coverage files found matching pattern: ", pattern)
    }
    
    verbose_print(sprintf("Found %d coverage files", length(coverage_files)))
    
    # Read and combine all files
    coverage_data <- bind_rows(lapply(coverage_files, function(file) {
        verbose_print(sprintf("Processing file: %s", basename(file)))
        
        read_delim(file, 
                  delim = "\t", 
                  col_names = c("feature", "coverage_depth", 
                              "no_bases_with_depth_eq_coverage", 
                              "size_bp", 
                              "fraction_bases_with_depth_eq_coverage"),
                  show_col_types = FALSE) %>%
            mutate(sample = tools::file_path_sans_ext(basename(file)))
    }))
    
    return(coverage_data)
}

# Function to generate summary statistics
generate_summary_stats <- function(coverage_data) {
    verbose_print("Calculating summary statistics...")
    
    coverage_summary <- coverage_data %>%
        group_by(sample) %>%
        summarise(
            mean_coverage = sum(coverage_depth * fraction_bases_with_depth_eq_coverage),
            median_coverage = median(coverage_depth),
            total_bases = first(size_bp),
            bases_covered = sum(no_bases_with_depth_eq_coverage[coverage_depth > 0]),
            percent_covered = bases_covered / total_bases * 100,
            x1_coverage = sum(fraction_bases_with_depth_eq_coverage[coverage_depth >= 1]),
            x5_coverage = sum(fraction_bases_with_depth_eq_coverage[coverage_depth >= 5]),
            x10_coverage = sum(fraction_bases_with_depth_eq_coverage[coverage_depth >= 10])
        )
    
    return(coverage_summary)
}

# Function to create visualizations
create_plots <- function(coverage_data, output_dir) {
    verbose_print("Generating plots...")
    
    # Coverage distribution plot
    p1 <- ggplot(coverage_data %>% 
                     filter(coverage_depth > 0, 
                            coverage_depth <= quantile(coverage_depth, 0.99))) +
        geom_line(aes(x = coverage_depth, 
                      y = fraction_bases_with_depth_eq_coverage, 
                      color = sample)) +
        scale_x_log10() +
        scale_y_log10() +
        labs(x = "Coverage Depth",
             y = "Fraction of Bases",
             title = "Coverage Distribution")
    
    # Cumulative coverage plot
    p2 <- coverage_data %>%
        group_by(sample) %>%
        arrange(coverage_depth) %>%
        mutate(cumulative_fraction = cumsum(fraction_bases_with_depth_eq_coverage)) %>%
        ggplot(aes(x = coverage_depth, y = cumulative_fraction, color = sample)) +
        geom_line() +
        scale_x_log10() +
        labs(x = "Coverage Depth",
             y = "Cumulative Fraction of Bases",
             title = "Cumulative Coverage Distribution")
    
    # Save plots
    ggsave(file.path(output_dir, "coverage_distribution.png"), 
           p1, width = 10, height = 6, units = "in")
    ggsave(file.path(output_dir, "cumulative_coverage.png"), 
           p2, width = 10, height = 6, units = "in")
}

fontsize <- 14

theme_set(
  theme_classic() +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size = fontsize, face='bold', color='black'), 
      axis.text.y = element_text(size = fontsize, face='bold', color='black'),
      axis.title.x = element_text(size = fontsize, face='bold', color='black'),
      plot.title = element_text(hjust=0.5,face = "bold", size= 26),
      legend.title= element_blank(),
      legend.position = 'none',
      legend.text=element_text(size = fontsize, face='bold', color='black'),
      axis.title.y = element_text(size = fontsize, face='bold', color='black'),
      strip.text = element_text(size = fontsize, face = 'bold', color = 'black'),
      strip.background.x = element_rect(fill = "white"),
      strip.background.y = element_rect(fill = "lightgrey"),
      plot.margin = margin(2, 2, 2, 2, "cm")
      ))

main <- function() {
    tryCatch({
        # Validate input
        if(is.null(opt$input_dir)) {
            stop("Input directory is required")
        }
        
        # Create output directory
        dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)
        
        # Combine coverage files
        coverage_data <- combine_coverage_files(opt$input_dir, opt$pattern)
        
        # Generate summary statistics
        coverage_summary <- generate_summary_stats(coverage_data)
        
        # Create visualizations
        create_plots(coverage_data, opt$output_dir)
        
        # Save results
        write_csv(coverage_data, 
                 file.path(opt$output_dir, "combined_coverage_stats.csv"))
        write_csv(coverage_summary, 
                 file.path(opt$output_dir, "coverage_summary.csv"))
        
        # Print summary
        cat("\nCoverage Summary:\n")
        print(coverage_summary)
        
        verbose_print("Analysis complete! Results saved in output directory.")
        
    }, error = function(e) {
        cat("ERROR: ", conditionMessage(e), "\n")
        quit(status = 1)
    })
}

# Execute main function
if(!interactive()) {
    main()
}
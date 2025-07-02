#!/bin/bash

#######################################
#                                     #
#              run_fastp              #
#                                     #
#######################################

##############################################################################
# Author: Nathaniel Cole (nc564@cornell.edu)                                 #
# GitHub: promicrobial (https://github.com/promicrobial)                     #   
# Date: 25-11-24                                                             #
# License: MIT                                                               #
# Version: 1.0                                                               #
#                                                                            #
# Description: fastp processing on R1/R2 or interleaved FASTQ inputs         #
#                                                                            #
# Dependencies:                                                              #
#   - fastp v0.23.2: For quality control and adapter trimming                #
#     Paper: https://doi.org/10.1093/bioinformatics/bty560                   #
#     Tool: https://github.com/OpenGene/fastp                                #
#                                                                            #
# Usage: ./run_fastp.sh <r1_file> <r2_file> <outdir> <base_name>                              #
#                                                                            #
# Last updated: 10-01-25                                                     #
##############################################################################

# Strict error handling
# -e: Exit immediately if a command exits with non-zero status
# -u: Treat unset variables as an error
# -o pipefail: Return value of a pipeline is the status of the last command to exit with a non-zero status

set -euo pipefail

################################################################################
# Help                                                                         #
################################################################################

# Display the help message and exit
usage() {
    cat << EOF
Usage: $(basename "$0") <r1_file> <r2_file> <outdir> <base_name>

This script runs fastp preprocessing on FASTQ files.

Arguments:
  r1_file    Path to the R1 FASTQ file (required)
  r2_file    Path to the R2 FASTQ file (use "NA" for single-end data)
  dir     Path output directory
  base_name  Base name for output files

Options:
  -h, --help    Show this help message and exit

Example:
  $(basename "$0") sample_R1.fastq.gz sample_R2.fastq.gz /path/to/tmp sample_name
  $(basename "$0") single_end.fastq.gz NA /path/to/tmp sample_name

Note:
  For single-end data, use "NA" as the r2_file argument.
EOF
}

# if help is requested
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    usage
    exit 0
fi

################################################################################
# Validations                                                                  #
################################################################################

# validate if correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Error: Incorrect number of arguments." >&2
    usage
    exit 1
fi

r1_file="$1"
r2_file="$2"
outdir="$3"
base_name="$4"

# validate input files exist
if [[ ! -f "$r1_file" ]]; then
    echo "Error: R1 file does not exist: $r1_file" >&2
    exit 1
fi

if [[ "$r2_file" != "NA" ]] && [[ ! -f "$r2_file" ]]; then
    echo "Error: R2 file does not exist: $r2_file" >&2
    exit 1
fi

# validate outdir exist, else create
if [[ ! -d "$outdir" ]]; then
    echo "Warning: Output directory, $outdir does not exist. Creating." >&2
    mkdir -p $outdir
fi

################################################################################
# Main script body                                                             #
################################################################################

fastp_out_r1="${outdir}/$(basename "$base_name")_fastp_R1.fastq.gz"
fastp_out_r2="${outdir}/$(basename "$base_name")_fastp_R2.fastq.gz"
fastp_json="${outdir}/$(basename "$base_name")_fastp.json"
fastp_html="${outdir}/$(basename "$base_name")_fastp.html"

if [[ -e "$r2_file" ]]; then
    # Paired-end fastp command
    fastp -i "$r1_file" -I "$r2_file" \
          -o "$fastp_out_r1" -O "$fastp_out_r2" \
          -j "$fastp_json" -h "$fastp_html" \
          --detect_adapter_for_pe \
          --correction --cut_right --thread 16 \
          --length_required 50 --qualified_quality_phred 20 \
          --unqualified_percent_limit 40 --cut_window_size 4
else
    # Single-end fastp command
    fastp -i "$r1_file" \
          -o "$fastp_out_r1" \
          -j "$fastp_json" -h "$fastp_html" \
          --cut_right --thread 8 \
          --length_required 50 --qualified_quality_phred 20 \
          --unqualified_percent_limit 40 --cut_window_size 4
fi

# Check if fastp was successful
if [ $? -ne 0 ]; then
    echo "Error: fastp preprocessing failed for $base_name" >&2
    exit 1
fi

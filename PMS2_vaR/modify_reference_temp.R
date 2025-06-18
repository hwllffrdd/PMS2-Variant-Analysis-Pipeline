# Description: Run the PMS2_vaR pipeline
# USAGE: Rscript runPMS2_vaR.R [-t tools_file] [-d datasets_file]

#libs
library(yaml)
library("optparse")
library("Biostrings")
library("universalmotif")
library(dplyr)
library(stringr)
# Check runPMS2_vaR.R is being called from PMS2_vaR folder

if (length(list.files(pattern = "modify_reference.R"))== 0){
  cat("Sorry, modify_reference.R .R should be called from PMS2_vaR.R  folder\n")
  quit()
}


# Build options list
#/media/emunte/Elements1/CNV_benchmark/fasta/GRCh37_Ensembl_67.fa
#/media/emunte/Elements1/PMS2_project/PMS2CL/Homo_sapiens_pms2CL.fa
option_list <- list(
  make_option(c("-r", "--reference"), type="character", default="",
              help="Full path to reference fasta(fa)", metavar="character"),
  make_option(c("-p", "--pms2CLfasta"), type="character", default="",
              help="Full path to PMS2CL fasta file", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default="",
              help="Output directory to store results", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
args <- parse_args(opt_parser);


# Get the line where chromosome 7 and 8 start in the reference genome
# Get the line where chromosome 7 and 8 start in the reference genome
cmd1 <- paste('grep -n ">chr7$"',  args$reference)
chr7 <- system(cmd1, intern=TRUE) %>% stringr::str_extract("[0-9]+") %>% as.numeric()
cmd2 <- paste('grep -n ">chr8$"',  args$reference)
chr8 <- system(cmd2, intern=TRUE) %>% stringr::str_extract("[0-9]+") %>% as.numeric()

# Create output directory
direct <- file.path(args$outputdir, "modified_reference")
dir.create(direct, showWarnings = FALSE)

# Get chr7 sequence
system(paste("awk 'NR==", chr7, ", NR==", chr8-1,"' ", args$reference, ">", file.path(direct, "chr7.fa")))

# Read the chr7 fasta file
genome <- readDNAStringSet(file.path(direct, "chr7.fa"))
ref2 <- paste0(genome$chr7, collapse = "")

# Define PMS2CL coordinates in hg38
base <- c(6731863, 6743066)
cat("Using PMS2CL coordinates:", base[1], "-", base[2], "\n")

# Calculate line numbers for the modification
start <- round(base[1] / 60 - 4)
end <- round(base[2] / 60 + 4)

# Get coordinates in the whole reference file
start.whole <- chr7 + start
end.whole <- chr7 + end
n.number <- end.whole - start.whole + 1

# Make the files
cmd1 <- paste0("awk 'NR <", start.whole, "' ", args$reference, " > ", file.path(args$outputdir, "reference_part1.fa"))
cmd2 <- paste0("awk 'NR >", end.whole, "' ", args$reference, " > ", file.path(args$outputdir, "reference_part3.fa"))
system(cmd1)
system(cmd2)

# Create N-sequence file
N.file <- paste0(rep(paste0(paste0(rep("N", 60), collapse = ""), "\n"), n.number), collapse = "")
N.file <- stringr::str_sub(N.file, 1, stringr::str_length(N.file)-1)
write.table(N.file, file.path(args$outputdir, "N.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)

# Merge all files
cmd3 <- paste("cat", file.path(args$outputdir, "reference_part1.fa"), file.path(args$outputdir, "N.txt"), file.path(args$outputdir, "reference_part3.fa"), ">", file.path(direct, paste0(tools::file_path_sans_ext(basename(args$reference)), "_without_PMS2CL.fa")))
system(cmd3)

# Cleanup
unlink(file.path(args$outputdir, "reference_part1.fa"))
unlink(file.path(args$outputdir, "reference_part3.fa"))
unlink(file.path(args$outputdir, "N.txt"))
unlink(file.path(direct, "chr7.fa"))

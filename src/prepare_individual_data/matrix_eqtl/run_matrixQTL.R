#!/usr/bin/env Rscript
# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(argparser, quietly=TRUE)
library(MatrixEQTL, quietly=TRUE)

p <- arg_parser("Run MatrixEQTL")
p <- add_argument(p, "gene_exp", help="")
p <- add_argument(p, "covars", help="")
p <- add_argument(p, "snps", help="")
p <- add_argument(p, "prefix", help="")
p <- add_argument(p, "--model", short='-m', default='linear', help="")
p <- add_argument(p, "--threshold", short='-t', default=1, help="")
p <- add_argument(p, "--output_dir", short="-o", default='.', help="")

argv <- parse_args(p)

stopifnot(argv$model %in% c('linear', 'anova', 'linear_cross'))

if(argv$model == 'linear') {
    useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
} else if (argv$model == 'anova') {
    useModel = modelANOVA
} else {
    useModel = modelLINEAR_CROSS
}

# Genotype file name
SNP_file_name = argv$snps 
#snps_location_file_name = argv$snps_loc 

expression_file_name = argv$gene_exp 
#gene_location_file_name = argv$gene_loc 

covariates_file_name = argv$covars 

# Output file name
output_file_name = file.path(argv$output_dir, paste(argv$prefix, '.assoc.txt', sep=''))

# Only associations significant at this level will be saved
pvOutputThreshold = argv$threshold 

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name = output_file_name,
pvOutputThreshold = pvOutputThreshold,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

## Results:

cat('Analysis for ', argv$prefix, ' done in: ', me$time.in.sec, ' seconds', '\n\n');

#!/usr/bin/env Rscript

# Author: James A. Perez 
# Description: Run colocalization analysis on a distinct eQTL and GWAS locus 

suppressMessages(library(argparser))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(data.table))
suppressMessages(library(tools))
suppressMessages(library(stringr))

random_string <- function(n=1, len=6)
{
    randomString <- c(1:n)
    for (i in 1:n)
    {
            randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
            len, replace=TRUE),
            collapse="")
    }
    return(randomString)
}

p <- arg_parser("Perform multivariable MR")
p <- add_argument(p, "exposureData", help="")
p <- add_argument(p, "outcomeData", help="")
p <- add_argument(p, "outcomeID", help="")
p <- add_argument(p, "prefix", help="")
p <- add_argument(p, "--outdir", short="-o", type="character", default="./", help="")

args <- parse_args(p)

if(!dir.exists(args$outdir)){
    dir.create(args$outdir)
}

outcome_dat = fread(args$outcomeData, sep='\t')
exposure_dat = fread(args$exposureData, sep='\t') 

# 1. Format exposure data
## keep major allele SNP codes for later use

snp_key = data.frame(SNP_ID=exposure_dat$SNP)
snp_key$SNP = sapply(exposure_dat$SNP, function(x) paste(strsplit(x=x, split=':')[[1]][1],  strsplit(x=x, split=':')[[1]][2], sep=':')) 
snp_key$BP = sapply(exposure_dat$SNP, function(x) paste(strsplit(x=x, split=':')[[1]][2])) 

exposure_dat$SNP = snp_key$SNP 
exposure_dat = format_data(exposure_dat, type='exposure', snp_col='SNP', beta_col='beta', se_col='se', effect_allele_col='A1', other_allele_col='A2', pval='p', gene_col='gene', phenotype_col='gene')
cat(dim(exposure_dat)[1], ' associations in exposure data\n')

# 2. Format outcome data

outcome_dat_subset = outcome_dat[which(outcome_dat$Bp %in% snp_key$BP),] 
if(dim(outcome_dat_subset)[1]==0){
    cat('No SNPs in common with outcome...exiting')
    quit()
}

outcome_dat_subset$SNP = sapply(outcome_dat_subset$SNP, function(x) paste(strsplit(x=x, split=':')[[1]][1],  strsplit(x=x, split=':')[[1]][2], sep=':'))
outcome_dat_subset = outcome_dat_subset[with(outcome_dat_subset, order(Chr, Bp, Pvalue)),] 

if(any(duplicated(outcome_dat_subset$SNP))){
   cat('Duplicated SNPs detected. Removing the following...')
   print(outcome_dat_subset[duplicated(outcome_dat_subset$SNP),])
}

outcome_dat_subset = outcome_dat_subset[!duplicated(outcome_dat_subset[,c('Chr', 'Bp')]),]
outcome_dat_subset$phenotype = args$outcomeID
    
outcome_dat_subset = format_data(outcome_dat_subset, type='outcome', snp_col='SNP', beta_col='Beta', se_col='SE', effect_allele_col='Effect_allele', other_allele_col='Non_Effect_allele', pval='Pvalue', phenotype_col='phenotype')                         
    
# 3. Harmonize outcome and exposure IDs 
   
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat_subset, action=2)
dat <- subset(dat, mr_keep)

cat(dim(dat)[1], ' harmonised associations\n')

## keep a copy of full SNP code for plink input
dat$SNP_ID = snp_key[match(dat$SNP, snp_key$SNP), c('SNP_ID')]

outfile = file.path(args$outdir, paste0(args$prefix,'_harmonise.txt'))
write.table(dat, file=outfile, sep='\t', col.names=T, row.names=F, quote=F)

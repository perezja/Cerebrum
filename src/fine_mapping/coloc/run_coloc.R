#!/usr/bin/env Rscript

# Author: James A. Perez 
# Description: Run colocalization analysis on a distinct eQTL and GWAS locus 

suppressMessages(library(argparser))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(coloc))
suppressMessages(library(data.table))
suppressMessages(library(tools))
suppressMessages(library(stringr))

p <- arg_parser("Perform bayesian COLOC analysis.")
p <- add_argument(p, "locusFile", help="")
p <- add_argument(p, "eQTLFreqFile", help="")
p <- add_argument(p, "outcomeData", help="")
p <- add_argument(p, "outcomeID", help="")
p <- add_argument(p, "--outFile", short="-o", type="character", default="coloc.out", help="")

args <- parse_args(p)


# 1. Read association data

cat('Processing ', args$locusFile, '\n')
cat('Reading: ', args$outcomeData, '\n')
outcome_dat = fread(args$outcomeData, sep='\t')

cat('Reading: ', args$locusFile, '\n')
locus_dat = fread(args$locusFile, sep='\t')

# 2. Get MAF for SNPs 

maf = fread(args$eQTLFreqFile)

locus_dat = locus_dat[complete.cases(locus_dat),]
locus_dat = merge(locus_dat, maf, on='SNP') 
locus_cols = names(locus_dat)[!(names(locus_dat)==c('NCHROBS'))]
locus_dat = locus_dat[,..locus_cols]

# 3a. format GWAS data

cat('Formatting outcome data...', '\n')

## Remove allele from SNP codes
## (effect and non-effect alleles swapped in the SNP codes
## will prevent exposure and outcome SNPs to harmonise) 

outcome_dat_subset$SNP = sapply(outcome_dat_subset$SNP, function(x) paste(strsplit(x=x, split=':')[[1]][1],  strsplit(x=x, split=':')[[1]][2], sep=':'))

## Among multiallelic sites with duplicate chrom and position, keep the variants  
## with a higher significance in the GWAS outcome data 

outcome_dat_subset = outcome_dat[which(outcome_dat$Bp %in% locus_dat$BP),] 
if(dim(outcome_dat_subset)[1]==0){ 
    message('No SNPs in common after harmonization...skipping.')
    next
}
outcome_dat_subset = outcome_dat_subset[with(outcome_dat_subset, order(Chr, Bp, Pvalue)),] 

if(any(duplicated(outcome_dat_subset$SNP))){
    cat('Duplicated SNPs detected. Removing the following...')
    print(outcome_dat_subset[duplicated(outcome_dat_subset$SNP),])
}

outcome_dat_subset = outcome_dat_subset[!duplicated(outcome_dat_subset[,c('Chr', 'Bp')]),]
outcome_dat_subset$phenotype = args$outcomeID

## prepare dataframe for harmonization

outcome_dat_subset_ = format_data(outcome_dat_subset, type='outcome', snp_col='SNP', beta_col='Beta', se_col='SE', effect_allele_col='Effect_allele', other_allele_col='Non_Effect_allele', pval='Pvalue', phenotype_col='phenotype')                         
# 3b. format eQTL data

cat('Formatting eQTL data...', '\n')

## remove alleles from SNP codes

locus_dat$SNP = sapply(locus_dat$SNP, function(x) paste(strsplit(x=x, split=':')[[1]][1],  strsplit(x=x, split=':')[[1]][2], sep=':'))

names(locus_dat) <- c('Chr', 'SNP', 'Bp', 'id.exposure', 'beta.exposure', 'se.exposure', 'pval.exposure', 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure')

# 4. Harmonize GWAS and eQTL associations 
   
## exit if no SNPs in common 

if(!any(locus_dat$SNP %in% outcome_dat_subset_$SNP)){
    message(paste('Error: Locus has no SNPs in common...'))
    quit(status=1)
}

## haromise

coloc_dat = harmonise_data(exposure_dat = locus_dat, outcome_dat = outcome_dat_subset_, action=2)

# remove ambiguous palindromic alleles

coloc_dat = coloc_dat[!(coloc_dat$palindromic & coloc_dat$ambiguous),]
if(dim(coloc_dat)[1]==0){ 
    message('No SNPs in common after harmonization...skipping.')
    quit(status=1)
}

## reorder GWAS data according to haromised dataframe since we will use
## case proprtion and sample size columns for ```coloc.abf```
outcome_dat_subset <- outcome_dat_subset[match(coloc_dat$SNP, outcome_dat_subset$SNP),] 
stopifnot(all(outcome_dat_subset$SNP==coloc_dat$SNP))

# 4. Perform ABF Colocalization

res.tab <- coloc.abf(
        dataset1=list(beta=coloc_dat$beta.outcome, varbeta=coloc_dat$se.outcome^2, N=outcome_dat_subset$Nsamples, type='cc', s=outcome_dat_subset$Proportion_cases), 
        dataset2=list(beta=coloc_dat$beta.exposure, varbeta=coloc_dat$se.exposure^2, MAF=coloc_dat$eaf.exposure, sdY=1, type='quant'))

abf.res <- as.data.frame(as.list(res.tab[[1]])) 

locus_name = basename(args$locusFile)
ls_idx = which(locus_dat$pval.exposure==min(locus_dat$pval.exposure))
lead_snp = locus_dat[ls_idx,]$SNP[1]
lead_snp_genes = unique(locus_dat[ls_idx,]$`id.exposure`)

desc_row = list(locus_name=locus_name, lead_snp=lead_snp, lead_snp_gene=str_c(lead_snp_genes, collapse=','), outcome=args$outcomeID)
abf.res <- cbind(desc_row, abf.res)

coloc_res_outfile = args$outFile
coloc_harm_outfile = paste0(file.path(args$outFile), locus_name, '.harmonised.txt')
write.table(abf.res, file=coloc_res_outfile, sep='\t', col.names=T, row.names=F, quote=F)
write.table(coloc_dat, file=coloc_harm_outfile, sep='\t', col.names=T, row.names=F, quote=F)

print('done')

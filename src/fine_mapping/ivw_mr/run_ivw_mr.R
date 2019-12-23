#!/usr/bin/env Rscript

# Author: James A. Perez 
# Description: Run colocalization analysis on a distinct eQTL and GWAS locus 

suppressMessages(library(argparser))
suppressMessages(library(TwoSampleMR))
suppressMessages(library(data.table))
suppressMessages(library(tools))
suppressMessages(library(stringr))

p <- arg_parser("Perform multivariable MR")
p <- add_argument(p, "harmonised_dat", help="")
p <- add_argument(p, "prefix", help="")
p <- add_argument(p, "--sensitivity", flag=TRUE, help="")
p <- add_argument(p, "--outdir", short="-o", type="character", default="./", help="")

args <- parse_args(p)

dat = fread(args$harmonised_dat, sep='\t')

if(!dir.exists(args$outdir)){
    dir.create(args$outdir)
}

## 1. Obtain heterogeneity estimates

het_dat <- mr_heterogeneity(dat)
het_dat <- het_dat[het_dat$method=='Inverse variance weighted',]

# 2. Perform MR 
  
res <- mr(dat)

res <- merge(res, het_dat, by=c('id.exposure', 'id.outcome', 'outcome', 'exposure', 'method'), all.x=T)
res <- res[res$method %in% c('Wald ratio', 'Inverse variance weighted'), !(names(res) %in% c('id.outcome', 'id.exposure'))] 


mr_outfile = file.path(args$outdir, paste0(args$prefix,'_mr.txt'))
write.table(res, file=mr_outfile, sep='\t', col.names=T, row.names=F, quote=F)

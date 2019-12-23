library(argparser, quietly=TRUE)
library(data.table)

p <- arg_parser("Order participants genotypes to match expression matrix.")  
p <- add_argument(p, "genotype_dir", help="")
p <- add_argument(p, "expression_dir", help="")
p <- add_argument(p, "--output_dir", short="-o", help="")
args <- parse_args(p)

#rids = c('BM22', 'TCTX', 'DLPFC', 'PFCTX', 'PAR')
rids = c('TCTX', 'PFCTX', 'PAR')

stopifnot(dir.exists(args$genotype_dir) | dir.exists(args$expression_dir))

expFiles <- list.files(args$expression_dir)
expFiles <- expFiles[grepl('*expression.txt.gz$', expFiles)]
prefix = sapply(expFiles, function(x) gsub(".expression.txt.gz", "", x, fixed=TRUE))
expFiles <- expFiles[sapply(prefix, function(x) x %in% rids)]
expPrefix = sapply(expFiles, function(x) gsub(".expression.txt.gz", "", x, fixed=TRUE))

genoFiles <- list.files(args$genotype_dir)
genoFiles <- genoFiles[grepl('*.snps.txt$', genoFiles)]
prefix = sapply(genoFiles, function(x) gsub(".snps.txt", "", x, fixed=TRUE))
genoFiles <- genoFiles[sapply(prefix, function(x) x %in% rids)]
genoPrefix = sapply(genoFiles, function(x) gsub(".snps.txt", "", x, fixed=TRUE))

expFiles <- expFiles[match(genoPrefix, expPrefix)]

for(i in 1:length(rids)){

print(expFiles[i]); print(genoFiles[i])

exp = fread(file.path(args$expression_dir, expFiles[i]), sep='\t') 
exp = as.data.frame(exp)
pids = names(exp)[2:ncol(exp)]
rm(exp)

snps = fread(file.path(args$genotype_dir, genoFiles[i]), sep='\t')
snps = as.data.frame(snps)

label = gsub(".snps.txt", "", genoFiles[i], fixed=TRUE)
print(label)
if(label=='PAR'){
names(snps)[2:ncol(snps)] = gsub(pattern='_((?=MAP)|(?=DIAN)|(?=NIALOAD)).+', replacement='', x=names(snps)[2:ncol(snps)], perl=TRUE)
} else if(label=='TCTX') {
names(snps)[2:ncol(snps)] = gsub(pattern='(?<=Omni25Exome).+', replacement='', x=names(snps)[2:ncol(snps)], perl=TRUE)
} else {
names(snps)[2:ncol(snps)] = gsub(pattern='(?<=HiSeqWGS).+', replacement='', x=names(snps)[2:ncol(snps)], perl=TRUE)
}
snpName = snps$SNP
snps = snps[,2:ncol(snps)]
snps <- snps[,c(match(pids, names(snps)))]

stopifnot(names(snps) == pids)

snps = cbind(SNP=snpName, snps)
rownames(snps) <- NULL

outfile = file.path(args$output_dir, paste0(label, '.snps.txt')) 
write.table(snps, file=outfile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
print(paste0('wrote to: ', outfile))

}

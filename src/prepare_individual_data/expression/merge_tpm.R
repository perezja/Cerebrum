#!/usr/bin/env Rscript
# Author: James A. Perez

library(argparser, quietly=TRUE)

WriteTable <- function(data, filename, index.name) {
    datafile <- gzfile(filename, open = "wt")
    on.exit(close(datafile))
    header <- c(index.name, colnames(data))

    # write shape as first line
    writeLines(paste0(dim(data), collapse="\t"), con=datafile, sep="\n")
    # then write header
    writeLines(paste0(header, collapse="\t"), con=datafile, sep="\n")
    write.table(format(data, digits=6), datafile, sep="\t", col.names=FALSE, quote=FALSE)
}

p <- arg_parser("Merge quantified files from Salmon")
p <- add_argument(p, "paths", help="Newline-delimited file listing paths.")
p <- add_argument(p, "sample_ids", help="Newline-delimited file listing sample ids.")
p <- add_argument(p, "gtf", help="Path to gtf file.")
p <- add_argument(p, "--outfile", short="-o", help="Output file name.", default='tpm.gct.gz')

argv <- parse_args(p)

library(tximport, quietly=TRUE)
library(rtracklayer, quietly=TRUE)

sample_ids = scan(argv$sample_ids, character())
paths = scan(argv$paths, character())

merge_tpm <- function(files, sample_ids, tx2gene)
{
    txi <- tximport(files=files,
   	            type='salmon',
 		    tx2gene=tx2gene,
                    countsFromAbundance='no')
    txi_mat <- summarizeToGene(txi,
	           tx2gene,
		   ignoreTxVersion=FALSE,
		   ignoreAfterBar=FALSE,
		   countsFromAbundance='no')
    tpm <- data.frame(txi_mat[['abundance']])
    colnames(tpm) <- c(sample_ids);
    return(tpm)

}

gff_obj <- readGFF(argv$gtf)
tx2gene = data.frame(gff_obj$transcript_id, gff_obj$gene_id)

df <- merge_tpm(paths, sample_ids, tx2gene)
rownames(df) <- sapply(rownames(df), function(x) sub('\\.[0-9]+','',x))
df <- df[order(row.names(df)),]

WriteTable(df, argv$outfile, 'gene_id') 

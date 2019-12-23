#!/usr/bin/env python3
# Author: James A. Perez 

import argparse
import subprocess
import pandas as pd
import numpy as np
import os
import json

def create_esd_files():
    """ Fragments gzipped SMR formatted eqtl association data into gene-level chunks.
    """
    
    of_prefix = os.path.join(args.esd_dir, args.prefix+'_')
    cmd = "zcat {0} | awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\";}} NR==1{{ $1=\"\"; hdr=substr($0,2);next }} !($1 in a) {{ print hdr > \"{1}\"$1\".esd\"; a[$1]++}} NR>1 {{ gene=$1; $1=\"\"; print substr($0,2) > \"{1}\"gene\".esd\"}}'".format(args.assoc_file, of_prefix)
    subprocess.check_call(cmd, shell=True)
   
parser = argparse.ArgumentParser(description='Prepare input for Metasoft.')
parser.add_argument('assoc_file', help='Path to eQTL-association file in ESD format with added gene column.')
parser.add_argument('esd_dir', help='Path to eQTL-association file in ESD format with added gene column.')
parser.add_argument('gene_bed', help='TSS 4-BED file (chr, start, end, gene).')
parser.add_argument('ensembl2gene', help='File mapping ensembl to gene symbols.')
parser.add_argument('prefix', help='Prefix for esd and flist files.')
parser.add_argument('--use_esd', action='store_true', help='Use the a pre-built esd directory of exposure files.')
parser.add_argument('-o', '--output_dir', default='.', help='Directory for output flist file.')
args = parser.parse_args() 

args.esd_dir = os.path.abspath(args.esd_dir)

def main():

    if not os.path.exists(args.esd_dir):
        os.makedirs(args.esd_dir)
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    if not args.use_esd:
        print('Making *.esd files...')
        create_esd_files()

    else:
        print('skipping *.esd build...')

    def gene_id(fn):
        return(fn.split('_')[-1].split('.')[0])
    
    gene_ids = [gene_id(x) for x in os.listdir(args.esd_dir)]
    if not gene_ids[0].startswith('ENSG'):
        print('ERROR: gene ids not properly being extracted from directory: {}'.format(args.esd_dir))
        quit()

    with open(args.ensembl2gene) as fp:
        ensbl2gene = {row.split('\t')[0]: row.split('\t')[1] for row in fp.read().strip().split('\n')}

    # gene bed must have header 
    gene_bed = pd.read_csv(args.gene_bed, sep='\t', index_col='gene_id', skiprows=1, names=['chr', 'tss', 'tssp1', 'gene_id', 'strand'])
    gene_bed = gene_bed[gene_bed.index.isin(gene_ids)]
    gene_bed = gene_bed[~gene_bed.index.duplicated()]

    if gene_bed['chr'].str.contains('chr').any():
        print('stripping "chr" prefix from chromosome names...')
        gene_bed['chr'] = gene_bed['chr'].strip('chr')

    # biomaRt gives strand as 1/-1
    gene_bed['orientation'] = gene_bed['strand'].apply(lambda x: '+' if x==1 else '-')

    assert gene_bed.shape[0] == len(gene_ids)
    
    # flist: Chr, ProbeID, GeneticDistance, ProbeBp, Gene, Orientation, PathOfEsd

    flist = pd.DataFrame(0, index=gene_ids, columns=['Chr', 'ProbeID', 'GeneticDistance', 'ProbeBp', 'Gene', 'Orientation', 'PathOfEsd'])
    for gene_id, gseries in gene_bed.iterrows():

        esd_fp = os.path.join(args.esd_dir, args.prefix+'_'+gene_id+'.esd')

        assert os.path.exists(esd_fp)

        if gene_id not in ensbl2gene.keys():
            gene_name = 'NA'
        else:
            gene_name = ensbl2gene[gene_id]

        gene_chr = gseries['chr']
        gene_bp = gseries['tss']
        gene_strand = gseries['orientation']

	# all tss coordinates should have already been flipped for (-) stand gene tss's
        flist.loc[gene_id] = [gene_chr, gene_id, np.nan, gene_bp, gene_name, gene_strand, esd_fp]
        
    flist_outfile = os.path.join(args.output_dir, args.prefix+'.flist')
    flist.to_csv(flist_outfile, sep='\t', na_rep=0, index=False)

    print('wrote *.flist to: {}'.format(flist_outfile))

if __name__ == "__main__":
    main()


#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import glob
from collections import defaultdict
import json
import os
import gzip

parser = argparse.ArgumentParser(prog='.')
parser.add_argument('cis_metasoft', type=str, help='METASOFT file with all cis associations')
parser.add_argument('index_snps', type=str, help='GWAS Index SNP file.')
parser.add_argument('prefix', type=str,  help='Prefix for outfile.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')
parser.add_argument('-c', '--chunksize', type=int, default=1e6, help=".")

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

def locus_name(prefix, locus_row):

    chrom=str(locus_row['CHR'])
    st=str(locus_row['STARTBP'])
    en=str(locus_row['ENDBP'])

    return('.'.join([prefix,chrom,st,en]))

def segment_by_locus(cis_file, gwas_df):

    for i,eqtl_chunk_df in enumerate(pd.read_csv(cis_file, sep='\t', chunksize=args.chunksize)):

        print('\rProcessing chunk {}'.format(i+1), end='', flush=True)
    
        eqtl_chunk_df['CHR'] = eqtl_chunk_df['VARIANT_ID'].map(lambda x: int(x.split(',')[1].split(':')[0]))
        eqtl_chunk_df['BP'] = eqtl_chunk_df['VARIANT_ID'].map(lambda x: int(x.split(',')[1].split(':')[1]))
        eqtl_chunk_df['SNP'] = eqtl_chunk_df['VARIANT_ID'].map(lambda x: str(x.split(',')[1]))
        eqtl_chunk_df['GENE'] = eqtl_chunk_df['VARIANT_ID'].map(lambda x: x.split(',')[0])

        for j,gwas_locus in gwas_df.iterrows():

            mask = (
                (eqtl_chunk_df['BP'] < gwas_locus['ENDBP']) &
                (eqtl_chunk_df['BP'] > gwas_locus['STARTBP']) &
                (eqtl_chunk_df['CHR'] == gwas_locus['CHR'] )
            ).values
            
            assert(len(mask) == eqtl_chunk_df.shape[0])

            eqtl_coloc_df = eqtl_chunk_df[mask]

            if eqtl_coloc_df.empty:
                continue

            eqtl_coloc_df = eqtl_coloc_df[['CHR', 'BP', 'SNP', 'GENE', 'BETA_FE', 'STD_FE', 'PVALUE_FE']]  
            eqtl_coloc_df.rename({'BETA_FE':'BETA', 'STD_FE':'SE', 'PVALUE_FE':'P'}, inplace=True, axis=1)

            locus_outfile = os.path.join(args.output_dir, locus_name(args.prefix, gwas_locus))

            if os.path.exists(locus_outfile):
                eqtl_coloc_df.to_csv(locus_outfile, sep='\t', index=False, header=False, mode='a')
            else:
                eqtl_coloc_df.to_csv(locus_outfile, sep='\t', index=False, header=True)

def main():

    gwas_df = pd.read_csv(args.index_snps, sep='\t')
    segment_by_locus(args.cis_metasoft, gwas_df)
    print('wrote locus files to: {}'.format(args.output_dir))

if __name__ == '__main__':
    main()

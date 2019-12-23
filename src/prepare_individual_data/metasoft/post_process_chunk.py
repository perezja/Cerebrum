#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import glob
from collections import defaultdict
import json
import os
import gzip

parser = argparse.ArgumentParser(prog='Stratify meta-analysis results by pvalue and label cis/trans associations.')
parser.add_argument('metasoft_file', type=str, help='')
parser.add_argument('tss_bed', type=str, help='TSS bed which defines reference for cis-qtls.')
parser.add_argument('prefix', type=str, help='Prefix for outfile.')
parser.add_argument('-p', '--pthreshold', default=5e-8, type=float, help='pvalue threshold.')
parser.add_argument('-d', '--cis_distance', type=int, default=1000000, help='Window size to call cis QTLs.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')
parser.add_argument('-c', '--chunksize', type=int, default=1e6, help=".")
parser.add_argument('--debug', action='store_true', help='')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

def stratify_by_distance(gene2tss, metasoft_file, header):

    # concatenate chunks
    print('Loading chunks for {}'.format(args.prefix))

    re2_outfile = os.path.join(args.output_dir, '.'.join([args.prefix,'re2.txt']))
    fe_outfile = os.path.join(args.output_dir, '.'.join([args.prefix,'fe.txt']))
    joint_outfile = os.path.join(args.output_dir, '.'.join([args.prefix,'joint.txt']))

    with open(re2_outfile, 'w+') as fp:
        fp.write('\t'.join(header) + '\n') 
    with open(fe_outfile, 'w+') as fp:
        fp.write('\t'.join(header) + '\n') 
    with open(joint_outfile, 'w+') as fp:
        fp.write('\t'.join(header) + '\n') 

    cols = ['VARIANT_ID', 'NSTUDY', 'PVALUE_FE', 'BETA_FE', 'STD_FE', 'PVALUE_RE', 'BETA_RE', 'STD_RE', \
    'PVALUE_RE2', 'STAT1_RE2', 'STAT2_RE2', 'PVALUE_BE', 'I_SQUARE', 'Q', 'PVALUE_Q', 'TAU_SQUARE']
    for i,chunk_df in enumerate(pd.read_csv(args.metasoft_file, sep='\t', chunksize=args.chunksize, skiprows=1, names=cols, usecols=cols)):

        print('\rProcessing chunk {}'.format(i+1), end='', flush=True)
    
        chunk_df['CHR'] = chunk_df['VARIANT_ID'].map(lambda x: str(x.split(',')[1].split(':')[0]))
        chunk_df['BP'] = chunk_df['VARIANT_ID'].map(lambda x: int(x.split(',')[1].split(':')[1]))
        chunk_df['GENE'] = chunk_df['VARIANT_ID'].map(lambda x: x.split(',')[0])

        for gene_id, grouped_df in chunk_df.groupby('GENE'):

            mask = (
                (abs(grouped_df['BP'] - gene2tss[gene_id]['tss']) <= args.cis_distance) &
                (grouped_df['CHR'] == gene2tss[gene_id]['chr'] )
            ).values
    
            assert(len(mask) == grouped_df.shape[0])

            grouped_df['TYPE'] = np.where(mask, 'CIS', 'TRANS')

            df = grouped_df[grouped_df['PVALUE_FE']<args.pthreshold]
            df[header].to_csv(fe_outfile, sep='\t', index=False, header=False, mode='a')

            df = grouped_df[grouped_df['PVALUE_RE2']<args.pthreshold]
            df[header].to_csv(re2_outfile, sep='\t', index=False, header=False, mode='a')

            df = grouped_df[(grouped_df['PVALUE_FE']<args.pthreshold) | (grouped_df['PVALUE_RE2']<args.pthreshold)]
            df[header].to_csv(joint_outfile, sep='\t', index=False, header=False, mode='a')

def main():

    genes = defaultdict(lambda: dict)
    with open(args.tss_bed) as fp:
        #fp.readline() # header
        for row in fp.read().strip().split('\n'): 
            genes[row.split('\t')[3]] = {'chr':str(row.split('\t')[0]), 'tss': int(row.split('\t')[1])}
        
    # prepare header
    # remove study p and m values

    header = ['VARIANT_ID', 'TYPE', 'NSTUDY', 'PVALUE_FE', 'BETA_FE', 'STD_FE', 'PVALUE_RE', 'BETA_RE', 'STD_RE', \
    'PVALUE_RE2', 'STAT1_RE2', 'STAT2_RE2', 'PVALUE_BE', 'I_SQUARE', 'Q', 'PVALUE_Q', 'TAU_SQUARE']

    stratify_by_distance(genes, args.metasoft_file, header)

    print('\nFinished post processing.')

if __name__ == '__main__':
    main()

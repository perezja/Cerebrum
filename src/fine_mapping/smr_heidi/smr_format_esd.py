import argparse
import pandas as pd
import numpy as np
import gzip
import os

parser = argparse.ArgumentParser(prog='Format METASOFT into .esd format for SMR.')
parser.add_argument('assoc_file', type=str, help='METASOFT association file.')
parser.add_argument('maf', type=str, help='MAF file from Plink (using "--freq").')
parser.add_argument('prefix', type=str, help='Outfile prefix.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')
parser.add_argument('-c', '--chunk_size', type=int, default=1000000, help='Chunksize for metasoft_file')
parser.add_argument('--debug', action='store_true', help='')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

def format_esd(assoc_df, maf_df):

    assoc_df = assoc_df[['RSID', 'BETA_RE', 'STD_RE', 'PVALUE_RE2']]

    rsid = assoc_df['RSID']
    assoc_df['SNP'] = rsid.map(lambda x: x.split(',')[1])
    assoc_df['gene'] = rsid.map(lambda x: x.split(',')[0])
    assoc_df['Chr'] = rsid.map(lambda x: x.split(',')[1].split(':')[0])
    assoc_df['Bp'] = rsid.map(lambda x: x.split(',')[1].split(':')[1])

    assoc_df = assoc_df.merge(maf_df, on='SNP') 

    assoc_df = assoc_df[['gene', 'Chr', 'SNP', 'Bp', 'A1', 'A2', 'MAF', 'BETA_RE', 'STD_RE', 'PVALUE_RE2']]
    assoc_df.columns = ['gene', 'Chr','SNP', 'Bp', 'A1', 'A2', 'Freq', 'Beta_re', 'se_re', 'pval_re2']

    return(assoc_df)

def main():

    smr_outfile = os.path.join(args.output_dir, args.prefix+'.txt.gz')
    header = ['gene', 'Chr','SNP', 'Bp', 'A1', 'A2', 'Freq', 'Beta_re', 'se_re', 'pval_re2']

    with gzip.open(smr_outfile, 'wt') as fp:
        fp.write('\t'.join(header) + '\n') 

    maf = pd.read_csv(args.maf, delim_whitespace=True)

    for i,chunk_df in enumerate(pd.read_csv(args.assoc_file, sep='\t', chunksize=args.chunk_size)):

        print('\rProcessing chunk {}'.format(i+1), end='', flush=True)

        esd_df = format_esd(chunk_df, maf)

        esd_df = esd_df.drop_duplicates()
        
        if esd_df.shape[0]==0:
            print('chunk '+str(i), ' has 0 rows.')

        esd_df.to_csv(smr_outfile, sep='\t', index=False, header=False, compression='gzip', mode='a')

if __name__ == '__main__':
    main()

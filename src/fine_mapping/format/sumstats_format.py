import argparse
import pandas as pd
import numpy as np
import gzip
import os

parser = argparse.ArgumentParser(prog='Format METASOFT into .esd format for SMR.')
parser.add_argument('assoc_file', type=str, help='METASOFT association file.')
parser.add_argument('maf', type=str, help='MAF file from Plink (using "--freq").')
parser.add_argument('ensembl2gene', type=str, help='MAF file from Plink (using "--freq").')
parser.add_argument('prefix', type=str, help='Outfile prefix.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')
parser.add_argument('-c', '--chunk_size', type=int, default=1000000, help='Chunksize for metasoft_file')
parser.add_argument('--debug', action='store_true', help='')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

def format_esd(assoc_df, maf_df):

    assoc_df = assoc_df[['VARIANT_ID', 'TYPE', 'BETA_FE', 'STD_FE', 'PVALUE_FE']]

    rsid = assoc_df['VARIANT_ID']
    assoc_df['SNP'] = rsid.map(lambda x: x.split(',')[1])
    assoc_df['gene'] = rsid.map(lambda x: x.split(',')[0])
    assoc_df['Chr'] = rsid.map(lambda x: x.split(',')[1].split(':')[0])
    assoc_df['Bp'] = rsid.map(lambda x: x.split(',')[1].split(':')[1])

    assoc_df = assoc_df.merge(maf_df, on='SNP') 

    assoc_df = assoc_df[['gene', 'Chr', 'SNP', 'TYPE', 'Bp', 'A1', 'A2', 'MAF', 'BETA_FE', 'STD_FE', 'PVALUE_FE']]

    return(assoc_df)

def main():

    smr_outfile = os.path.join(args.output_dir, args.prefix+'.esd')
    sumstats_outfile = os.path.join(args.output_dir, args.prefix+'.sumstats.txt')

    header = ['Chr','SNP','Bp','A1','A2','Freq','Beta','se','p']
    with open(smr_outfile, 'wt') as fp:
        fp.write('\t'.join(header) + '\n') 

    header = ['SNP','gene','name', 'type','A1', 'A2', 'maf', 'beta', 'se', 'p']
    with open(sumstats_outfile, 'wt') as fp:
        fp.write('\t'.join(header) + '\n') 

    with open(args.ensembl2gene) as fp:
        e2g = {line.split('\t')[0]: line.split('\t')[1] for line in fp.read().strip().split('\n')}

    maf = pd.read_csv(args.maf, delim_whitespace=True)

    for i,chunk_df in enumerate(pd.read_csv(args.assoc_file, sep='\t', chunksize=args.chunk_size)):

        print('\rProcessing chunk {}'.format(i+1), end='', flush=True)

        assoc_df = format_esd(chunk_df, maf)
        assoc_df['name'] = assoc_df['gene'].map(lambda x: e2g[x])

        format_df = assoc_df[['SNP', 'gene', 'name', 'TYPE', 'A1', 'A2', 'MAF', 'BETA_FE', 'STD_FE', 'PVALUE_FE']]
        format_df.to_csv(sumstats_outfile, sep='\t', index=False, header=False, mode='a')

        format_df = assoc_df[['Chr', 'SNP', 'Bp', 'A1', 'A2', 'MAF', 'BETA_FE', 'STD_FE', 'PVALUE_FE']]
        format_df.to_csv(smr_outfile, sep='\t', index=False, header=False, mode='a')

if __name__ == '__main__':
    main()

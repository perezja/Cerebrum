#!/usr/bin/env python3
import argparse
import subprocess
import pandas as pd
import tempfile
import shutil
import numpy as np
import re
import os

def parse_args():

    parser = argparse.ArgumentParser(description='Prepare GWAS loci for colocalization analysis.')
    parser.add_argument('gwas_file', help='File to summary statistics eQTL file')
    parser.add_argument('gene_bed', help='BED of grch37 gene TSS coordinates.')
    parser.add_argument('sumstats_file', help='Genome-wide significant eQTLs in summary statistics.')
    parser.add_argument('prefix', help='Prefix for output file.')
    parser.add_argument('--gene_mb', default=1, type=int, help='Distance of eGenes from SNP.')
    parser.add_argument('--eqtl_kb', default=100, type=int, help='Distance of eQTL to lead GWAS SNP.')
    parser.add_argument('-o', '--output_dir', help='Output directory.')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    args.output_dir = os.path.abspath(args.output_dir)

    return(args)
    
def intra_rep_peak_merge(narrow_peak, tmpdir):

    # keep peak name and mean_signal (cols=4,7) 

    prefix = os.path.split(narrow_peak)[1].split('_')[0]
    outfile = os.path.join(tmpdir, prefix + 'intra_rep_merge.narrowPeak')

    cmd='tail -n+2 ' + narrow_peak \
      + ' | sort -k1,1 -k2,2n' \
      + ' | bedtools merge -i -' \
      + ' -d ' + str(args.window_size) \
      + ' -c 4,5,7' \
      + ' -o collapse,mean,mean' \
      + ' -delim ' + '","' \
      + ' > ' + outfile 

    if debug:
        print('* '+cmd)

    subprocess.check_call(cmd, shell=True)

def run_bedtools(gwas_bed, gene_bed, outfile):

    cmd= 'bedtools window -a {}'.format(gwas_bed) \
      + ' -b ' + gene_bed \
      + ' -w ' + str(args.gene_mb*1e6) \
      + ' > ' + outfile

    subprocess.check_call(cmd, shell=True)

    df = pd.read_csv(outfile, sep='\t', header=None, names=['Chr', 'Bp', 'stop', 'SNP', 'gene_chr', 'tss', 'tssp1', 'gene', 'strand']) 

    df = df[['SNP', 'gene']]

    return(df)

def main(args):

    # 1. Find genes +/- 1Mb from GWAS SNP

    gwas_df = pd.read_csv(args.gwas_file, sep='\t')
    gwas_df['stop'] = gwas_df['Bp'] + 1
    gwas_df = gwas_df[['Chr', 'Bp', 'stop', 'SNP']]

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    infile = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
    outfile = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)

    gwas_df.to_csv(infile.name, sep='\t', index=False, header=False)
    gg_df = run_bedtools(infile.name, args.gene_bed, outfile.name)

    # 2. Find GWAS SNPs with >= 1 eGene in locus

    eqtl_df = pd.read_csv(args.sumstats_file, sep='\t')
    eqtl_df = eqtl_df[eqtl_df['type']=='CIS']

    eqtl_df = eqtl_df[eqtl_df['gene'].isin(gg_df['gene'])]

    print('{} total eGenes in GWAS loci'.format(len(eqtl_df['gene'].unique())))

    total_loci = len(gg_df['SNP'].unique())
    gg_df = gg_df[gg_df['gene'].isin(eqtl_df['gene'])]

    print('Amounting to {}/{} total loci'.format(len(gg_df['SNP'].unique()),total_loci))

    # 3. Create loci for GWAS SNPs with eGene

    gl_df = pd.DataFrame()
    for gwas_snp in gg_df['SNP'].unique():
        pos = int(gwas_snp.split(':')[1])
        gl_df = gl_df.append(pd.Series({'CHR':gwas_snp.split(':')[0], 'SNP':gwas_snp, 'BP':pos, 'STARTBP':(pos-args.eqtl_kb*1e3), 'ENDBP':(pos+args.eqtl_kb*1e3)}), ignore_index=True) 


    gl_df = gl_df[['CHR', 'SNP', 'BP', 'STARTBP', 'ENDBP']]
    gl_df = gl_df.astype({'BP':int, 'STARTBP':int, 'ENDBP':int})

    outfile = os.path.join(args.output_dir, args.prefix+'.index')
    gl_df.to_csv(outfile, sep='\t', index=False)

    print('wrote index SNP interval file to: {}'.format(outfile)) 
    
    shutil.rmtree(tmpdir)

if __name__ == '__main__':
    args = parse_args()
    main(args)

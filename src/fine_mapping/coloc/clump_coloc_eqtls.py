#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import subprocess
from collections import defaultdict
import tempfile
import glob
import os
import re

parser = argparse.ArgumentParser(prog='.')
parser.add_argument('bfile', type=str, help='METASOFT file with all cis associations')
parser.add_argument('coloc_res_file', type=str, help='GWAS Index SNP file.')
parser.add_argument('locus_dir', type=str, help='Directory with locus association inputs to coloc.')
parser.add_argument('eqtl_sumstats', type=str, help='')
parser.add_argument('-p', '--coloc_threshold', default=0.7, type=float, help='PP.H4 threshold for calling a colocalized locus.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')
parser.add_argument('-c', '--chunksize', type=int, default=1e6, help=".")

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

def run_ld_clumping(bfile, assoc_file, prefix, output_dir, p1='5e-8', p2='5e-8', pval_field='p', snp_field='SNP', kb='100', r2='0.8'):

    outdir = os.path.join(output_dir, 'plink_clumping')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outprefix = os.path.join(outdir, prefix)

    cmd = 'plink1.9 --bfile {} --clump {}'.format(bfile, assoc_file) + \
          ' --clump-p1 {} --clump-p2 {}'.format(p1, p2) + \
          ' --clump-field {} --clump-snp-field {}'.format(pval_field, snp_field) + \
          ' --clump-kb {} --clump-r2 {}'.format(kb, r2) + \
          ' --out {}'.format(outprefix)

    subprocess.check_call(cmd, shell=True)

    return(outprefix+'.clumped')

def clump_loci(locus_labels, eqtl_df, cr_df, tmpdir):

    locus_files = [(locus,os.path.join(args.locus_dir, locus)) for locus in locus_labels]

    coloc_df = pd.DataFrame()

    for locus,locus_file in locus_files:

        locus_df = pd.read_csv(locus_file, sep='\t') 
        locus_eqtl_df = eqtl_df[eqtl_df['SNP'].isin(locus_df['SNP'])]
        plinkin = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)

        locus_eqtl_df.to_csv(plinkin.name, sep='\t', index=False)

        clumped_file = run_ld_clumping(args.bfile, plinkin.name, 'coloc_eqtls', tmpdir)
        cl_df = pd.read_csv(clumped_file, delim_whitespace=True)

        top_locus_snp = cl_df.iloc[1,:] 
        coloc_snps = [re.split(r'\([0-9+]\)$', snp)[0] for snp in top_locus_snp.SP2.strip().split(',')]
        coloc_snps = coloc_snps + [top_locus_snp.SNP]

        eqtls = eqtl_df[eqtl_df['SNP'].isin(coloc_snps)]

        coloc_snps = [':'.join(i.split(':')[:2]) for i in (coloc_snps)]

        locus_desc = cr_df[cr_df['locus_name']==locus]
        locus_desc.reset_index(inplace=True)

        eqtls = eqtls.assign(locus_name=locus_desc['locus_name'][0]) 
        eqtls = eqtls.assign(trait=locus_desc['trait'][0]) 
        eqtls = eqtls.assign(PPH4=locus_desc['PP.H4.abf'][0]) 

        coloc_df = coloc_df.append(eqtls)
 
    return(coloc_df.dropna()) 

def aggregate_loci(locus_dir, locus_labels):
    
    locus_list = list()
    for fp in glob.glob(os.path.join(locus_dir,'cortex_region*')):
        if os.path.split(fp)[1] in locus_labels:
            locus_list.append(pd.read_csv(fp, sep='\t'))    
    print('Aggregated association data from {} loci.'.format(len(locus_list)))
        
    df = pd.concat(locus_list)
    return(df.dropna()) 

def main():
    
    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    eqtl_df = pd.read_csv(args.eqtl_sumstats, sep='\t')
    cr_df = pd.read_csv(args.coloc_res_file, sep='\t')

    cr_df = cr_df[cr_df['PP.H4.abf']>args.coloc_threshold]
    print('{} colocalised loci detected.'.format(cr_df.shape[0]))

    coloc_loci = cr_df['locus_name'].tolist()
    coloc_df = clump_loci(coloc_loci, eqtl_df, cr_df, tmpdir)
    
    outfile = os.path.join(args.output_dir, 'neurodeg_coloc_clumped.txt')
    coloc_df.to_csv(outfile, sep='\t', index=False)

    print('wrote to: {}'.format(outfile))

if __name__ == '__main__':
    main()

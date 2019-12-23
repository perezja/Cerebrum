import argparse
import pandas as pd
import numpy as np
import gzip
import os

parser = argparse.ArgumentParser(prog='Format GWAS data for COLOC and MVMR.')
parser.add_argument('gwas_file', type=str, help='METASOFT association file.')
parser.add_argument('study', type=str, choices=['kunkle', 'metal_pd', 'rheenen'], help='GWAS study name.')
parser.add_argument('prefix', type=str, help='Outfile prefix.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')

args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

args = parser.parse_args()

df = pd.read_csv(args.gwas_file, delim_whitespace=True)

gwas_outfile = os.path.join(args.output_dir, args.prefix+'.outcome.txt')

if args.study == 'kunkle':

    df['SNP'] = df.apply(lambda x: str(x['Chromosome']) + ':' + str(x['Position']) + ':' + x['Effect_allele'] + ':' + x['Non_Effect_allele'], axis=1)
    df = df.rename({'Chromosome':'Chr', 'Position':'Bp'}, axis=1)
    df.drop('MarkerName', inplace=True, axis=1)

    df['Proportion_cases'] = 21982/41944 
    df['Nsamples']=21982+41944
#    df['Proportion_cases'] = 21982/(41944+21982)
    df = df[['Chr', 'Bp', 'SNP', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue', 'Proportion_cases', 'Nsamples']]

    df.to_csv(gwas_outfile, sep='\t', header=True, index=False)


elif args.study == 'metal_pd':

    df['Proportion_cases'] = df.apply(lambda x: x['TotalN_CASE']/x['TotalN_CONTROL'], axis=1)
    df['Nsamples']=df.apply(lambda x: x['TotalN_CASE']+x['TotalN_CONTROL'], axis=1)
#    df['Proportion_cases'] = df.apply(lambda x: x['TotalN_CASE']/(x['TotalN_CONTROL']+x['TotalN_CASE']), axis=1)
    df = df[df['Proportion_cases'] < 1]
    df = df.rename({'ANNO_CHR':'Chr', 'ANNO_BP':'Bp', 'Allele1':'Effect_allele', 'Allele2':'Non_Effect_allele', 'Effect':'Beta', 'StdErr':'SE', 'P-value':'Pvalue'}, axis=1) 
    df['SNP'] = df.apply(lambda x: str(x['Chr']) + ':' + str(x['Bp']) + ':' + x['Effect_allele'] + ':' + x['Non_Effect_allele'], axis=1) 
    df = df[['Chr', 'Bp', 'SNP', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue', 'Proportion_cases', 'Nsamples']]
    df.to_csv(gwas_outfile, sep='\t', header=True, index=False)

else:
    # rheenan

    df = df.dropna()
    df['snp'] = df.apply(lambda x: ':'.join([str(x['chr']), str(x['bp']), str(x['a1']), str(x['a2'])]), axis=1)
    
    Ncases=12577; Ncontrols=23475
    #Pcase = Ncases / Ncontrols
    Pcase= Ncases/(Ncases+Ncontrols)
    df['Nsamples']=Ncases+Ncontrols


    # need to convert LMM beta to normal OR
    #OR=df['b'].map(lambda x: ((Pcase + x)/(1 - Pcase -  x)) / (Pcase/(1 - Pcase)))

    # the same approximation gives you the SE

    #SE=df['se'].map(lambda x: ((Pcase + x)/(1 - Pcase -  x)) / (Pcase/(1 - Pcase)))

    df['Proportion_cases'] = Pcase 
    #df['b'] = OR
    #df['se'] = SE

    df = df.rename({'chr':'Chr', 'bp':'Bp', 'snp':'SNP', 'a1':'Effect_allele', 'a2':'Non_Effect_allele', 'b':'Beta', 'se':'SE', 'p':'Pvalue'}, axis=1)
    df = df[['Chr', 'Bp', 'SNP', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue', 'Proportion_cases', 'Nsamples']]
    df.to_csv(gwas_outfile, sep='\t', header=True, index=False)


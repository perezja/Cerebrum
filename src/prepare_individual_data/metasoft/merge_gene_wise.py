#!/usr/bin/env python3
# Author: James A. Perez 

import argparse
import numpy as np
import pandas as pd
from datetime import datetime
import subprocess
import string
import json
import random
import os
import gzip

parser = argparse.ArgumentParser(description='Prepare input for Metasoft.')
parser.add_argument('gene_list', help='File listing genes to include in Metasoft input chunk.')
parser.add_argument('gene_json', help='JSON file containing path "has_gene" attribute.')
parser.add_argument('-o', '--output_dir', default='./', help='Directory for output files.')
args = parser.parse_args() 

def main():

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    with open(args.gene_json) as fp:
        gene_json = json.load(fp)
    with open(args.gene_list) as fp:
        gene_list = fp.read().strip().split('\n') 

    def get_variant_ids(gene):

        variant_ids = set()
        for rid,fp in gene_json[gene].items():

            cmd = "gunzip -c {} | awk -F\"\\t\" '{{print $2\",\"$1 }}' - ".format(fp)
            pids = subprocess.check_output(cmd, shell=True).decode('utf8').strip().split('\n')

            variant_ids.update(pids)

        return(list(variant_ids))


    metasoft_df = pd.DataFrame()

    region_ids = list(gene_json[list(gene_json.keys())[0]].keys())
    region_ids.sort()

    for n, gene in enumerate(gene_list):

        variant_ids = get_variant_ids(gene)

        ms_df = pd.DataFrame(np.nan, index=variant_ids,  columns=[j for i in region_ids for j in [i+'_slope', i+'_slope_se']], dtype=np.float64)

        for rid in region_ids:

            rg_df = pd.read_csv(gene_json[gene][rid], sep='\t', header=None, names=['SNP','gene','slope','t-stat', 'p-value', 'FDR'], usecols=['SNP','gene','slope','t-stat','p-value'])
            # remove assocations with zero effect size (zero t-statistic)
            rg_df = rg_df[rg_df['t-stat']>0]

            rg_df.index = rg_df['gene']+','+rg_df['SNP']
            rg_df['slope_se'] = rg_df.apply(lambda x: x['slope'] / x['t-stat'], axis=1)

            ms_df[[rid+'_slope', rid+'_slope_se']] = ms_df.join(rg_df, how='outer')[['slope','slope_se']]

        if metasoft_df.empty:
            metasoft_df = ms_df
        else:
            metasoft_df = metasoft_df.append(ms_df)

        print("[ {} ] * gene '{}' ({}/{}) processed.".format(datetime.now().strftime("%b %d %H:%M:%S"), gene, str(n+1), str(len(gene_list))))

    metasoft_df.index.name = 'variant_id'

    def random_string(l):
        letters = string.ascii_letters
        return ''.join(random.choice(letters) for i in range(l))

    unique_id = random_string(6) 
    ms_outfile = os.path.join(args.output_dir, 'ms_chunk_'+unique_id+'.txt.gz')

    if not metasoft_df.empty:
        with gzip.open(ms_outfile, 'wt', compresslevel=1) as f:
            metasoft_df.to_csv(f, sep='\t', na_rep='NA')

if __name__ == "__main__":
    main()


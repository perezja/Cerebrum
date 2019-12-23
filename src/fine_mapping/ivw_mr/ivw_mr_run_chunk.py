#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import numpy as np
import tempfile
import shutil
import queue
import json
import time
import logging
from common import *
from collections import defaultdict
import os
import re
import sys


def parse_args():

    parser = argparse.ArgumentParser(description='Get pruned list of variants to keep.')
    parser.add_argument('bfile', help='Plink bfile containing variant data for all participants.')
    parser.add_argument('exposure_data', help='Chunk of exposure association data.') 
    parser.add_argument('outcome_data', help='Outcome data for TwoSampleMR.')
    parser.add_argument('outcome_id', help='Name of outcome measure.')
    parser.add_argument('--prefix', default='out', help='') 
    parser.add_argument('--subset', help='Subset of exposures.')
    parser.add_argument('--chunksize', type=int, default=50, help='Subset of exposures.')
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET','DEBUG','INFO', 'WARNING','CRITICAL',
                                                        'ERROR','CRITICAL'], help='Log level')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    args.output_dir = os.path.abspath(args.output_dir)

    log.setLevel(args.log_level)
    log.info(sys.argv)    

    return(args)

def main(args):
         
    harm_outdir = os.path.join(args.output_dir, 'harmonised_data')
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if not os.path.exists(harm_outdir):
        os.makedirs(harm_outdir)

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)

    # Harmonise all exposures with outcome dataframe

    harmonised_file = harmonise_snps(args.exposure_data, args.outcome_data, args.outcome_id, args.prefix, harm_outdir) 

    # Run LD clumping on harmonised SNPs

    clumped_dat = run_ld_clumping(args.bfile, harmonised_file, args.prefix, tmpdir)

    # Subset harmonised dataframe with independent index SNPs

    log.info('loading independent SNPs...')
    with open(clumped_dat) as fp:
            # header
        fp.readline()
            # SNP = col3
        indep_snps = [line.strip().split()[2] for line in fp.read().strip().split('\n')]

    dat = pd.read_csv(harmonised_file, sep='\t') 
    n_pre = dat.shape[0]

    dat = dat[dat['SNP_ID'].isin(indep_snps)] 
    log.info('{} SNPs filtered. {} remaining.'.format(n_pre, dat.shape[0]))

    ivw_dat_outfile = os.path.join(tmpdir, args.prefix + '_ivw_input.txt')
    dat.to_csv(ivw_dat_outfile, sep='\t', index=False) 

    # Run IVR MR

    arg_list = [ivw_dat_outfile, args.prefix, '--sensitivity', '-o', args.output_dir] 
    run_twosamplemr(arg_list)

    log.info('Finished TwoSampleMR analysis.')

    shutil.rmtree(tmpdir)
    os.remove(args.exposure_data)

if __name__ == '__main__':
    args = parse_args()
    main(args)

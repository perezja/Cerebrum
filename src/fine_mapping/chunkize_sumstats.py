#!/usr/bin/env python3
# Author: James A. Perez 

import argparse
import subprocess
import pandas as pd
import numpy as np
import os
import json

def chunkize_by_genes():
    """ Fragments eqtl association data into gene-level chunks.
    """
    if args.sumstats_file.endswith('.gz'):
        cmd = "zcat{0} | "
    else:
        cmd = ""
    
    of_prefix = os.path.join(args.output_dir, args.prefix+'.')
    cmd += "awk -F\"\\t\" 'BEGIN{{OFS=\"\\t\";}} NR==1{{hdr=$0;next }} !($2 in a) {{ print hdr > \"{1}\"$2\".sumstats\"; a[$2]++}} NR>1 {{ gene=$2; print $0  > \"{1}\"gene\".sumstats\"}}' {0} ".format(args.sumstats_file, of_prefix)
    subprocess.check_call(cmd, shell=True)
   
parser = argparse.ArgumentParser(description='Chunkize a summary statistics eQTL file for downstream processing.')
parser.add_argument('sumstats_file', help='Path to eQTL sumstats file.')
parser.add_argument('prefix', help='Prefix for chunk files.')
parser.add_argument('-o', '--output_dir', default='chunks', help='Directory for output chunk files.')
args = parser.parse_args() 

args.output_dir = os.path.abspath(args.output_dir)

def main():

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    chunkize_by_genes()

    print('wrote chunks to: {}'.format(args.output_dir))

if __name__ == "__main__":
    main()


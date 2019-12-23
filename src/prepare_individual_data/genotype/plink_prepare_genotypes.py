#!/bin/env python3
import argparse
import subprocess
from datetime import datetime
import os

parser = argparse.ArgumentParser(description='Prepare genotypes for QTL analysis.')
parser.add_argument('bfile', help='Plink bfile containing variant data for all participants.')
parser.add_argument('region_participant_ids', help='Directory with <region>_id_list.txt files to subset the bfile.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

def id(x):
    return(os.path.split(x)[1].split('_')[0])
def fp(x):
    return(os.path.join(args.region_participant_ids, x))

region_id_list = [ (id(i), fp(i)) for i in os.listdir(args.region_participant_ids) if i.endswith('.txt')] 

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

for rid, fn in region_id_list:

      
    print("[ {} ] Processing '{}'".format(datetime.now().strftime('%b %d %H:%M:%S'), rid)) 

    prefix = os.path.join(args.output_dir, rid)
    cmd = 'plink2 --bfile ' + args.bfile \
        + ' --recode A-transpose' \
        + ' --keep-fam '+ fn \
        + ' --out ' + prefix 

    subprocess.check_call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    cmd = "cut -f2,7- " + prefix + ".traw" + " > " + prefix + ".snps.txt" 

    subprocess.check_call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    os.remove(prefix+".traw")


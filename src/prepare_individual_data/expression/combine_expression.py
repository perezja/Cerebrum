#!/usr/bin/env python3
import numpy as np
import pandas as pd
import subprocess
import argparse
from datetime import datetime
import simplejson as json
import tempfile
import shutil
import glob
import gzip
import re
import os

class CombineExpression():
    def __init__(self, tpm_counts_json):
        with open(tpm_counts_json, 'r') as f:
            GE = json.load(f)

        self.CTS = dict()
        self.TPM = dict()

        for study, regions in GE.items():
            for rid, paths in regions.items():

                print("* loading paths for '{}' in '{}'.".format(rid, study))
                assert(all([os.path.exists(i) for i in list(paths.values())]))

                self.CTS[rid] = self.__get_counts_files(paths.get("counts"))
                self.TPM[rid] = self.__get_tpm_files(paths.get("tpm"))

        self.__tmp = tempfile.mkdtemp(dir=args.temp_dir)
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)

    def __get_counts_files(self, path):
        cts = glob.glob(os.path.join(path,'*','*.ReadsPerGene.out.tab'))
        cts = pd.Series(cts, index=[i.split('/')[-2] for i in cts])
        cts.sort_index(inplace=True)

        return(cts)
    def __get_tpm_files(self, path):
        tpm = glob.glob(os.path.join(path,'*','quant.sf'))
        tpm = pd.Series(tpm, index=[i.split('/')[-2] for i in tpm])
        tpm.sort_index(inplace=True)

        return(tpm)

    def __merge_region_counts(self, rid, paths):

        skip=lambda x: x in range(4) # header from STAR gct
        usecols=[0,1] # cols={gene_id, all_counts}
        sample_ids = paths.index.values

        df = pd.read_csv(paths[0], sep='\t', skiprows=skip, header=None, usecols=usecols, names=['gene_id', sample_ids[0]], index_col=0)
        df.index = df.index.str.replace(r'\.[0-9]+','')

        if df[sample_ids[0]].dtype == np.float64:
            dtype = np.float32
        elif df[sample_ids[0]].dtype == np.int64:
            dtype = np.int32
        else:
            dtype = df[sample_ids[0]].dtype

        gct_df = pd.DataFrame(0, index=df.index, columns=list(sample_ids), dtype=dtype)
        gct_df[sample_ids[0]] = df[sample_ids[0]].astype(dtype)
        for k, (i,p) in enumerate(zip(sample_ids[1:], paths[1:])):

            print("\rProcessing '{}': {}/{}".format(rid, (k+2), len(paths)), end='', flush=True)
            df = pd.read_csv(p, sep='\t', skiprows=skip, header=None, usecols=usecols, names=['gene_id', i], index_col=0)
            df.index = df.index.str.replace(r'\.[0-9]+','')
            gct_df[i] = df[i]
            
        print()

        with gzip.open(os.path.join(args.output_dir, rid+'.gct.gz'), 'wt', compresslevel=6) as f:
            f.write('{0}\t{1}\n'.format(gct_df.shape[0], gct_df.shape[1]))
            gct_df.to_csv(f, sep='\t', float_format='%.6g')

    def __merge_region_tpm(self, rid, paths):

        tpms = os.path.join(self.__tmp, rid+'_paths.txt')
        sample_ids = os.path.join(self.__tmp, rid+'_ids.txt')

        with open(tpms, 'w+') as f:
            f.write('\n'.join([i for i in paths]))
        with open(sample_ids, 'w+') as f:
            f.write('\n'.join([i for i in paths.index.values]))

        exc = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'merge_tpm.R')
        outfile = os.path.join(args.output_dir, rid+'.gtt.gz')
        cmd = ' '.join([exc, tpms, sample_ids, args.gtf, '-o', outfile]) 

        subprocess.check_call(cmd, shell=True)

    def merge_region_counts(self):

        print("[ {} ] Merging Counts.".format(datetime.now().strftime("%b %d %H:%M:%S")))

        for rid, paths in self.CTS.items():
            self.__merge_region_counts(rid, paths)

    def merge_region_tpm(self):

        print("[ {} ] Merging TPM.".format(datetime.now().strftime("%b %d %H:%M:%S")))

        for k, (i, p) in enumerate(self.TPM.items()):

            print("\rProcessing '{}': {}/{}".format(i, k+1, len(self.TPM.keys())),end='',flush=True)
            self.__merge_region_tpm(i, p)

    def __del__(self):
        if os.path.exists(self.__tmp):
            shutil.rmtree(self.__tmp)

parser = argparse.ArgumentParser(description='Run pipeline from RNA-Seq JSON file')
parser.add_argument('json', type=str, help='Path to the JSON file')
parser.add_argument('gtf', type=str, help='Path to the gtf file')
parser.add_argument('--mode', required=True, choices=['cts','tpm'], type=str, help='Path to the gtf file')
parser.add_argument('--temp_dir',type=str, default='./',help='Temporary directory')

parser.add_argument('-o','--output_dir',type=str, default='.',help='File with sample subset to process')

args = parser.parse_args()


def main():
    ce = CombineExpression(args.json)
    if args.mode=='cts':
        ce.merge_region_counts()
    else:
        ce.merge_region_tpm()

if __name__ == "__main__":
    main()

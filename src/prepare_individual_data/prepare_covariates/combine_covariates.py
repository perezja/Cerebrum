#!/usr/bin/env python3
# Author: Francois Aguet
import pandas as pd
import numpy as np
import argparse
import subprocess
import os
from datetime import datetime
import contextlib

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

class PrepareCovariates():
    def __init__(self, region_participant_lookup, joint_phenotype_file):

        with open(region_participant_lookup) as f:
            rp_list = f.read().strip().split('\n')

        rp_list = [(rp.split('\t')[0], rp.split('\t')[1]) for rp in rp_list]

        self.region_ids = list(set([rp[1] for rp in rp_list])) 
        self.region_ids.sort()

        print('\nRead in participants for {} regions.'.format(len(self.region_ids)))

        self.RP = dict()
        for rid in self.region_ids:
            self.RP[rid] = [rp[0] for rp in rp_list if rp[1] == rid]

        self.joint_pheno_df = pd.read_csv(joint_phenotype_file, sep='\t', index_col=0)  

        self.expression_covariates = dict() 

    def run_peer(self, rid):
        
        exc = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'run_PEER.R')
        outdir = os.path.join(args.output_dir, 'peer', rid)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        print("[ {} ] Calculating PEER factors for '{}'".format(datetime.now().strftime('%b %d %H:%M:%S'), rid))
        rg_exp_file = os.path.join(args.exp_dir, rid+".expression.bed.gz") 

        with cd(outdir):

            subprocess.check_call(' '.join([exc, rg_exp_file, rid, args.n]), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.expression_covariates[rid] = os.path.join(outdir, rid+'.PEER_covariates.txt')

    def combine_covariates(self):

        for rid in self.region_ids:

            self.run_peer(rid)

            print('Combining phenotype and expression covariates...')

            mask = self.joint_pheno_df.index.isin(self.RP[rid])
            phenotypes_df = self.joint_pheno_df[mask]

            assert(len(self.RP[rid])==phenotypes_df.shape[0])

            phenotypes_df = phenotypes_df.transpose()

            expression_df = pd.read_csv(self.expression_covariates[rid], sep='\t', index_col=0, dtype=str)
            combined_df = pd.concat([phenotypes_df[expression_df.columns], expression_df], axis=0)
            combined_df.index.name = "ID"

            # from Francois Aguet (1) 
            # identify and drop colinear covariates

            C = combined_df.astype(np.float64).T
            Q,R = np.linalg.qr(C-np.mean(C, axis=0))
            colinear_ix = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
            if np.any(colinear_ix):
                print('Colinear covariates detected:')
                for i in C.columns[colinear_ix]:
                    print("  * dropped '{}'".format(i))
                    combined_df = combined_df.loc[~colinear_ix]

            combined_df.to_csv(os.path.join(args.output_dir, rid+'.combined_covariates.txt'), sep='\t')

        print('done.')


parser = argparse.ArgumentParser(description='Combine covariates for region-specific participants into a single matrix.')
parser.add_argument('joint_phenotype_file', help='Full path to joint phenotype file for all participants across regions.')
parser.add_argument('region_participant_lookup', help='Full path to lookup table linking all participants to a given region.')
parser.add_argument('exp_dir', help='Directory with final expression tables.')
parser.add_argument('-n', type=str, default='15', help='Number of hidden confounders to include.')
parser.add_argument('-o', '--output_dir', default='combined_covariates', help='Output directory')
args = parser.parse_args()
 
def main():
    
    pc = PrepareCovariates(args.region_participant_lookup, args.joint_phenotype_file)
    pc.combine_covariates()

if __name__ == "__main__":
    main()


# (1) https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/combine_covariates.py

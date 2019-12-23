#!/usr/bin/env python3
# Author: James A. Perez 
# (adapted from Francois Aguet)

import numpy as np
import pandas as pd
import gzip
import subprocess
import scipy.stats as stats
from datetime import datetime
import argparse
import json
import os
import re

import rnaseqnorm

class NormalizeExpression():
    def __init__(self, gct_dir, gtt_dir, sample_participant_lookup, sample_id_list=None):

        with open(sample_participant_lookup) as f:
            sp_list = f.read().strip().split('\n')

        """ ['sample_id'] = 'participant_id' """
        self.SP = { sp.split('\t')[0]: sp.split('\t')[1] for sp in sp_list}

        if sample_id_list:
            with open(sample_id_list, 'r') as f:
                sample_id_list = f.read().strip().split('\n')

            self.SP = { sid: self.SP[sid] for sid in sample_id_list }

        gct_files = [os.path.join(gct_dir, i) for i in os.listdir(gct_dir) ]
        gct_files.sort(key=lambda x: os.path.split(x)[1].split('.')[0])
        gtt_files = [os.path.join(gtt_dir, i) for i in os.listdir(gtt_dir) ]
        gtt_files.sort(key=lambda x: os.path.split(x)[1].split('.')[0])

        self.regions = [os.path.split(i)[1].split('.')[0] for i in gct_files]

        self.GE_files = dict()
        for (rid, gct, gtt) in zip(self.regions, gct_files, gtt_files):
            self.GE_files[rid] = {"gct": gct, "gtt":gtt}

        self.excluded_samples = dict()

        """ container to hold expression data during processes"""
        self.counts_df = None
        self.tpm_df = None
        self.norm_df = None
        self.bed_df = None

        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)

    def __read_expression(self, rid):

        (gct_file, gtt_file) = self.GE_files[rid].values()

        assert(os.path.split(gct_file)[1].split('.')[-2] == 'gct')
        assert(os.path.split(gtt_file)[1].split('.')[-2] == 'gtt')

        print("[ {} ] Reading expression for '{}'.".format(datetime.now().strftime('%b %d %H:%M:%S'), rid))
        counts_df = pd.read_csv(gct_file, sep='\t', skiprows=1, compression='gzip', index_col=0)
        tpm_df = pd.read_csv(gct_file, sep='\t', skiprows=1, compression='gzip', index_col=0)

        ix = np.intersect1d(counts_df.index, tpm_df.index) 
        counts_df = counts_df.loc[ix]
        tpm_df = tpm_df.loc[ix]

        assert(counts_df.index.tolist() == tpm_df.index.tolist())

        mask = counts_df.columns.isin(self.SP) 

        self.excluded_samples[rid] = list(counts_df.columns[[not i for i in mask]])

        print("  * ({}/{}) matching samples in sample_to_participants_lookup.".format(np.sum(mask), counts_df.shape[1]))

        counts_df = counts_df.loc[:, mask]
        tpm_df = tpm_df.loc[:, mask]

        assert(counts_df.columns.tolist() == tpm_df.columns.tolist())

        self.counts_df = counts_df
        self.tpm_df = tpm_df
        
        return(0)

    def __gtf_to_bed(self, annotation_gtf, feature='gene', exclude_chrs=[]):
        """
        Parse genes from GTF, create placeholder DataFrame for BED output
        """
        chrom = []
        start = []
        end = []
        gene_id = []
        with open(annotation_gtf, 'r') as gtf:
            for row in gtf:
                row = row.strip().split('\t')
                if row[0][0]=='#' or row[2]!=feature: continue # skip header
                chrom.append(row[0])
    
                # TSS: gene start (0-based coordinates for BED)
                if row[6]=='+':
                    start.append(np.int64(row[3])-1)
                    end.append(np.int64(row[3]))
                elif row[6]=='-':
                    start.append(np.int64(row[4])-1)  # last base of gene
                    end.append(np.int64(row[4]))
                else:
                    raise ValueError('Strand not specified.')
    
                gene_id_elem = row[8].split(';',1)[0].split(' ')[1].replace('"','')
                gene_id_elem = re.sub(r'\.[0-9]+','',gene_id_elem) 
                gene_id.append(gene_id_elem)
    
        bed_df = pd.DataFrame(data={'chr':chrom, 'start':start, 'end':end, 'gene_id':gene_id}, columns=['chr', 'start', 'end', 'gene_id'], index=gene_id)
        # drop rows corresponding to excluded chromosomes
        mask = np.ones(len(chrom), dtype=bool)
        for k in exclude_chrs:
            mask = mask & (bed_df['chr']!=k)
        bed_df = bed_df[mask]

        self.bed_df = bed_df

        return(0)

   
    def __prepare_bed(self, chr_subset=None, ignoreVersion=True):

        bed_df = pd.merge(self.bed_df, self.norm_df, left_index=True, right_index=True)
    
        # sort by start position
        bed_df = bed_df.groupby('chr', sort=False, group_keys=False).apply(lambda x: x.sort_values('start'))
    
        if chr_subset is not None:
            # subset chrs from VCF
            bed_df = bed_df[bed_df.chr.isin(chr_subset)]

        self.bed_df = bed_df

        return(0)
    
    def __write_bed(self, output_name):
        """
        Write DataFrame to BED format
        """
        assert self.bed_df.columns[0]=='chr' and self.bed_df.columns[1]=='start' and self.bed_df.columns[2]=='end'
        # header must be commented in BED format
        header = self.bed_df.columns.values.copy()
        header[0] = '#chr'
        self.bed_df.to_csv(output_name, sep='\t', index=False, header=header)
        subprocess.check_call('bgzip -f '+output_name, shell=True, executable='/bin/bash')
        subprocess.check_call('tabix -f '+output_name+'.gz', shell=True, executable='/bin/bash')

        return(0)

    def __apply_normalization(self, sample_frac_threshold=0.2, count_threshold=6, tpm_threshold=0.1, mode='tmm'):
        """
        Genes are thresholded based on the following expression rules:
          TPM >= tpm_threshold in >= sample_frac_threshold*samples
          read counts >= count_threshold in sample_frac_threshold*samples
        
        vcf_lookup: lookup table mapping sample IDs to VCF IDs
        
        Between-sample normalization modes:
          tmm: TMM from edgeR
          qn:  quantile normalization
        """
    
        ns = self.counts_df.shape[1]
    
        # expression thresholds
        mask = (
            (np.sum(self.tpm_df>=tpm_threshold,axis=1)>=sample_frac_threshold*ns) &
            (np.sum(self.counts_df>=count_threshold,axis=1)>=sample_frac_threshold*ns)
        ).values
    
        # apply normalization
        if mode.lower()=='tmm':
            tmm_counts_df = rnaseqnorm.edgeR_cpm(self.counts_df, normalized_lib_sizes=True)
            self.norm_df = rnaseqnorm.inverse_normal_transform(tmm_counts_df[mask])
        elif mode.lower()=='qn':
            qn_df = rnaseqnorm.normalize_quantiles(self.tpm_df.loc[mask])
            self.norm_df = rnaseqnorm.inverse_normal_transform(qn_df)
        else:
            raise ValueError('Unsupported mode {}'.format(mode))

        assert(self.norm_df is not None)

        return(0)

    def __normalize(self, rid):

        self.__read_expression(rid)

        print('Normalizing data ({})'.format(args.normalization_method), flush=True)

        self.__apply_normalization(sample_frac_threshold=args.sample_frac_threshold, count_threshold=args.count_threshold, tpm_threshold=args.tpm_threshold, mode=args.normalization_method)
        print('  * {} genes in input tables.'.format(self.counts_df.shape[0]), flush=True)
        print('  * {} genes remain after thresholding.'.format(self.norm_df.shape[0]), flush=True)
    
        # change sample IDs to participant IDs
        self.norm_df.rename(columns=self.SP, inplace=True)
    
        print('Writing GE file.', flush=True)
        self.norm_df.to_csv(os.path.join(args.output_dir, rid+'.expression.txt.gz'), compression='gzip', sep='\t')
 
        if args.annotation_gtf:
            self.__gtf_to_bed(args.annotation_gtf)
            self.__prepare_bed()
            print('Writing BED file.', flush=True)
            self.__write_bed(os.path.join(args.output_dir, rid+'.expression.bed'))

    def __reset_containers(self):

        self.counts_df = None
        self.tpm_df = None
        self.norm_df = None
        self.bed_df = None

    def write_excluded_samples(self):

        print("\nWriting excluded samples as json.")

        with open(os.path.join(args.output_dir, 'excluded_samples.json'), 'w+') as f:
           json.dump(self.excluded_samples, f)
 
    def normalize(self):

        for rid in self.regions:
            self.__normalize(rid)
            self.__reset_containers()
   
parser = argparse.ArgumentParser(description='Generate normalized expression BED files for eQTL analyses.')
parser.add_argument('gct_dir', help='GCT file with expression in normalized units, e.g., TPM or FPKM.')
parser.add_argument('gtt_dir', help='GCT file with read counts')
parser.add_argument('sample_participant_lookup', help='Lookup table linking a desired subset of samples to participants.')
parser.add_argument('--sample_id_list', default=None, help='Subset of samples to normalize.')
parser.add_argument('--annotation_gtf', help='Normalized expression table also written to BED if GTF annotation provided.')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
parser.add_argument('--tpm_threshold', type=np.double, default=0.1, help='Selects genes with > expression_threshold expression in at least sample_frac_threshold')
parser.add_argument('--count_threshold', type=np.int32, default=6, help='Selects genes with >= count_threshold reads in at least sample_frac_threshold samples')
parser.add_argument('--sample_frac_threshold', type=np.double, default=0.2, help='Minimum fraction of samples that must satisfy thresholds')
parser.add_argument('--normalization_method', default='tmm', help='Normalization method: TMM or quantile normalization (qn)')
args = parser.parse_args()

def main():
    ne = NormalizeExpression(args.gct_dir, args.gtt_dir, args.sample_participant_lookup, args.sample_id_list)
    ne.normalize()
    ne.write_excluded_samples()

    print("\n[ {} ] Finished normalization.".format(datetime.now().strftime('%b %d %H:%M:%S')))
 

if __name__=="__main__":
    main()



#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import numpy as np
from collections import defaultdict
import logging
import tempfile
import shutil
import queue
import json
import time
import os
import sys

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)

log = logging.getLogger(__name__)

def parse_args():

    parser = argparse.ArgumentParser(description='Prepare METASOFT input')
    parser.add_argument('gene_chunk_dir', help='Directory with study-specific gene chunk directories.')
    parser.add_argument('annot_gene_list', type=str, help='ENSEMBL to gene symbol map.')
    parser.add_argument('-c', '--chunksize', type=int, default=25, help='Number of genes to merge per blade on MGI.')
    parser.add_argument('--log-level', default='INFO', choices=['NOTSET','DEBUG','INFO', 'WARNING','CRITICAL',
                                                        'ERROR','CRITICAL'], help='Log level')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory.')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    args.output_dir = os.path.abspath(args.output_dir)

    log.setLevel('INFO')
    log.info(sys.argv)    

    return(args)

def bsub(cmd, job_name, mem=60, gtmp=4, docker_image="apollodorus/bioinf:pr", queue="research-hpc"):
    """ Creates a bsub command for LSF job submission.
    Args:
        cmd: command to be run.
        queue: queue to submit job
        job_name: name of job
        mem: RAM space to request from scheduler (MB)
        gtmp: temporary space to request from scheduler (GB)
        docker_image: [username]/[repos]:[tag]
    Returns:
        string: bsub command string.
    """

    mem = str(mem)
    gtmp = str(gtmp)
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -Is -q "+queue \
        + " -J "+job_name \
        + " -M "+mem+"000000" \
        + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
        + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
        + " -a 'docker("+docker_image+")'" \
        + " /bin/bash -c '"+cmd+"'" 

    
    return(bsub)

def loop_chunks(l, n, d):

    # l: list, n: chunksize, d: defaultdict(lambda: dict())

    for num, i in enumerate(range(0, len(l), n)): 

        chunk_files = l[i:i+n]

        chunk_name = 'chunk_'+str(num+1)

        d[chunk_name] = {'return_code':None, 'genes':chunk_files, 'log':None}

        yield(chunk_name, chunk_files)

def merge_gene_chunk(gene_list, gene_json, chunk_name, gene_outdir, process_container, process_dict, log_fp):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'merge_gene_wise.py')
    cmd = " ".join([exc, gene_list, gene_json, '-o', gene_outdir])

    with open(log_fp, 'w+') as fp:
        po = subprocess.Popen(bsub(cmd, 'build_metasoft'), stdout=fp, stderr=subprocess.STDOUT, shell=True)
        process_container.put((chunk_name, po))

def main(args):
         
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    def region_name(fp):
        return(os.path.split(fp)[1])

    region_list = ['DLPFC', 'PAR', 'TCTX', 'PFCTX', 'BM22']
## DEBUG
#    region_list = ['TYPE1', 'TYPE2', 'TYPE3'] 

    # 1. Get a tuple list of region names and their corresponding directories 

    region_gene_dir_list= [(region_name(i), os.path.join(os.path.abspath(args.gene_chunk_dir),i)) for i in os.listdir(args.gene_chunk_dir) if region_name(i) in region_list]

    assert all([i in region_list for i,j in region_gene_dir_list])

    # 2. Get intersection of genes across regions

    genes = set()

    ## a. First find union of genes 

    for _ ,region_gene_dir in region_gene_dir_list:
        genes.update([os.path.splitext(gene)[0] for gene in os.listdir(region_gene_dir) if gene.startswith('ENSG')])

    ## b. Next obtain intersection across all regions from the union

    for _ ,region_gene_dir in region_gene_dir_list:
        genes.intersection_update([os.path.splitext(gene)[0] for gene in os.listdir(region_gene_dir) if gene.startswith('ENSG')])
    
    genes = list(genes)
    log.info('{} genes form intersection across regions {}'.format(str(len(genes)), ' '.join(region_list)))
    
    with open(args.annot_gene_list) as fp:
        annot_genes = [row.split('\t')[0] for row in fp.read().strip().split('\n')]

    unfiltered_num = len(genes)
    genes = [i for i in genes if i in annot_genes]

    log.info('{} genes filtered.'.format(str(unfiltered_num - len(genes))))
    log.info('{} genes remaining.'.format(len(genes)))

    ## c. Get dictionary of paths to gene-level association data per region 

    region_gene_paths = defaultdict(lambda: dict())
    for gene in genes:
        region_gene_paths[gene] = {region_id: os.path.join(region_gene_dir, gene + '.gz') for region_id, region_gene_dir in region_gene_dir_list}
        assert all([os.path.exists(fp) for fp in region_gene_paths[gene].values()])

    gene_json = os.path.join(args.output_dir, 'region_gene_paths.json')
    with open(gene_json, 'w+') as fp:
        json.dump(dict(region_gene_paths), fp)

    # 3. Run parallelization of gene-wise merging across all regions

    tmpdir = tempfile.mkdtemp(dir=args.output_dir)
    
    log_file_dir = os.path.join(args.output_dir, 'logs') 
    merge_outdir = os.path.join(args.output_dir, 'region_merged_gene_chunks')

    chunk_dict = defaultdict(lambda: dict())
    chunk_processes = queue.Queue()

    if not os.path.exists(log_file_dir):
        os.makedirs(log_file_dir)

    for k,(chunk_name, gene_list) in enumerate(loop_chunks(genes, args.chunksize, chunk_dict)):

        log_outfile = os.path.join(log_file_dir, chunk_name +'.log')
        chunk_dict[chunk_name]['log'] = log_outfile

        gl_path = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)
        with open(gl_path.name, 'w') as fp:
            fp.write('\n'.join(gene_list))

        merge_gene_chunk(gl_path.name, gene_json, chunk_name, merge_outdir, chunk_processes, chunk_dict, log_outfile)

if __name__ == '__main__':

    args = parse_args()
    main(args)

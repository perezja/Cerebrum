import subprocess
import logging
import sys
import re
import os

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)

log = logging.getLogger(__name__)

def bsub(cmd, job_name, log_fp, docker_image="apollodorus/eqtl-coloc:v1", mem=8, gtmp=2, queue="research-hpc"):
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
    bsub = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q "+queue \
        + " -J "+job_name \
        + " -M "+mem+"000000" \
        + " -o "+log_fp \
        + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
        + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
        + " -a 'docker("+docker_image+")'" \
        + " /bin/bash -c '"+cmd+"'"

    return(bsub)

def harmonise_snps(exposure_files, outcome_dat, outcome_id, prefix, output_dir):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'harmonise_snps.R')
    cmd = ' '.join([exc, exposure_files, outcome_dat, outcome_id, prefix, '-o', output_dir])

    subprocess.check_call(cmd, shell=True)

    return(os.path.join(output_dir, prefix+'_harmonise.txt'))

def run_ld_clumping(bfile, assoc_file, prefix, output_dir, p1='5e-8', p2='5e-8', pval_field='pval.exposure', snp_field='SNP_ID', kb='500', r2='0.1'):

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

def run_twosamplemr(arg_list):

    exc = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_ivw_mr.R')
    cmd = ' '.join([exc] + arg_list)

    return(subprocess.check_call(cmd, stderr=subprocess.STDOUT, shell=True))

def loop_chunks(df, n, d):

    # l: list, n: chunksize, d: defaultdict(lambda: dict())
    
    df = df.sort_values(by='gene')

    l = df['gene'].unique().tolist()

    for num, i in enumerate(range(0, len(l), n)): 

        genes = l[i:i+n]

        chunk_df = df[df['gene'].isin(genes)] 

        chunk_name = 'chunk'+str(num+1)

        d[chunk_name] = {'return_code':None, 'genes':genes, 'log':None}

        yield(chunk_name, chunk_df)

def merge_chunks(chunk_dir_loc, header, outfile):
    """ process to concatenate all finished chunk output files."""

    with open(outfile, mode='wt') as fp:
        fp.write('\t'.join(header)+'\n')

    # find: full path to all files in specified directory
    # tail: cat all but header line
    cmd = "find {} -maxdepth 1 -type f | xargs -I{{}} tail -n+2 {{}} | cat >> {}".format(chunk_dir_loc, outfile) 

    subprocess.check_call(bsub(cmd, 'merging'), shell=True)

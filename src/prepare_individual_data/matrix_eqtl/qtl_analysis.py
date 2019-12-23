import argparse
import os, re
from datetime import datetime
import concurrent.futures
import json

from parallelize import Parallelize

class RunMatrixQTL(Parallelize):
    """Run MatrixQTL for all brain regions.

    Attributes: 
        region_id: defines which region the instance is processing
        cmd_array: a tuple array holding region id (rid) and corresponding 
        suffix: outfile suffix

            MatrixQTL command.  tuple: (rid,  cmd)
        header: string used when concatenating chunk output if LSF used.
    """
    
    def __init__(self, arg_tuple):
        Parallelize.__init__(self, args.chunksize, args.output_dir, args.tmp_dir)
         
        self.cmd_array = arg_tuple 
        self.suffix = ".assoc.txt.gz"

        self.header = ['SNP','gene','beta','t-stat','p-value','FDR']

    @property
    def cmd_array(self):
        return(self.__cmd_array)
    @cmd_array.setter
    def cmd_array(self, arg_tuple):
        """Sets cmd_array and region_id attributes."""

        cmd_array = []

        exc = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                'run_matrixQTL.R')

        self.region_id = arg_tuple[0] 

        expression_file = arg_tuple[1]
        covariates_file = arg_tuple[2]
        snps_file = arg_tuple[3]

        print('creating chunks for {}'.format(expression_file))
        for chunk_name, chunk_file in self.create_chunks(expression_file):   # base class method

            job_name = self.region_id + '_' + chunk_name

            cmd = ' '.join([exc, chunk_file, covariates_file, snps_file, job_name]) 
            cmd += ' --output_dir ' + self.chunk_outdir   # base class attribute

            if args.lsf:
                cmd_array.append((job_name, self.bsub(cmd, job_name, mem=75, gtmp=12)))
            else:
                cmd_array.append((job_name, cmd))

        self.__cmd_array = cmd_array

    def run(self):

        print("[ {} ] Started processing region '{}'".format(datetime.now().strftime('%b %d %H:%M:%S'),self.region_id))

        super().run(self.cmd_array)
 #       super().merge_chunks(self.header, self.region_id, self.suffix)
 #       super().merge_logs(self.region_id)

        print("[ {} ] Finished region '{}'".format(datetime.now().strftime('%b %d %H:%M:%S'),self.region_id))

    def __del__(self):

        super().write_job_monitor(self.region_id)

def arg_array(paths):
    """ Generator for fetching command arguments for 'run_matrixQTL.R' 

    Args: 
        paths: dictionary with input directory paths

    Yields: 
        tuple: (region id, expression file, covariance file, snps file)
    """

    def id(i):
        return(i.split('.')[0])  
    def fp(path, fn):
        return(os.path.join(path, fn))

    exp_dir = paths.get('exp_dir')
    exp_files = [(id(i), fp(exp_dir,i)) for i in os.listdir(exp_dir) if i.endswith('.txt.gz')]
    exp_files.sort(key=lambda x: x[0])

    cov_dir = paths.get('cov_dir')
    cov_files = [(id(i), fp(cov_dir,i)) for i in os.listdir(cov_dir) if i.endswith('.txt')]
    cov_files.sort(key=lambda x: x[0])

    snp_dir = paths.get('snp_dir')
    snp_files = [(id(i), fp(snp_dir,i)) for i in os.listdir(snp_dir) if i.endswith('.txt')]
    snp_files.sort(key=lambda x: x[0])

    assert len(exp_files) == len(cov_files) == len(snp_files)

    for exp_tuple, cov_tuple, snp_tuple in zip(exp_files, cov_files, snp_files): 
        yield (exp_tuple[0], exp_tuple[1], cov_tuple[1], snp_tuple[1]) 
    
def process_region(arg_tuple):
    """Put class instantiation in function call for multiprocessing. """

    mqtl = RunMatrixQTL(arg_tuple)
    mqtl.run()

parser = argparse.ArgumentParser(description='Run Matrix-EQTL in Parallel.')
parser.add_argument('json', help='JSON file containing paths to expression, covariate and genotype directories.')
parser.add_argument('chunksize', type=int, help='Number of chunks for splitting expression tables.')
parser.add_argument('--regions', nargs='+', help='Subset of regions to run.')
parser.add_argument('--lsf', action='store_true', help='Flag specifying whether to run in parallel on LSF.')
parser.add_argument('--max_workers', type=int, default=600, help='Flag specifying whether to run in parallel on LSF.')
parser.add_argument('--tmp_dir', help='Directory for temporary files.')
parser.add_argument('-o', '--output_dir', default='./qtls', help='Directory for output files.')
args = parser.parse_args()

if not args.tmp_dir:
    args.tmp_dir = args.output_dir

def main():

    with open(args.json) as fp:
        paths = json.load(fp)

    executor = concurrent.futures.ProcessPoolExecutor(max_workers=args.max_workers)
    if args.regions:
        futures = [executor.submit(process_region, arg_tuple) for arg_tuple in arg_array(paths) if arg_tuple[0] in args.regions] 
    else:
        futures = [executor.submit(process_region, arg_tuple) for arg_tuple in arg_array(paths)] 

    concurrent.futures.wait(futures)

if __name__ == "__main__":
    main()

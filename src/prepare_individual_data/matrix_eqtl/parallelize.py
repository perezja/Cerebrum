import pandas as pd
import subprocess
import tempfile
import shutil
import queue
import json
import gzip
import time
import os

class Parallelize():
    """A base class for parallelizing jobs in LSF.

    Attributes:

        chunksize: size of chunks to process from input dataframe

        outdir: Output directory for final merged output 
        tmpdir: Temporary directory for chunk output
        chunk_indir: Directory for input dataframe chunks 
        chunk_outdir: Directory for processed output chunks
        chunk_logdir: Directory for chunk logs

        job_monitor: dictionary holding return code status and log file
            for each chunk
    """
    def __init__(self, chunksize, outdir, tmpdir):

        self.chunksize = chunksize 
        self.tmpdir = tmpdir 
        self.outdir = outdir

        self.job_monitor = dict()

    @property
    def tmpdir(self):
        return(self.__tmpdir)
    @tmpdir.setter
    def tmpdir(self, tmpdir):
        """ Sets chunk directory attributes.""" 

        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)

        self.__tmpdir = tempfile.mkdtemp(dir=tmpdir) 

        self.chunk_indir = os.path.join(self.__tmpdir, 'chunk_indir')
        self.chunk_outdir = os.path.join(self.__tmpdir, 'chunk_outdir')
        self.chunk_logdir = os.path.join(self.__tmpdir, 'chunk_logs')

        chunk_dirs = [self.chunk_indir, self.chunk_outdir, self.chunk_logdir]

        for path in chunk_dirs:
            os.mkdir(path)

    @property
    def outdir(self):
        return(self.__outdir)
    @outdir.setter
    def outdir(self, x):
        if not os.path.exists(x):
            os.makedirs(x)
        self.__outdir = x  

    def loop_chunks(self, l, n):
        """ A generator for yielding array chunks. 
        Args:
            l: array-like
            n: chunk size

        Yields: 
            list: subset of 'l'     
        """
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def create_chunks(self, df_path):
        """ A generator for yielding df chunks.
        
        Args: 
            df_path: file path to pandas dataframe

        Returns: 
            tuple: (chunk_name, chunk_file_path).
        """

        if os.path.splitext(df_path)[1] == '.gz':
            self.df = pd.read_csv(df_path, sep='\t', compression='gzip', index_col=0) 
        else:
            self.df = pd.read_csv(df, sep='\t', index_col=0) 

        rows = self.df.shape[0]
        n = self.chunksize
        for num, i in enumerate(range(0, rows, n)):

            chunk_name = 'chunk'+str(num+1)
            chunk_path = os.path.join(self.chunk_indir, chunk_name+'.txt')

            df = self.df.iloc[i:i+n]
            df.to_csv(chunk_path, sep='\t')

            yield (chunk_name, chunk_path) 

    def bsub(self, cmd, job_name, queue="research-hpc", mem=65, gtmp=20, docker_image="apollodorus/brain-eqtl:eqtl-analysis"):
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
        bsub = "bsub -q "+queue \
            + " -J "+job_name \
            + " -M "+mem+"000000" \
            + " -R 'select[mem>"+mem+"000 && gtmp > "+gtmp+"]" \
            + " rusage[mem="+mem+"000, gtmp="+gtmp+"]'" \
            + " -a 'docker("+docker_image+")'" \
            + " /bin/bash -c '"+cmd+"'" 

        return(bsub)

    def run(self, cmd_array):
        """run an array of commands and wait until jobs completed.
        Args: 
            cmd_array: tuple list (job_name, cmd).
        Returns:
            None

        Updates: job_monitor dictionary as chunk runs finish.

        Raises: Exception if chunk fails with non-zero return code. 
        """

        child_processes = queue.Queue()

        for job_name, cmd in cmd_array:

            log_path = os.path.join(self.chunk_logdir, job_name+".log")
            log = open(log_path, "w+")

            self.job_monitor[job_name] = {"return_code": None, "cmd":cmd, "log":log_path} 

            po = subprocess.Popen(cmd, stdout=log, stderr=log, shell=True) 
            child_processes.put((job_name, po))

        while not child_processes.empty():

            job_name, po = child_processes.get() 
            
            while po.poll() is None:
                time.sleep(0.5)

            self.job_monitor[job_name]["return_code"] = po.returncode 

            if po.returncode == 0:

                n_jobs = len(self.job_monitor.keys())
                n_complete =  n_jobs - child_processes.qsize()

            else:

                returncode = self.job_monitor[job_name]["return_code"]
                with open(self.job_monitor[job_name]["log"], "r") as f:
                    print(f.read())
 
                raise Exception(" {} failed with exit code {}.".format(job_name, returncode))

    def write_job_monitor(self, prefix):
        """ write out job monitoring dictionary for every chunk. """

        with open(os.path.join(self.outdir, prefix+'_chunks.json'), 'w+') as fp:
            json.dump(self.job_monitor, fp)

    def merge_chunks(self, header, prefix, suffix):
        """ process to concatenate all finished chunk output files."""
        
        fn = os.path.join(self.outdir, prefix + suffix)

        with gzip.open(fn, mode='wt') as gz:
            gz.write('\t'.join(header)+'\n')

        # find: full path to all files in specified directory
        # tail: cat all but header line
        cmd = "find {} -maxdepth 1 -type f | sort -V | xargs -I{{}} tail -n +2 {{}} | gzip -c >> {}".format(self.chunk_outdir, fn) 

        subprocess.check_call(cmd, shell=True)

    def merge_logs(self, prefix):
        """ process to concatenate all chunk log files."""
        
        fn = os.path.join(self.outdir, prefix+'.logs.txt.gz')

        cmd = "find {} -maxdepth 1 -type f | sort -V | xargs -I{{}} cat {{}} | gzip -c >> {}".format(self.chunk_logdir, fn)  
        subprocess.check_call(cmd, shell=True)

    def __del__(self):

        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)



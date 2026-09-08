#coding=utf-8
from __future__ import print_function
import argparse
import logging
import os
from collections import defaultdict
example_text = '''example:
    python this.py 
'''
parser = argparse.ArgumentParser(description="The script is to .",
                                formatter_class= argparse.RawTextHelpFormatter,
                                usage = '%(prog)s [-h]',                         
                                epilog=example_text)

parser.add_argument('--sample_infomation','-i',type=str,help="the sample information file (sep by tab), colnames: sample_name    fq1 fq2",required= True,metavar='')
parser.add_argument('--config_yaml', '-c',type= str,help="the config.yaml file template",required= True,metavar='')
parser.add_argument('--run_snakemake_sh', '-r',type= str,help="the run_snakemake.sh file template",required= True,metavar='')
parser.add_argument('--path', '-p',type= str,help="the main output directory of the pipeline",required= True,metavar='')
parser.add_argument('--set_cores','-j',type=str,help="set the number of cores (threads) for snakemake",required= True,metavar='')
args = parser.parse_args()

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def create_dirs(dirlist):
    if not isinstance(dirlist, list):
        dirlist = [dirlist]
    for dirname in dirlist:
        if not os.path.isdir(dirname):
            logger.info("Creating directory %s" % (dirname))
            os.makedirs(dirname)
        else:
            logger.info("Directory %s exists" % (dirname))

def delete_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)
        logger.info(f"File '{file_path}' has been deleted and will be recreated next time.")
    else:
        logger.info(f"File '{file_path}' does not exist and will be created.")
            
def sample_info_dict(sample_infomation_file):
    """
    convert the sample information file to a dictionary
    """
    samp_info_dict = defaultdict(list)
    with open(sample_infomation_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                lines = line.strip().split('\t')
                samp_info_dict[lines[0].strip()].extend([lines[1].strip(), lines[2].strip()] )
        return samp_info_dict

def change_config_file(path, samplename, fq1, fq2, run_snakemake_sh_file, config_yaml_file, cores):
    config_yaml_file_path = os.path.join(path, 'config.yaml')
    _run_name = samplename + '_' + 'run_snakemake.sh'
    run_snakemake_sh_file_path = os.path.join(path, _run_name)
    with open(run_snakemake_sh_file_path, 'w') as outf:
        with open(run_snakemake_sh_file, 'r') as inpf:
            for line in inpf:
                if line.strip(" ").startswith("pwd"):
                    print('pwd='+str(path),file = outf)
                elif line.strip(" ").startswith("threads"):
                    print('threads='+str(cores),file = outf)
                else:
                    outf.write(line)
                    
    with open(config_yaml_file_path, 'w') as outf:
        with open(config_yaml_file, 'r') as configf:
            for line in configf:
                if line.strip(" ").startswith("id"):
                    print('    id: "'+samplename+'"',file = outf)
                elif line.strip(" ").startswith("raw_fq1"):
                    print('    raw_fq1: '+str(fq1),file = outf)
                elif line.strip(" ").startswith("raw_fq2"):
                    print('    raw_fq2: '+str(fq2),file = outf)
                elif line.strip(" ").startswith("threads"):
                    print('threads: '+str(cores),file = outf)
                else:
                    outf.write(line)
    

def create_all_dirs_files(samp_info_dict, path, config_yaml_file, run_snakemake_sh_file, cores):
    '''
    splitting the script into two parts by depth of bam files, 
    inorder to run the pipeline faster by considering how many jobs took to run in parallel.
    '''
    output_dir = os.path.join(path, 'output/01.fq2methRatio')
    shell_dir = os.path.join(path, 'shell')
    
    create_dirs([output_dir, shell_dir])
    
    delete_file(os.path.join(shell_dir, 'all2run.sh'))
    # create all file in output
    for sample in samp_info_dict:
        sample_dir = os.path.join(output_dir, sample)
        create_dirs([sample_dir])
        fq1 = samp_info_dict[sample][0]
        fq2 = samp_info_dict[sample][1]

        change_config_file(sample_dir, sample, fq1, fq2, run_snakemake_sh_file, config_yaml_file, cores)
    # create all file in shell
        with open(os.path.join(shell_dir, 'all2run.sh'), 'a') as outf:
            _run_name = os.path.join(sample_dir, str(sample) + '_' + 'run_snakemake.sh')
            print(f'bash {_run_name}', file = outf)
            
if __name__ == '__main__':
    sample_infomation_file = args.sample_infomation
    config_yaml_file = args.config_yaml
    run_snakemake_sh_file = args.run_snakemake_sh
    cores = args.set_cores
    path = os.path.abspath(args.path)
    samp_info_dict = sample_info_dict(sample_infomation_file)
    create_all_dirs_files(samp_info_dict, path, config_yaml_file, run_snakemake_sh_file, cores)
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

parser.add_argument('--sample_infomation','-i',type=str,help="the sample information file separated by tab: sample_name    R1.fastq.gz    R2.fastq.gz   **",required= True,metavar='')
parser.add_argument('--config_yaml', '-c',type= str,help="the config.yaml file template",required= True,metavar='')
parser.add_argument('--run_snakemake_sh', '-r',type= str,help="the run_snakemake.sh file template",required= True,metavar='')
parser.add_argument('--output_dir', '-p',type= str,help="the path of the output directory",required= True,metavar='')
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
            
# def sample_info_dict(sample_infomation_file):
#     """
#     convert the sample information file to a dictionary
#     """
#     samp_info_list = []
#     with open(sample_infomation_file, 'r') as f:
#         for line in f:
#             if not line.startswith('#'):
#                 sampleName, R1, R2, *_ = line.strip().split('\t')
#                 samp_info_list.append([sampleName, R1, R2])
#         return samp_info_list

def change_config_file(path, sampleInfoList, run_snakemake_sh_file, config_yaml_file):
    config_yaml_file_path = os.path.join(path, 'config.yaml')
    _run_name = sampleInfoList[0] + '_' + 'run_snakemake.sh'
    run_snakemake_sh_file_path = os.path.join(path, _run_name)
    with open(run_snakemake_sh_file_path, 'w') as outf, open(run_snakemake_sh_file, 'r') as inpf:
        for line in inpf:
            if line.strip(" ").startswith("pwd"):
                print('pwd='+str(path),file = outf)
            else:
                outf.write(line)
    with open(config_yaml_file_path, 'w') as outf, open(config_yaml_file, 'r') as configf:
        for line in configf:
            if line.strip(" ").startswith("id"):
                print('    id: "'+sampleInfoList[0]+'"',file = outf)
            elif line.strip(" ").startswith("clean_fq1_path"):
                print('    clean_fq1_path: '+str(sampleInfoList[1]),file = outf)
            elif line.strip(" ").startswith("clean_fq2_path"):
                print('    clean_fq2_path: '+str(sampleInfoList[2]),file = outf)
            else:
                outf.write(line)
    
def create_all_dirs_files(sample_infomation_file, pwd, config_yaml_file, run_snakemake_sh_file):
    '''
    splitting the script into two parts by depth of bam files, 
    inorder to run the pipeline faster by considering how many jobs took to run in parallel.
    '''
    output_dir = os.path.join(pwd, 'output/')
    shell_dir = os.path.join(pwd, 'bin/output_sh')
    log_dir = os.path.join(pwd, 'bin/output_sh/log')
    create_dirs([log_dir])
    
    with open(sample_infomation_file, 'r') as f, open(os.path.join(shell_dir, 'all_run.sh'), 'a') as outf:
        for line in f:
            if not line.startswith('#'):
                sample, R1, R2, *_ = line.strip().split('\t')
            sample_dir = os.path.join(output_dir, sample)
            create_dirs([sample_dir])
            change_config_file(sample_dir, sample, run_snakemake_sh_file, config_yaml_file)
            # create all file in shell
            _run_name = os.path.join(sample_dir, str(sample) + '_' + 'run_snakemake.sh')
            print(f'(/APP/u22/x86/bin/time -f "MaxMemory_KB:%M\\nRealTime_sec:%e\\nUserTime_sec:%U" sh {_run_name} && echo "{sample} Done!" || echo "{sample} Undone") > {log_dir}/{sample}.log 2>&1', file = outf)
            
if __name__ == '__main__':
    delete_file(os.path.join(args.PWD, 'bin/output_sh/all_run.sh'))
    pwd = os.path.abspath(args.PWD)
    create_all_dirs_files(args.sample_infomation, pwd, args.config_yaml, args.run_snakemake_sh)

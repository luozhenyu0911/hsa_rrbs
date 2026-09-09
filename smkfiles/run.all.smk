# include the config file
configfile: "config.yaml"

# define a function to return target files based on config settings
def run_all_input(wildcards):

    run_all_files = []
    # run_all_files.append("01.qc/{}.multiqc.html".format(config['samples']['id']))
    run_all_files.append("04.metrics/{}_lambda_PE_report.txt".format(config['samples']['id']))
    run_all_files.append("04.metrics/{}_conversion_rate.txt".format(config['samples']['id']))
    run_all_files.append("{}_sum_info.txt".format(config['samples']['id']))


    # run_all_files.append("03.methylation/{}_pe.bismark.cov.gz".format(config['samples']['id']))
    # run_all_files.append("03.methylation/{}_pe_splitting_report.txt".format(config['samples']['id']))
    # run_all_files.append("03.methylation/{}.CG_5x.bedgraph.gz".format(config['samples']['id']))
    # run_all_files.append("04.metrics/{}.CG_depth.stat.txt".format(config['samples']['id']))
    # run_all_files.append("02.bismark_bt2/{}_pe.cram".format(config['samples']['id']))


    return run_all_files


# rule run all, the files above are the targets for snakemake
rule run_all:
    input:
        run_all_input
		
smk_path = config['params']['smk_path']
include: smk_path+"/map.smk"
include: smk_path+"/metrics.smk"

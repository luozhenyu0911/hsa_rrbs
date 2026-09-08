
from pathlib import Path

def modify_config(config_file, output_file, sample_name, column):
    """
    """
    context = Path(output_file).parts[-2]
    if context == 'CG':
        inputf = "/XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/input/Control_vs_Case_CG_input.txt"
    elif context == 'CHG':
        inputf = "/XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/input/Control_vs_Case_CHG_input.txt"
    elif context == 'CHH':
        inputf = "/XYFS02/HDD_POOL/gzfezx_shhli/gzfezx_shhlixy_1/Projects/20260514_S256rrbs/input/Control_vs_Case_CHH_input.txt"

    with open(config_file, 'r') as f, open(output_file, 'w') as outputf:
        for line in f:
            if line.startswith('column'):
                print(f'column={str(column)}', file=outputf)
            elif line.startswith('name'):
                print(f'name={sample_name}', file=outputf)
            elif line.startswith('inputfile'):
                print(f'inputfile={inputf}', file=outputf)
            else:
                outputf.write(line)

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Modify config file for methylation pipeline')
    parser.add_argument('-c','--config_file',type=str, required=True, help='Path to config file')
    parser.add_argument('-o','--output_file',type=str, required=True, help='Path to output file')
    parser.add_argument('-n','--sample_name',type=str, required=True, help='Sample name')
    parser.add_argument('-col','--column',type=str, required=True, help='Column name')
    args = parser.parse_args()
    modify_config(args.config_file, args.output_file, args.sample_name, args.column)

if __name__ == '__main__':
    main()
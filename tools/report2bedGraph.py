import argparse
import gzip

def main():
    parser = argparse.ArgumentParser(description='Convert a report file to a bedGraph file.')
    parser.add_argument('--report_file', '-r', required=True, help='Path to the report file.')
    parser.add_argument('--bedGraph_file', '-b', required=True, help='Path to the bedGraph file to be created.')
    parser.add_argument('--min_count', '-c', type=int, default=5, help='Minimum count threshold for bedGraph entries. Default is 5.')
    args = parser.parse_args()

    # Read the report file
    with gzip.open(args.report_file, 'rt') as inf, gzip.open(args.bedGraph_file, 'wt') as outf:
        for line in inf:
            # <chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>
            chromosome, start, end, methylation_percentage, methylated_count, non_methylated_count = line.strip().split('\t')
            total = int(methylated_count) + int(non_methylated_count)
            methylation_rate = int(methylated_count) / total if total > 0 else 0
            
           # Write to bedGraph file if the methylation rate is above the threshold
            if (total >= args.min_count) and methylation_rate > 0 :
                outf.write(f'{chromosome}\t{int(start)-1}\t{int(end)}\t{methylation_rate}\n')

if __name__ == '__main__':
    main()
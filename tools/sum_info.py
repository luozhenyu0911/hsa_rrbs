import argparse

def main():
    
    parser = argparse.ArgumentParser(description="Summarize information from a CSV file.")
    parser.add_argument("--sampleNames", '-s', type=str, required=True, help="the sample names")
    parser.add_argument("--output",'-o', type=str, required=True, help="Path to the output CSV file.")
    args = parser.parse_args()
    
    with open(f'{args.output}', 'w') as outputf:
        print(f'Sample Name\t{args.sampleNames}', file=outputf)
        
        # trims
        with open(f'01.qc/{args.sampleNames}.R2.fq.gz_trimming_report.txt', 'r') as f:
            for line in f:
                if line.startswith('Total reads processed'):
                    total_reads = line.split(':')[1].strip()
                    print(f'Total reads\t{total_reads}', file=outputf)
                elif line.startswith('Reads written (passing filters)'):
                    reads_written = line.split(':')[1].strip()
                    print(f'Reads (passing filters)\t{reads_written}', file=outputf)
                    
        # alignment
        with open(f'02.bismark_bt2/{args.sampleNames}_PE_report.txt', 'r') as f:
            for line in f:
                if line.startswith('Sequence pairs analysed in total'):
                    total_pairs = line.split(':')[1].strip()
                    print(f'Total Sequence pairs\t{total_pairs}', file=outputf)
                elif line.startswith('Number of paired-end alignments with a unique best hit:'):
                    mapped_pairs = line.split(':')[1].strip()
                    print(f'Mapped Sequence pairs\t{mapped_pairs}', file=outputf)
                elif line.startswith('Mapping efficiency:'):
                    mapped_efficiency = line.split(':')[1].strip()
                    print(f'Mapping efficiency\t{mapped_efficiency}', file=outputf)
        
        # methylation
        with open(f'03.methylation/{args.sampleNames}_pe_splitting_report.txt', 'r') as f:
            for line in f:
                if line.startswith("Total methylated C's in CpG context:"):
                    methylated_C = line.split(':')[1].strip()
                    print(f'Total methylated C in CpG\t{methylated_C}', file=outputf)
                elif line.startswith('C methylated in CpG context:'):
                    methylated_C_ratio = line.split(':')[1].strip()
                    print(f'C methylated in CpG\t{methylated_C_ratio}', file=outputf)
                elif line.startswith("Total methylated C's in CHG context:"):
                    methylated_CHG = line.split(':')[1].strip()
                    print(f'Total methylated C in CHG\t{methylated_CHG}', file=outputf)
                elif line.startswith('C methylated in CHG context:'):
                    methylated_CHG_ratio = line.split(':')[1].strip()
                    print(f'C methylated in CHG\t{methylated_CHG_ratio}', file=outputf)
                elif line.startswith("Total methylated C's in CHH context:"):
                    methylated_CHH = line.split(':')[1].strip()
                    print(f'Total methylated C in CHH\t{methylated_CHH}', file=outputf)
                elif line.startswith('C methylated in CHH context:'):
                    methylated_CHH_ratio = line.split(':')[1].strip()
                    print(f'C methylated in CHH\t{methylated_CHH_ratio}', file=outputf)
        
        # conversion rate
        with open(f'04.metrics/{args.sampleNames}_conversion_rate.txt', 'r') as f:
            for line in f:
               if line.startswith('Conversion rate of total C base:'):
                    conversion_rate = line.split(':')[1].strip()
                    print(f'Conversion rate\t{conversion_rate}', file=outputf)
        
        # coverage
        with open(f'04.metrics/{args.sampleNames}.CG_depth.stat.txt', 'r') as f:
            for line in f:
                if line.startswith('all CpG Sites\tTotal_Sites'):
                    Total_Sites = line.split('\t')[2].strip()
                    print(f'Total CpG Sites\t{Total_Sites}', file=outputf)
                elif line.startswith('all CpG Sites\tAverage_Depth'):
                    Average_Depth = line.split('\t')[2].strip()
                    print(f'Average CpG Depth\t{Average_Depth}', file=outputf)
                elif line.startswith('all CpG Sites\t>=5x'):
                    gt5x = line.split('\t')[2].strip()
                    print(f'>=5x\t{gt5x}', file=outputf)
                elif line.startswith('CpG Sites(DP>5)\tTotal_Sites'):
                    Total_Sites = line.split('\t')[2].strip()
                    print(f'Total CpG Sites(>=5x)\t{Total_Sites}', file=outputf)
                elif line.startswith('CpG Sites(DP>5)\tAverage_Depth'):
                    Average_Depth = line.split('\t')[2].strip()
                    print(f'Average CpG Depth(>=5x)\t{Average_Depth}', file=outputf)

if __name__ == "__main__":
    main()

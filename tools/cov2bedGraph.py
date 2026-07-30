import argparse
import gzip

def main():
    parser = argparse.ArgumentParser(description='Convert a report file to a bedGraph file and calculate depth statistics.')
    parser.add_argument('--report_file', '-r', required=True, help='Path to the report file.')
    parser.add_argument('--bedGraph_file', '-b', required=True, help='Path to the bedGraph file to be created.')
    parser.add_argument('--stats_file', '-s', required=True, help='Path to the output statistics file.')
    parser.add_argument('--min_count', '-c', type=int, default=5, help='Minimum count threshold for bedGraph entries. Default is 5.')
    args = parser.parse_args()

    # 定义需要统计的深度阈值列表
    depth_thresholds = [5, 10, 15, 20, 25, 30]

    # 初始化统计字典
    # raw: 原始文件, filtered: 过滤后的文件
    stats = {
        'raw': {'total_sites': 0, 'total_depth': 0, 'threshold_counts': {t: 0 for t in depth_thresholds}},
        'filtered': {'total_sites': 0, 'total_depth': 0, 'threshold_counts': {t: 0 for t in depth_thresholds}}
    }

    # 读写文件
    with gzip.open(args.report_file, 'rt') as inf, gzip.open(args.bedGraph_file, 'wt') as outf:
        for line in inf:
            # 解析每一行
            chromosome, start, end, methylation_percentage, methylated_count, non_methylated_count = line.strip().split('\t')
            total = int(methylated_count) + int(non_methylated_count)
            methylation_rate = int(methylated_count) / total if total > 0 else 0
            
            # --- 1. 统计原始文件数据 ---
            stats['raw']['total_sites'] += 1
            stats['raw']['total_depth'] += total
            for t in depth_thresholds:
                if total >= t:
                    stats['raw']['threshold_counts'][t] += 1

            # --- 2. 过滤逻辑与统计过滤后数据 ---
            # 这里的过滤条件：深度 >= min_count 且 甲基化率 > 0
            if (total >= args.min_count) and methylation_rate > 0:
                # 写入 bedGraph (注意：保持您原本的 0-based start 转换)
                outf.write(f'{chromosome}\t{int(start)-1}\t{int(end)}\t{methylation_rate}\n')
                
                # 统计过滤后的数据
                stats['filtered']['total_sites'] += 1
                stats['filtered']['total_depth'] += total
                for t in depth_thresholds:
                    if total >= t:
                        stats['filtered']['threshold_counts'][t] += 1

    # --- 3. 计算并输出统计结果 ---
    with open(args.stats_file, 'w') as sf:
        sf.write("Category\tMetric\tValue\n")
        
        for group in ['raw', 'filtered']:
            group_name = "Raw_Data" if group == 'raw' else "Filtered_Data"
            total_sites = stats[group]['total_sites']
            total_depth = stats[group]['total_depth']
            
            # 计算平均深度
            avg_depth = total_depth / total_sites if total_sites > 0 else 0
            sf.write(f"{group_name}\tTotal_Sites\t{total_sites}\n")
            sf.write(f"{group_name}\tAverage_Depth\t{avg_depth:.2f}\n")
            
            # 计算各深度阈值的位点数和比例
            for t in depth_thresholds:
                count = stats[group]['threshold_counts'][t]
                percentage = (count / total_sites * 100) if total_sites > 0 else 0
                sf.write(f"{group_name}\t>={t}x\t{count}({percentage:.2f}%)\n")
            
            sf.write("-" * 40 + "\n") # 加上分割线方便阅读

    print(f"Report file processed successfully. Output files created: {args.stats_file}")

if __name__ == '__main__':
    main()
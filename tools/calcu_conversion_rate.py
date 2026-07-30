import sys

def safe_to_int(name, value):
    """尝试将值转为整数，成功则打印变量名:值，失败则报错"""
    if not value:
        raise ValueError(f"{name}: 变量为空，无法转换")
    
    try:
        result = int(value)
        print(f"{name}: {result}")
        return result
    except (ValueError, TypeError):
        raise ValueError(f"{name}: 无法将 '{value}' 转换为整数")
    
if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    if len(sys.argv) != 3:
        print("Usage: python calcu_conversion_rate.py input_file output_file")
        sys.exit(1)
    with open(input_file, "r") as inf, open(output_file, "w") as outf:
        for line in inf:
            line = line.strip()
            if line.startswith("Total methylated C's in CpG context:"):
                methylated_CpG = line.split(":")[1].strip()
            elif line.startswith("Total methylated C's in CHG context:"):
                methylated_CHG = line.split(":")[1].strip()
            elif line.startswith("Total methylated C's in CHH context:"):
                methylated_CHH = line.split(":")[1].strip()
            # unmethylated
            elif line.startswith("Total unmethylated C's in CpG context"):
                unmethylated_CpG = line.split(":")[1].strip()
            elif line.startswith("Total unmethylated C's in CHG context:"):
                unmethylated_CHG = line.split(":")[1].strip()
            elif line.startswith("Total unmethylated C's in CHH context:"):
                unmethylated_CHH = line.split(":")[1].strip()
        for name, var in [("methylated_CpG", methylated_CpG), ("methylated_CHG", methylated_CHG), ("methylated_CHH", methylated_CHH), ("unmethylated_CpG", unmethylated_CpG), ("unmethylated_CHG", unmethylated_CHG), ("unmethylated_CHH", unmethylated_CHH)]:
            safe_to_int(name, var)
        # conversion rate
        conversion_rate = 1 - (int(methylated_CpG) + int(methylated_CHG) + int(methylated_CHH))/(int(methylated_CpG)+ int(methylated_CHG) + int(methylated_CHH) + int(unmethylated_CpG) + int(unmethylated_CHG) + int(unmethylated_CHH))
        outf.write(f"Conversion rate of total C base: {conversion_rate}\n")
        # conversion rate
        



import sys

def load_mappings(mapping_file):
    mappings = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            original, replacement = line.strip().split('\t')
            mappings[original] = replacement
    return mappings

def replace_names_in_fasta(fasta_file, mappings, output_file):
    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                seq_name = line[1:].strip().split()[0]  # 仅获取序列名称的第一部分:line[1:]：去掉行首的 > 字符
                if seq_name in mappings:
                    f_out.write(f'>{mappings[seq_name]}\n')
                else:
                    f_out.write(line)
            else:
                f_out.write(line)

def main(mapping_file, fasta_file, output_file):
    mappings = load_mappings(mapping_file)
    replace_names_in_fasta(fasta_file, mappings, output_file)

if __name__ == "__main__":  #Python中的一种惯例，表示如果这个脚本是直接运行的，那么__name__变量将被设置为"__main__"
    if len(sys.argv) != 4:
        print("Usage: python replace_fasta_names.py <mapping_file> <fasta_file> <output_file>")
    else:
        mapping_file = sys.argv[1]
        fasta_file = sys.argv[2]
        output_file = sys.argv[3]
        main(mapping_file, fasta_file, output_file)

#!/usr/bin/env python
import pandas as pd
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import re
import argparse
import sys
from tqdm import tqdm

parser = argparse.ArgumentParser(
    description=f"RP8的IVD项目：自动化blast和操作比对结果\
    For example: python {sys.argv[0]} -i input.file -p sequence.path -o output -m max_target_seqs")
parser.add_argument("-i",type=str,
    help="由实验同时填写的输入文档，必须包含：病原种拉丁名和生产编号两列")
parser.add_argument("-p",type=str,
    help="一代测序的存放位置")   
parser.add_argument("-m",type=int,default=100,
    help="指定blast的max_target_seqs参数,默认是100"),
parser.add_argument("-o",type=str,
    help="输出结果文件名称")
args = parser.parse_args()



def get_path_and_seq(df, path):
    if '生产编号' and '病原种拉丁名' not in df.columns:
        raise ValueError("DataFrame 必须包含 '病原种拉丁名' 和 '生产编号' 列")  
    
    df['测序结果Seq文件名'] = ''
    df['产物序列'] = ''
    df['实际产物大小'] = 0

    for idx, row in df.iterrows():       # 遍历 DataFrame 中的每个生产编号
        prod_id = str(row['生产编号'])
        
        # 查找包含生产编号的文件
        found_file = None
        for filename in os.listdir(path):
            if re.search(prod_id, filename) and filename.endswith('.seq'):
                found_file = filename
                break
        
        if found_file:
            file_path = os.path.join(path, found_file)
            
            # 读取文件内容并更新 DataFrame
            with open(file_path, 'r') as file:
                sequence = file.read().strip()  # 读取序列数据并去除首尾空白 
                df.at[idx, '测序结果Seq文件名'] = found_file
                df.at[idx, '产物序列'] = sequence
                df.at[idx, '实际产物大小'] = len(sequence)
        else:
            df.at[idx, '测序结果Seq文件名'] = '未找到文件'
            df.at[idx, '产物序列'] = ''
            df.at[idx, '实际产物大小'] = 0
    return df



def run_local_blast(sequence, db_path, output_file):
    # 将序列写入临时查询文件
    query_file = 'temp_query.fasta'
    with open(query_file, 'w') as f:
        f.write(f'>query\n{sequence}\n')
   
    # 创建 BLAST 命令行
    blast_cline = NcbiblastnCommandline(
        query=query_file, 
        db=db_path, 
        evalue=1e-5, 
        outfmt="6 sacc sscinames qlen qcovs pident evalue",  # 指定输出格式及字段
        out=output_file, 
        max_target_seqs=args.m,  #j
        num_threads=20
    )
 
    stdout, stderr = blast_cline()
    os.remove(query_file)

    if stderr:
        print("BLAST stderr:", stderr)
    

def parse_blast_results(output_file):
    # 解析 BLAST 输出
    df = pd.read_csv(output_file, sep='\t', header=None, names=[
        'ACC号', '比对物种标题', '比对序列长度', '对比覆盖度', '对比相似度','Evalue'])
    df['pid_cov'] = (df['对比覆盖度'] * df['对比相似度'] / 100).round(2)
    
    group_stats = df.groupby('比对物种标题').agg(
        count=('ACC号', 'size'),
        max_pid_cov=('pid_cov', 'max')
    ).reset_index()
    
    total_count = group_stats['count'].sum()
    group_stats['percentage'] = (group_stats['count'] / total_count * 100).round().astype(int).astype(str) + '%'

    def format_row(row):
        return f"{row['比对物种标题']}|count:{row['count']}({row['percentage']})|pid_cov:{row['max_pid_cov']}"
    formatted_rows = group_stats.apply(format_row, axis=1).tolist()
    # 填充到最多3项
    while len(formatted_rows) < 3:
        formatted_rows.append('NA')
    formatted_rows = formatted_rows[:3]
    result_string = ';'.join(formatted_rows)
    
    return result_string



def parse_speci_blast_results(output_file,specie):
    df = pd.read_csv(output_file, sep='\t', header=None, names=[
    'ACC号', '比对物种标题', '比对序列长度', '对比覆盖度', '对比相似度','Evalue'])

    filtered_df = df[df['比对物种标题'].str.contains(specie, case=False, na=False)]
    if filtered_df.empty:
        cov = df.loc[df['Evalue'].idxmin(), '对比覆盖度']  # .idxmax() 返回具有最大值的索引位置
        ident = df.loc[df['Evalue'].idxmin(),"对比相似度"]
        result = f'Cov:{cov}|Ident:{ident}|非特异'
    else:
        cov = filtered_df.loc[filtered_df['Evalue'].idxmin(), '对比覆盖度']  # .idxmax() 返回具有最大值的索引位置
        ident = filtered_df.loc[filtered_df['Evalue'].idxmin(),"对比相似度"]
        result = f'Cov:{cov}|Ident:{ident}|特异'
    return result




def blast_sequences(df, db_path, output_folder):
    # 确保输出文件夹存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # 处理每一条序列
    # for idx, row in df.iterrows(): df.shape[0] 返回df有多少行
    for idx, row in tqdm(df.iterrows(), total=df.shape[0], desc="Processing Sequences"):
        sequence = row['产物序列']
        query_id = row['生产编号']
        specie = row['病原种拉丁名']

        
        ##原始长度比对NT数据库
        output_file1 = os.path.join(output_folder, f'{query_id}_blast_results1.txt')
        run_local_blast(sequence, db_path, output_file1)
        blast_results = parse_blast_results(output_file1)
        speci_blast_result = parse_speci_blast_results(output_file1,specie)

        df.at[idx, 'NT库比对结果'] = blast_results
        df.at[idx, '指定物种结果'] = speci_blast_result

        ##截取序列后比对NT数据库
        if len(sequence) > 60:  # 确保序列长度足够
            trimmed_sequence = sequence[40:-20]
        else:
            trimmed_sequence = sequence 
        output_file2 = os.path.join(output_folder, f'{query_id}_blast_results2.txt')
        run_local_blast(trimmed_sequence, db_path, output_file2)
        blast_results = parse_blast_results(output_file2)
        speci_blast_result = parse_speci_blast_results(output_file2,specie)

        df.at[idx, '截取产物序列'] = trimmed_sequence
        df.at[idx, '截取后NT库比对结果'] = blast_results
        df.at[idx, '截取后指定物种结果'] = speci_blast_result

    return df






if __name__ == "__main__":
    df1 = pd.read_excel(args.i, sheet_name="输入模板")
    df1 = get_path_and_seq(df1, args.p)
    db_path = "/data_tngs/common_database/genome/nt/nt"  # 本地 NT 数据库路径
    output_folder = "./blast_results"
    
    final_result = blast_sequences(df1, db_path, output_folder) 
    final_result.to_excel(args.o + ".xlsx",index=True)


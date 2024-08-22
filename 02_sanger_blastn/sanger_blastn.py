#!/usr/bin/env python
import pandas as pd
import os
# from Bio.Blast.Applications import NcbiblastnCommandline
# from Bio.Blast import NCBIXML
import re
import argparse
import sys
import subprocess

parser = argparse.ArgumentParser(
    description=f"Sanger产物序列blast到NT库，并整理分析比对结果\
    For example: python {sys.argv[0]} -i input.file -p sequence.path -o output -m max_target_seqs")
parser.add_argument("-i",type=str,
    help="输入模板文件，必须包含：生产编号两列")
parser.add_argument("-p",type=str,
    help="一代结果的存放目录，仅读取其中的seq文件")   
parser.add_argument("-m",type=str,default=100,
    help="指定blast的max_target_seqs参数,默认是100"),
parser.add_argument("-o",type=str,
    help="输出目录，不存在则自动创建")
parser.add_argument("-x",type=str,default=1,
    help="指定是否需要重新跑blast：1是需要；0是不需要，默认是1")
parser.add_argument("-s",type=str,
    help="指定配置文件：一代引物对应表.xlsx，会根据生产编号，索引出病原中文名和病原拉丁名")
args = parser.parse_args()

# args.i = "./muban2.xls"
# args.p = "./raw_sanger"
# args.m = 100
# args.o = "./test"
# args.x = "0"
# args.s = "./一代引物对应表.xlsx"



def get_path_and_seq(df,df_s, path):
    if '生产编号'  not in df.columns:
        raise ValueError("DataFrame 必须包含 '生产编号' 列")  
    
    df['病原种中文名'] = None
    df['病原种拉丁名'] = None
    
    for i, row in df.iterrows():
        production_number = row['引物名称']
        match = df_s[df_s['引物名'].str.contains(production_number, case=False, na=False)]
        
        if not match.empty:
            df.at[i, '病原种中文名'] = match.iloc[0]['病原']  
            df.at[i, '病原种拉丁名'] = match.iloc[0]['拉丁名']
    
    df['测序结果Seq文件名'] = ''
    df['产物序列'] = ''
    df['实际产物大小'] = 0

    num1 = 0
    lack_pid_list = []
    for idx, row in df.iterrows():  
        prod_id = str(row['生产编号'])
        
        # 查找包含生产编号的文件
        found_file = None
        for filename in os.listdir(path):
            if re.search(prod_id, filename) and filename.endswith('.seq'):
                found_file = filename
                break
        
        if found_file:
            file_path = os.path.join(path, found_file)
            num1 += 1
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
            lack_pid_list.append(prod_id)
    
    print("模板文件一共有记录:%d, 在结果文件夹找到Seq文件:%d, 未找到:%d, 为:%s"%(df.shape[0],num1,len(lack_pid_list),",".join(lack_pid_list)))
    return df



def run_local_blast(sequence, db_path, output_file):
    # 创建 BLAST 命令行: 使用subprocess 调用本地 blastn
    blast_command = [
        'blastn',
        '-query', sequence,
        '-db', db_path,
        '-evalue', '1e-5',
        '-outfmt', '6 qseqid sacc sscinames qlen qcovs pident evalue',  # 指定输出格式及字段
        '-out', output_file,
        '-max_target_seqs', str(args.m),
        '-num_threads', '20'
    ]
 
    try:
        subprocess.run(blast_command, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print("BLAST stderr:", e.stderr)
        print("BLAST returncode:", e.returncode) 


def parse_blast_results(df_raw,output_file,spec_row):
    # 解析 BLAST 输出：根据query ID进行分类解析
    df_all = pd.read_csv(output_file, sep='\t', header=None, names=[
        '生产编号','ACC号', '比对物种标题', '比对序列长度', '对比覆盖度', '对比相似度','Evalue'])
    
    # 过滤条件列表
    filter_terms = ['synthetic construct', '其他过滤条件1', '其他过滤条件2']
    # 过滤掉“比对物种标题”中包含任何过滤条件的行
    df_all = df_all[~df_all['比对物种标题'].str.contains('|'.join(filter_terms))] #~ 取反

    
    for idx, row in df_raw.iterrows():
        query_id = row['生产编号']
        df = df_all.loc[df_all['生产编号'] == query_id].copy()  

        if df.shape[0] == 0:
            df_raw.at[idx,spec_row] = "无比对产物"
        else : 
            df['pid_cov'] = (df['对比覆盖度'] * df['对比相似度'] / 100).round(2)
            group_stats = df.groupby('比对物种标题').agg(
                count=('ACC号', 'size'),
                max_pid_cov=('pid_cov', 'max')
            ).reset_index()
            
            total_count = group_stats['count'].sum()
            group_stats['percentage'] = (group_stats['count'] / total_count * 100).round().astype(int).astype(str) + '%'

            #先对group_stats按照max_pid_cov，后按照count进行降序
            group_stats_sorted = group_stats.sort_values(by=['max_pid_cov', 'count'], ascending=[False, False])
            def format_row(row):
                return f"{row['比对物种标题']}|count:{row['count']}({row['percentage']})|pid_cov:{row['max_pid_cov']}"
            formatted_rows = group_stats_sorted.apply(format_row, axis=1).tolist()
            # 填充到最多3项
            while len(formatted_rows) < 3:
                formatted_rows.append('NA')
            formatted_rows = formatted_rows[:3]
            result_string = '、'.join(formatted_rows)
            
            df_raw.at[idx,spec_row] = result_string  #提供正确的行索引和列名








def parse_speci_blast_results(df_raw,output_file,spec_row):
    df_all = pd.read_csv(output_file, sep='\t', header=None, names=[
    '生产编号','ACC号', '比对物种标题', '比对序列长度', '对比覆盖度', '对比相似度','Evalue'])

    for idx,row in df_raw.iterrows():
        query_id = row['生产编号']
        specie = row['病原种拉丁名']
        df = df_all.loc[df_all['生产编号'] == query_id]  

        if df.shape[0] == 0:
            df_raw.at[idx,spec_row] = "无比对产物"
        else:
            filtered_df = df[df['比对物种标题'].str.contains(specie, case=False, na=False)]
            if filtered_df.empty:
                cov = df.loc[df['Evalue'].idxmin(), '对比覆盖度']  # .idxmax() 返回具有最大值的索引位置
                ident = df.loc[df['Evalue'].idxmin(),"对比相似度"]
                result = f'非特异：未比对到该物种'
            else:
                cov = filtered_df.loc[filtered_df['Evalue'].idxmin(), '对比覆盖度']  # .idxmax() 返回具有最大值的索引位置
                ident = filtered_df.loc[filtered_df['Evalue'].idxmin(),"对比相似度"]
                result = f'Cov:{cov.round(2)}|Ident:{ident.round(2)}|特异'

            df_raw.at[idx,spec_row] = result 

            

def blast_sequences(df, db_path, output_folder):
    # 确保输出文件夹存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    #批量处理所有原始长度序列
    query_file_path = os.path.join(output_folder, 'queries.fasta')
    with open(query_file_path, 'w') as f:
        for idx, row in df.iterrows():
            sequence = row['产物序列']
            query_id = row['生产编号']
            f.write(f">{query_id}\n{sequence}\n")
    ##原始长度比对NT数据库
    output_file1 = os.path.join(output_folder, 'blast_results1.txt')
    if args.x == "1" :
        run_local_blast(query_file_path, db_path, output_file1)
    elif args.x == "0" and not os.path.exists(output_file1):
        print("请检出是否有跑好的blast结果文件，否则请修改 --x 参数为 1")
    elif args.x == "0" and os.path.exists(output_file1):
        pass
    parse_blast_results(df,output_file1,'NT库比对结果')  ##对输出的结果进行解析1
    parse_speci_blast_results(df,output_file1,"指定物种结果") ##对输出的结果进行解析1



    #批量处理裁剪之后的序列
    query_file_path2 = os.path.join(output_folder, 'queries2.fasta')
    with open(query_file_path2, 'w') as f:
        for idx, row in df.iterrows():
            sequence = row['产物序列']
            if len(sequence) > 60:  # 确保序列长度足够
                trimmed_sequence = sequence[40:-20]
            else:
                trimmed_sequence = sequence 
            query_id = row['生产编号']
            df.at[idx,"截取产物序列"] = trimmed_sequence
            f.write(f">{query_id}\n{trimmed_sequence}\n")
    ##裁剪序列比对NT数据库
    output_file2 = os.path.join(output_folder, 'blast_results2.txt')
    if args.x == "1" :
        run_local_blast(query_file_path2, db_path, output_file2)
    elif args.x == "0" and not os.path.exists(output_file1):
        print("请检出是否有跑好的blast结果文件，否则请修改 --x 参数为 1")
    elif args.x == "0" and os.path.exists(output_file1):
        pass
    parse_blast_results(df,output_file2,'截取后NT库比对结果')  
    parse_speci_blast_results(df,output_file2,"截取后指定物种结果") 

    return df   



if __name__ == "__main__":
    muban = args.i
    if not muban or not os.path.exists(muban):
        exit("错误，文件不存在:%s"%(muban))
    if muban.endswith("xlsx"):
        df1 = pd.read_excel(muban, sheet_name=0)
    else:
        df1 = pd.read_csv(muban,sep="\t")

    
    df_s = pd.read_excel(args.s,sheet_name = "Sheet3")
    df1 = get_path_and_seq(df1,df_s,args.p)
    db_path = "/data_tngs/common_database/genome/nt/nt"  # 本地 NT 数据库路径
    
    out_dir = os.path.abspath(args.o)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    output_folder = "%s/blast_results"%(out_dir)
    final_result = blast_sequences(df1, db_path, output_folder) 
    out_name = "%s/result_blast.xls"%(out_dir)
    final_result.to_csv(out_name,index=None,sep="\t")
    print("done! 结果文件:%s"%(out_name))
    
    


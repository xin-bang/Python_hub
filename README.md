# Python Hub

本项目是将日常工作中用到的有用的python代码作为记录，一来是方便以后遇到相同的问题能够直接调用；二来是加强对python代码的学习：

## 01_replace_fasta_name

```         
代码逻辑足够简洁，且处理序列相关的问题，仅使用到python的sys包，没有调用biopython（当然主要操作对象不是碱基序列）
```

#### 功能：

根据name.txt的信息，对序列文件的序列名称进行名称替换

-   输入文件1：name.txt：包含序列的原始名称、需要替换的名称
-   输入文件2：待处理的序列文件

#### 使用：

``` python
Usage: python replace_fasta_names.py <mapping_file> <fasta_file> <output_file
python3 replace_fasta_names.py name.txt 240708_TPMN00173_0339_A000H5VGV7_肺支_dTGM070_cluster.fasta 240708_TPMN00173_0339_A000H5VGV7_肺支_dTGM070_rename.fasta
```

## 02_sanger_blastn

```         
指定多个def module；调用本地nt库，使用Bio的接口进行blast运行，同时对blast结果进行整理操作。一些python的写法比较老道
```

#### 功能

-   根据输入文件索引相应的序列文件，
-   对相应的序列文件进行blast 比对，去100个比对结果，
-   从比对结果中进一步处理，按照相应的格式进行排序：

```         
例如：
    Adenovirus|count:12(23%)|pid_cov:97
    Adenovirus|count:12(8%)|pid_cov:96
    Adenovirus|count:12(2%)|pid_cov:90
< 按照比对HSP对应的物种拉丁名进行group_by；count为该物种的HSP占比情况；pid_cov为ident和cov的趁机，取最大值>
```

#### 使用

```         
usage: sanger_blastn.py [-h] [-i I] [-p P] [-m M] [-o O] [-x X] [-s S]

Sanger产物序列blast到NT库，并整理分析比对结果 For example: python sanger_blastn.py -i input.file -p sequence.path -o output -m max_target_seqs

optional arguments:
  -h, --help  show this help message and exit
  -i I        输入模板文件，必须包含：生产编号两列
  -p P        一代结果的存放目录，仅读取其中的seq文件
  -m M        指定blast的max_target_seqs参数,默认是100
  -o O        输出目录，不存在则自动创建
  -x X        指定是否需要重新跑blast：1是需要；0是不需要，默认是1
  -s S        指定配置文件：一代引物对应表.xlsx，会根据生产编号，索引出病原中文名和病原拉丁名
```

# Python Hub

本项目是将日常工作中用到的有用的python代码作为记录，一来是方便以后遇到相同的问题能够直接调用；二来是加强对python代码的学习：

## 01_replace_fasta_name

    代码逻辑足够简洁，且处理序列相关的问题，仅使用到python的sys包，没有调用biopython（当然主要操作对象不是碱基序列）

#### 功能：

根据name.txt的信息，对序列文件的序列名称进行名称替换

* 输入文件1：name.txt：包含序列的原始名称、需要替换的名称
* 输入文件2：待处理的序列文件

#### 使用：

```python
Usage: python replace_fasta_names.py <mapping_file> <fasta_file> <output_file
python3 replace_fasta_names.py name.txt 240708_TPMN00173_0339_A000H5VGV7_肺支_dTGM070_cluster.fasta 240708_TPMN00173_0339_A000H5VGV7_肺支_dTGM070_rename.fasta
```

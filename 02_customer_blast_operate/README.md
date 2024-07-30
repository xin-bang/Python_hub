custome_blast_operate.py

## 说明：RP8的IVD项目：

#### 输入：
指定输入文件（实验同事写入）

#### 基本功能：
* 根据输入文件索引相应的序列文件，
*  对相应的序列文件进行blast 比对，去100个比对结果，
*  从比对结果中进一步处理，按照相应的格式进行排序：
```
例如：
    Adenovirus|count:12(23%)|pid_cov:97
    Adenovirus|count:12(8%)|pid_cov:96
    Adenovirus|count:12(2%)|pid_cov:90
< 按照比对HSP对应的物种拉丁名进行group_by；count为该物种的HSP占比情况；pid_cov为ident和cov的趁机，取最大值>
```

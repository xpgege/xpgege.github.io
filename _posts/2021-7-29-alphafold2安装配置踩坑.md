---
layout:     post
title:      alphafold2安装配置踩坑
subtitle:   alphafold2
date:       2021-7-29
author:     xpgege
header-img: img/post-web.jpg
catalog: true
tags:
   - 知识记录
---

​	2021年7月15日,deepmind在CASP14中大获全胜的alphafold2终于开源了，Demis Hassabis、John Jumper等人在Nature杂志上发表了文章[Highly accurate protein structure prediction with AlphaFold](https://www.nature.com/articles/s41586-021-03819-2)，这是一个基于神经网络的新模型，其预测的蛋白质结构能达到原子水平的准确度。

![CASP14 predictions](https://github.com/kuixu/alphafold/raw/main/imgs/casp14_predictions.gif)

Deepmind博客提供的结果，绿色表示实验结果，蓝色表示预测结果，在考虑到实验自身误差的前提下，alphafold预测的结果与实验结果相差无几。

**github**：https://github.com/deepmind/alphafold

------

下面我记录一下我安装alphafold2的踩坑，具体论文解析，我将后面在进一步了解



# 安装流程

 **slightly simplified version of AlphaFold**

可以进入注册登陆谷歌[colab](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb)，直接网页端进行预测

![截图录屏_选择区域_20210729151713](https://i.loli.net/2021/07/29/8vUygOhpKcDHuCP.png)

------



**whole version of AlphaFold**

1. Install [Docker](https://www.docker.com/).

   - Install [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) for GPU support.
   - Setup running [Docker as a non-root user](https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user).

2. Download genetic databases (see below).

3. Download model parameters (see below).

4. Check that AlphaFold will be able to use a GPU by running:

   ```
   docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
   ```

   The output of this command should show a list of your GPUs. If it doesn't, check if you followed all steps correctly when setting up the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) or take a look at the following [NVIDIA Docker issue](https://github.com/NVIDIA/nvidia-docker/issues/1447#issuecomment-801479573).



上面为官方提供的流程

------



## 1.git clone

```
git clone https://github.com/deepmind/alphafold
unzip alphafold-main.zip
cd alphafold-main
```

## 2.下载database

![截图录屏_选择区域_20210729152350](https://i.loli.net/2021/07/29/tiRev1oEfIbWj5G.png)

 **The total download size for the full databases is around 415 GB and the total size when unzipped is 2.2 TB.**

官方提供script脚本用于下载数据库，提前安装```aria2c ```和```rsync```

<img src="https://i.loli.net/2021/07/29/jyOZa3eSAkxvRId.png" alt="截图录屏_选择区域_20210729152957" style="zoom:150%;" />

```shell
mkdir database
sh scripts scripts/download_all_data.sh database
```

aria2c会创建一个临时文件，占据目标下载文件的空间，如mgy_clusters_2018_12.fa.gz，它会直接占据32.9G空间，然后再进行下载。

由于不知是否是国内网速问题，aria2c下载速度实在过慢，仅几百kb，而且十分不稳定



故我用**axel**代替**aria2c**，axel是一种多线程下载工具

将所有脚本中

```
aria2c "${SOURCE_URL}" --dir="${ROOT_DIR}"
替换成
axel -n 20 -N "${SOURCE_URL}" -o "${ROOT_DIR}"
```

axel 下载参数如https://blog.csdn.net/yangl2512/article/details/6833547



**我花了两天下载，主要是rsync同步速度过慢，而且Running rsync to fetch all mmCIF files (note that the rsync progress estimate might be inaccurate)...，这不是不准确，是相当不准确，第一天晚上我看到同步了99%，知道第二天晚上变成98%了.....**



## 3.docker 创建环境

```
docker build -f docker/Dockerfile -t alphafold .
#-f 指定dockerfile 路径
#-t 指定docker image的label
```

然而这一步需要2.41T空间来储存docker image，而docker 镜像一般储存在/var/lib/docker/，空间不够，所以要将docker 镜像储存目录迁移到更大的磁盘上。

参考：[docker磁盘空间不足解决办法 - 骑马挎枪打天下 - 博客园](https://www.cnblogs.com/linux123/p/12176784.html)

但是这一次我没有成功..（以前我做的时候成功了）

![截图录屏_选择区域_20210726160825](https://i.loli.net/2021/07/29/lxyvLzbGUHtBmEj.png)

因为部分在docker/lib/overlay2中的文件我无法操作，即使我sudo/su root

尝试了很久都找不到准确的答案

可能原因：

1.挂载/mnt/storage磁盘格式的问题

2.docker镜像内部权限的问题，外部无法访问内部文件

（希望知道的朋友能指点一二）



## 3*.no docker，but conda

在alphafold github issue中发现kuixu 不利用docker 的alphafold改进版本

https://github.com/kuixu/alphafold

```
git clone https://github.com/kuixu/alphafold
conda create -n af2 python=3.8 -y
conda activate af2

conda install -y -c nvidia cudnn==8.0.4
conda install -y -c bioconda hmmer hhsuite==3.3.0 kalign2

conda install -y -c conda-forge \
    openmm=7.5.1 \
    pdbfixer \
    pip

# python pkgs
pip3 install --upgrade pip \
    && pip3 install -r ./requirements.txt \
    && pip3 install --upgrade "jax[cuda111]" -f \
    https://storage.googleapis.com/jax-releases/jax_releases.html

# work_path=/path/to/alphafold-code
work_path=$(PWD)

# update openmm 
a=$(which python)
cd $(dirname $(dirname $a))/lib/python3.8/site-packages
patch -p0 < $work_path/docker/openmm.patch

```

python  及cuda对应版本如图所，示请根据自己gpu设备情况对 cudnn tensorflow版本 进行修改

![截图录屏_选择区域_20210729161407](https://i.loli.net/2021/07/29/r7ITlmXU8jqBs1M.png)





（其中**jax** 是google的一种新的AI框架，结合了tensorflow和pytorch的优点，如结合了gpu加速的numpy ，同样的数学运算，用JAX版的numpy可以加快30多倍，因为使用了CUDA GPU加速）

jax github:https://github.com/google/jax#installation

```
pip install --upgrade "jax[cuda111]" -f https://storage.googleapis.com/jax-releases/jax_releases.html

#For CUDA 11.1, 11.2, or 11.3, use cuda111. The same wheel should work for CUDA 11.x releases from 11.1 onwards.
For CUDA 11.0, use cuda110.
For CUDA 10.2, use cuda102.
For CUDA 10.1, use cuda101.
Older CUDA versions are not supported.
（我们low b的实验室还是1080Ti，我自己电脑显卡都2060Ti了，害，让老板买显卡去了）
```



## 4.运行alphafold

```
vim run_alphafold.py

修改如下内容
# Set to target of scripts/download_all_databases.sh
DOWNLOAD_DIR = '/path/to/database'

# Path to a directory that will store the results.
output_dir = '/path/to/output_dir'

uniclust30_database_path = os.path.join(
    #DOWNLOAD_DIR, 'uniclust30', 'UniRef30_2020_02','UniRef30_2020_02')
    DOWNLOAD_DIR, 'uniclust30', 'uniclust30_2018_08', 'uniclust30_2018_08')

```

在[AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)选了个长度较短的序列

```
wget https://www.uniprot.org/uniprot/P26647.fasta
python3 run_alphafold.py --fasta_paths=P26647.fasta --max_template_date=2020-05-14
# or simply
exp/run_local.sh P26647.fasta
```

**参数**：

```
	  'fasta_paths',
      'output_dir',
      'model_names',
      'data_dir',
      'preset', 
      'uniref90_database_path',
      'mgnify_database_path',
      'uniclust30_database_path',
      'bfd_database_path',
      'pdb70_database_path',
      'template_mmcif_dir',
      'max_template_date',
      'obsolete_pdbs_path',
```

control AlphaFold speed / quality tradeoff by adding either `--preset=full_dbs` or `--preset=casp14` to the run command. We provide the following presets:

- **casp14**: This preset uses the same settings as were used in CASP14. It runs with all genetic databases and with 8 ensemblings.（质量）

- **full_dbs**: The model in this preset is 8 times faster than the `casp14` preset with a very minor quality drop (-0.1 average GDT drop on CASP14 domains). It runs with all genetic databases and with no ensembling.（速度）

  

**输出：**

```
output_dir/
    features.pkl
    ranked_{0,1,2,3,4}.pdb
    ranking_debug.json
    relaxed_model_{1,2,3,4,5}.pdb
    result_model_{1,2,3,4,5}.pkl
    timings.json
    unrelaxed_model_{1,2,3,4,5}.pdb
    msas/
        bfd_uniclust_hits.a3m
        mgnify_hits.sto
        uniref90_hits.sto


```



- `features.pkl` – A `pickle` file containing the input feature Numpy arrays used by the models to produce the structures.
- `unrelaxed_model_*.pdb` – A PDB format text file containing the predicted structure, exactly as outputted by the model.
- `relaxed_model_*.pdb` – A PDB format text file containing the predicted structure, after performing an Amber relaxation procedure on the unrelaxed structure prediction, see Jumper et al. 2021, Suppl. Methods 1.8.6 for details.
- `ranked_*.pdb` – A PDB format text file containing the relaxed predicted structures, after reordering by model confidence. Here `ranked_0.pdb` should contain the prediction with the highest confidence, and `ranked_4.pdb` the prediction with the lowest confidence. To rank model confidence, we use predicted LDDT (pLDDT), see Jumper et al. 2021, Suppl. Methods 1.9.6 for details.
- `ranking_debug.json` – A JSON format text file containing the pLDDT values used to perform the model ranking, and a mapping back to the original model names.
- `timings.json` – A JSON format text file containing the times taken to run each section of the AlphaFold pipeline.
- `msas/` - A directory containing the files describing the various genetic tool hits that were used to construct the input MSA.
- `result_model_*.pkl` – A `pickle` file containing a nested dictionary of the various Numpy arrays directly produced by the model. In addition to the output of the structure module, this includes auxiliary outputs such as distograms and pLDDT scores. If using the pTM models then the pTM logits will also be contained in this file.




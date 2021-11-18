---
layout:     post
title:      Alphafold文献解析
subtitle:   Alphafold文献解析
date:       2021-4-10
author:     xpgege
header-img: img/post-web.jpg
catalog: true
tags:
    - 文献解析
---

# ***\*AlphaFold\****

 

长期以来，科学家一直对确定蛋白质的结构感兴趣，因为人们认为蛋白质的结构决定了其功能。一旦了解了蛋白质的形状，就可以猜测其在细胞中的作用，进而可以开发出与蛋白质的独特形状有关的药物。

![img](https://i.loli.net/2021/06/26/yMEJKYpO7zbNdPj.jpg) 

在过去的五十年中，研究人员已经能够使用诸如冷冻电子显微镜，核磁共振和X射线晶体学这样的实验技术在实验室中确定蛋白质的形状，但是每种方法都取决于大量的试验和错误，这可能需要花费多年的时间，并且每个蛋白质结构的成本高达数万或数十万美元。

这就是为什么生物学家想利用AI方法来替代这一艰难而艰苦的蛋白质过程的原因。不通过昂贵的实验来确定蛋白质，仅凭其序列即可预测蛋白质形状的能力，可以帮助加速研究。



# **1.** ***\*背景介绍\****

## **1.1.** ***\*氨基酸\****

 ![image-20211102114930640](https://i.loli.net/2021/11/02/9okln57vDmAB8RV.png)

根据高中生物知识我们知道，人体蛋白质由20种氨基酸，根据侧链性质的不同大致可以分为4类，带电（正电，负电氨基酸，不带电极性氨基酸，特殊氨基酸，疏水氨基酸，例如正电负电氨基酸会相互吸引和排至，疏水氨基酸残基会处在分子内部形成疏水[内核](https://baike.baidu.com/item/内核)，从而维系蛋白质的紧密三维结构。正是这些侧链化学性质不同，导致蛋白质三维结构的多样性。

这是任意一个氨基酸的基本结构。

![img](https://i.loli.net/2021/06/26/1HvhJLORY9rzCuS.jpg) 

定义：羧基作为氨基酸的主要官能团，与羧基C直接相连的C原子为α-C，与α-C直接相连的C原子为β-C（甘氨酸 R基为H没有β-C）



## **1.2.** ***\*蛋白质\****

一级结构：氨基酸残基的排列顺序

 

![img](https://i.loli.net/2021/06/26/UZ72hStO4fwsWCP.jpg) 

蛋白质的一级结构可以理解为一条线性的字符串，比如MSFIKTFSGKHFYYDKINKDDIVINDIAVSLSNICR。其基本组成单元是一个个的氨基酸，即一个个的字母。氨基酸有单字母表示和三字母表示，为了简洁，本文使用单字母表示，下图的例子是三字母表示。常见的氨基酸只有20种，所以一级结构的字符串通常只包含20种字母

![img](https://i.loli.net/2021/06/26/1BZb36ljHQU4kDC.jpg) 

20种氨基酸连接的方式为脱水缩合，即一个氨基酸的羧基-COOH和另一个氨基酸的氨基-NH2反应，丢掉一个H2O，形成一个肽键-CO-NH-，如下图所示。丢掉了羧基和氨基的氨基酸被称为氨基酸残基，这个名词很形象，氨基酸缺胳膊少腿，所以变成了“残”基。

 

二级结构：二级结构就是在一级结构的字符串的基础上，肽链怎样进行盘旋、折叠等变换，形成一种***\*局部\****的三维结构，这种局部的三维结构通常由氢键支撑。如α-helix，β-sheet，β-turn，random coil等八种结构

其中α螺旋的每个残基的-NH的H和临近的第4个残基的-CO的O形成氢键，由此支撑α螺旋的结构稳定性，如下图的箭头所指虚线

![img](https://i.loli.net/2021/06/26/XZilEkctNb8S2mK.jpg) 

β折叠则是两条肽链，平行排列，对应残基的-NH的H和-CO的O形成氢键，由此形成两股β折叠的结构，多股β折叠形成类似手风琴的样子。

![img](https://i.loli.net/2021/06/26/QJ84RVI5hxyWuZK.jpg) 

结构域( domain)：位于[超二级结构](https://baike.baidu.com/item/超二级结构/1479765)和三级结构间的一个层次，由100~ 200个氨基酸残基组成，各有独特的空间构象，通常domain会作为蛋白质的活性中心，便于与其他分子之间的相互作用

![img](https://i.loli.net/2021/06/26/izK7ukUotYn4q9a.jpg) 

三级结构：简单理解，三级结构就是把多个二级结构拼接到一起，折叠成一个完整的蛋白质三维结构，如下图所示。维持蛋白质三级结构的力比较多样，除了氢键之外，还有二硫键、金属键等。

![img](https://i.loli.net/2021/06/26/5Xgf2VMbkDImen7.jpg)![img](https://i.loli.net/2021/06/26/GWYjIKw4BZkE9ui.png)

![img](https://i.loli.net/2021/06/26/yKa62JsXEBktVFO.jpg) 

## **1.3.** 蛋白质结构预测问题描述

所以我们的输入是蛋白质一级结构序列

![img](file:////tmp/wps-xpgege/ksohtml/wpsiVwtKY.jpg) 

![img](file:////tmp/wps-xpgege/ksohtml/wpsjmxlJm.jpg) 

目标是蛋白质的三级结构，为蛋白质中每个原子指定3D坐标，对应于生物物理上合理的构象。

 

所以这是生物化学的基本挑战，可以说是最困难的预测任务之一，但目前蛋白质数据库已经积累了大量数据供以使用，但只有很少一部分已只结构的蛋白质储存在PDB数据库中，也就是深度学习中的label data，不仅已知序列，而且已知蛋白质分子中每个原子的三维坐标，而绝大多数蛋白质是仅已知序列的unlabel data。

![img](file:////tmp/wps-xpgege/ksohtml/wpsOFUfIK.jpg) 

那剩下大量的unlabel蛋白质数据我们应该如何使用呢？后面我将会讲到

 

 

## **1.4.** ***\*CASP\****

蛋白质三级结构预测，不得不提的是CASP这个比赛。CASP的全称是The Critical Assessment of protein Structure Prediction (CASP)，即蛋白质结构预测的关键评估，被誉为蛋白质结构预测的奥林匹克竞赛。CASP从1994年开始举办，每两年一届，最近的一届是2018年的CASP13。

每一届CASP比赛，都会提供大约100条未知结构的蛋白质序列，让所有参赛者进行结构预测，比赛结束之后，主办方会通过生化方法测定这些蛋白质的三维结构，然后和参赛者预测的结果进行比对，然后给出预测得分。

提供的蛋白质序列分为两类：一类序列和PDB数据库中已有结构的序列有相似性，由此可以基于模板预测，准确度比较高，这类算法称为Template-Based Modeling；另一类序列和PDB库已知结构的序列相似度很低，可以认为是全新的蛋白质，因为无法利用已有模板信息，需要进行从头测序（Free Modeling），目前的准确率相对较低。

CASP同时提供多种比赛项目，比如常规的结构预测（Regular targets）、数据辅助预测（Data-Assisted targets）和蛋白质接触面预测（Contact predictions）等，其中数据辅助预测中提供了核磁数据（NMR）、交联数据（XLMS）等

![img](file:////tmp/wps-xpgege/ksohtml/wpsvCfdH8.jpg) 

而 Deepmind（也就是谷歌的当年alphago设计团队）2018年在CASP13中提出了Alphafold参加第一次参赛就一鸣惊人，获得冠军并甩开第二名好多。

![img](file:////tmp/wps-xpgege/ksohtml/wps4NPcGw.jpg) 

而在2020年CASP14中，Deepmind再次提出Alphafold2，不仅取得了第一，远远超过他们之前的alphafold而且结构准确度甚至能达到90%，基本与实验水平，相差无几。下图是alphafold2的预测结果，其中绿色的为实验结果，蓝色的为计算结果，我们可以看到二者几乎保持overlap，说明alphafold的预测精度已经达到很高的水平。



由于Alphafold2论文还未公开，下面我仅对Alphafold1中的模型进行解析。

![img](https://i.loli.net/2021/06/26/CsWoteJ7d1mEO2Q.jpg) 

![img](https://i.loli.net/2021/06/26/Rr4EAaWLFZNpygu.png) 

 

论文及代码链接如下

代码：https://github.com/deepmind/deepmind-research/tree/master/alphafold_casp13

论文：

https://www.biorxiv.org/content/10.1101/846279v1.full.pdf

https://doi.org/10.1002/prot.25834

# **2.** ***\*Alphafold结构\****

了解了基本的生物学背景后，我们先来简单看看alphafold 的基本结构

结构主要分为3大块

1.对序列进行多序列比对和特征提取

2.预测预测成对氨基酸残基之间的距离分布

3.预测肽链骨架的扭角分布

4..对扭脚和距离分布生成的candidate structure进行模型评估

（因为自然界总有势能减小的趋势，选择势能最低的结构作为最后预测结果）

所以需要构建三个神经网络分别预测distance，torsion和candidate structure

 

## **2.1.** ***\*多序列比对和特征提取\****

前面提到，大量的unlabel蛋白质数据我们应该如何使用呢？我们可以利用生物信息学中一个常用工具MSA Multiple sequence alignment 多序列比对，将一个quary senquence S 输入到所有生物的蛋白质序列数据库中进行多序列比对，找到与之对应的相似序列，称作同源序列（可能具有共同的基因祖先）

![img](https://i.loli.net/2021/06/26/NLzh9fkVeBZCxKd.jpg) 

Reference：Improved contact prediction in proteins: Using pseudolikelihoods to infer Potts models

 

 

所以有个最简单的情况 如果这个相似序列的蛋白质结构已知了，我们将它的结构迁移到新的蛋白质结构预测中来，也就是Template-Based Modeling

而对于Free Modeling，我们也可以获得它们之间的进化关系。生物大分子有个特性，就是它们在进化中是保守的，因为生物体内无时无刻不在发生着突变，如果突变导致蛋白质的结构发生剧烈的改变，这对生物体来说是个灾难，所以对于大多数蛋白质，他都会保持功能的稳定性，比如图中的两个点发生了接触，如果S突变为F，为了保持蛋白质功能的稳定，在生物进化过程中H也会突变为W，使得这两个位点始终保持接触

所以我们可以把多序列比对的结果也作为特征，加入到序列原本的特征中进行训练，从而利用到unlabel蛋白质提供的信息。

 

![img](https://i.loli.net/2021/06/26/fdLBbxXJ4CmocHp.jpg) 

 

数据：PDB中的蛋白质结构数据，通过CATH的35%蛋白质序列相似度进行聚类和去冗余，	 得到31247条domain，其中 train 29,427  and  test sets 1,820 

利用HHblits和PSI-Blast两种多序列比对工具，与 Uniclust30数据库进行比对

特征：

Number of HHblits alignments (scalar).

 

Sequence-length features: 

1-hot amino acid type (21 features); 

PSI-BLAST (21 features),

HHblits profile (22 features), 

non-gapped profile (21 features)

HHblits bias, HMM profile (30 features)

Potts model bias (22 features)

deletion probability (1 feature)

residue index

 

Sequence-length-squared features: Potts model parameters (484 features)

Frobenius norm (1 feature)

gap matrix (1 feature)

 

 

## **2.2.** ***\*距离分布预测\****

### **2.2.1.** ***\*预测方法\****

假设有n条蛋白质，每条序列长度为L（每个蛋白质长度L不同），将序列进行多序列比对，将多序列比对结果和序列自身信息，有m个FEATURE，得到n*L*m的矩阵

将每个FEATURE之间作协相关得到n*L*L*m的矩阵作为input

后面神经网络的处理就和image-recognition有点相似

将这个L*L的sequence矩阵（每个蛋白质长度L不同），分割成若干个64*64的不重合的区块如下图所示，既能降低内存占用，也能起到data argumentation，降低过拟合的作用。

Alphafold对每个sequence进行四次分割，分别输入到四个超参略微不同的神经网络中，最后对结果进行加权平均。

![img](https://i.loli.net/2021/06/26/sdNwaFT5Hbr7fuA.jpg)![img](https://i.loli.net/2021/06/26/bZPs2Yo3X5fNgay.jpg) 

[（https://predictioncenter.org/casp13/doc/presentations/Pred_CASP13-DeepLearning-AlphaFold-Senior.pdf）](https://predictioncenter.org/casp13/doc/presentations/Pred_CASP13-DeepLearning-AlphaFold-Senior.pdf)

输出的数据为任意两个氨基酸之间的Cβ原子之间的距离（range 2–22 Å）的概率，将这个距离等分成64份，每个bin就是一个小区间，所以最终结果转化为一个64分类问题，对应图中最后的结果是一个n*L*L*64的一个矩阵，每一个样本点就是，在该位置成对氨基酸之间距离在该区间的概率分布。

 

 

![img](https://i.loli.net/2021/06/26/IKhpSbR1AnFJfzW.jpg) 

举个简单例子

我们来看这个蛋白质与29号残基的距离，假设截取的氨基酸的序列为S和FYPDTWL  根据原始数据坐标信息，算得他们之间的距离为S与FYPDTWL直接的距离分别为1，2，3，4，5，6，7Å

为了方便说明，假设1-7Å分成8份

每个区间为[0-1),[1-2),....[6-7)，[7-8)

S与F，1Å 在[1-2)区间里，S与Y，2Å在 [2-3)里.....

故y_label S-F 为[0,1,0,0,0,0,0,0]，y_label S-Y为[0,0,1,0,0,0,0,0]

Y_predS-F可能为[0.1,0.5,0.1,0.1,0.1,0.05,0.05,0]，[0.1,0,0.7,0.05,0.05,0.05,0.05,0]

故Y_pred 即为预测的离散概率分布

同理，每一对氨基酸都对应一个点，同时对应一个离散概率分布

就有n*L*L个概率分布，故Y.shape=n*L*L*64

 

[（https://predictioncenter.org/casp13/doc/presentations/Pred_CASP13-DeepLearning-AlphaFold-Senior.pdf）](https://predictioncenter.org/casp13/doc/presentations/Pred_CASP13-DeepLearning-AlphaFold-Senior.pdf)

我们来看这个蛋白质中1-41号氨基酸残基与29号氨基酸残基之间的预测情况。

画出对应的分布直方图，由于这个分布是不连续的，如果有同学选修了数值分析，应该知道，可以利用cubic spline（三次样条插值）对曲线进行拟合，得到右图。

![img](https://i.loli.net/2021/06/26/2RiHMycuvJ73qFW.jpg)![img](https://i.loli.net/2021/06/26/cZ9N37nAH1UaOiv.jpg) 

 

然后对拟合的曲线进行积分，得到平均距离d*，如果d*<8Å,我们可以认为这两个氨基酸相互接触。左图中，红色的线表示，两个氨基酸残基实际距离，黑线表示8Å，我们可以看到在29号氨基酸附近，红线基本处于分布图最高点附近，分布都预测的比较好，而离29号较远的分布预测稍差。

### **2.2.2.** ***\*神经网络结构\****

了解距离分布预测的流程后，我们再来回过头看神经网络的设计

普通卷积

![img](https://i.loli.net/2021/06/26/Cu7FfvlEcBxdmAw.png)![img](https://i.loli.net/2021/06/26/BmPTlzCSxq19iXN.jpg) 

l 卷积层：dilated conv，即空洞卷积

 

Reference：https://www.jianshu.com/p/f743bd9041b3

![img](https://i.loli.net/2021/06/26/4m1pdKNOGbxRSjr.jpg) 

(a) 图对应3x3的1-dilated conv，和普通的卷积操作一样，

(b) 图对应3x3的2-dilated conv，实际的卷积kernel size还是3x3，但是空洞为1，也就是对于一个7x7的图像patch，只有9个红色的点和3x3的kernel发生卷积操作，其余的点略过。也可以理解为kernel的size为7x7，但是只有图中的9个点的权重不为0，其余都为0。 可以看到虽然kernel size只有3x3，但是这个卷积的感受野已经增大到了7x7（如果考虑到这个2-dilated conv的前一层是一个1-dilated conv的话，那么每个红点就是1-dilated的卷积输出，所以感受野为3x3，所以1-dilated和2-dilated合起来就能达到7x7的conv）

(c) 图是4-dilated conv操作，同理跟在两个1-dilated和2-dilated conv的后面，能达到15x15的感受野。

(d) 对比传统的conv操作，3层3x3的卷积加起来，stride为1的话，只能达到(kernel-1)*layer+1=7的感受野，也就是和层数layer成线性关系，而dilated conv的感受野是指数级的增长。

dilated的好处是不做pooling损失信息的情况下，加大了感受野，让每个卷积输出都包含较大范围的信息。

 

l 激活函数：ELU![img](https://i.loli.net/2021/06/26/7SJ2rnjZoex94HA.jpg)

![img](https://i.loli.net/2021/06/26/Tz6SU5BYueR1Wwa.jpg) 

 

融合了sigmoid和ReLU，左侧具有软饱和性，右侧无饱和性，右侧线性部分使得ELU能够缓解梯度消失，而左侧软饱能够让ELU对输入变化或噪声更鲁棒。

ELU的输出均值接近于零，所以收敛速度更快。

 

l 结构图

 

![img](https://i.loli.net/2021/06/26/ZdJqmVoU5aefPb1.png)![img](https://i.loli.net/2021/06/26/olxNhOGUDjRQP9a.jpg)          ![img](https://i.loli.net/2021/06/26/m2pfyCuWUIEdV6J.jpg)

整个神经网络每个层都由左图的残差卷积块构成，每个残差卷积块的Dilation依次为1，2，4，8循环，共220层网络结构

 

l 超参数

• 7 groups of 4 blocks with 256 channels, cycling through dilations 

1, 2, 4, 8.

• 48 groups of 4 blocks with 128 channels, cycling through dilations 

1, 2, 4, 8.

• Optimization: synchronized stochastic gradient descent

• Batch size: batch of 4 crops on each of 8 GPU workers.

• 0.85 dropout keep probability.

• Nonlinearity: ELU.

• Learning rate: 0.06.

• Auxiliary loss weights: secondary structure: 0.005; accessible surface area: 0.001. These auxiliary losses were cut by a factor 10 after 

100 000 steps.

• Learning rate decayed by 50% at 150,000, 200,000, 250,000 and 

350,000 steps.

• Training time: about 5 days for 600,000 steps.



## **2.3.** ***\*扭角分布预测\****

### **2.3.1.** ***\*扭脚\****

我先来讲解一下什么是肽链的扭脚，我们可以看到每个氨基酸残基有三个骨架原子，一个羧基C，一个αC和一个N原子，由于肽键作为一个酰胺键，键长介于单键与双键之间，拥有一部分双键的性质，不能自由转动，所以就和两边两个αC共六个原子构成了肽的基本单元，肽平面，也叫酰胺平面。（参考视频：[Ramachandran Principle: Protein Atomic Clashes vs. Phi & Psi - YouTube](https://www.youtube.com/watch?v=Q1ftYq13XKk)）

所以我们只要找到每个氨基酸残基的扭角，就可以确定一个蛋白质骨架的结构。而αC和羧基C相连的单键的扭角称为φ, αC和氨基N相连的单键的扭角为ψ，所以长度为L的蛋白质有2L个扭角。

![img](https://i.loli.net/2021/06/26/DTroJsI6XQzuwSW.png)![img](https://i.loli.net/2021/06/26/5agyexlOXK87T6D.jpg) 

图片来源：[Image:Dihedral-angles-anim.gif - Proteopedia, life in 3D](https://proteopedia.org/wiki/index.php/Image:Dihedral-angles-anim.gif)

那我们怎么求这个扭角呢？

这个扭角其实也就是高中所学的二面角，我们训练数据已知每个原子的坐标，例如求下图中的phi，只要构建两个平面，123号原子一个平面，234号原子一个平面，求出两个平面的法向量m、n，然后就可以得到每两个氨基酸残基之间扭角。

![img](https://i.loli.net/2021/06/26/Vs6871Gh9FdNR4E.jpg) 



### **2.3.2.** ***\*神经网络结构\****

Alphafold使用变分自编码器（参考https://zhuanlan.zhihu.com/p/351805989）端到端训练骨架扭脚。结构参照Deepmind的另一篇文章[DRAW: A Recurrent Neural Network For Image Generation](https://arxiv.org/abs/1502.04623)，DRAW网络结合了模仿人眼偏好注意力机制可迭代的构造复杂图像，进行顺序自动编码，如下图所示，在MNIST数据集中，和传统的VAEs网络不同的是，图像的生成不是一次性全部生成，而是随着时间仅仅关注红色矩形框内的pixel，使得每一部分的图像的生成是相互独立的，这和人类数字的书写习惯基本保持一致。

![img](https://i.loli.net/2021/06/26/gUeEKmY4FVBaNiH.jpg)![img](https://i.loli.net/2021/06/26/YsXlTApZn12DxPK.jpg) 

与这项工作相似，Alphafold也不是一次性预测整个肽链的扭脚的分布，而是crop成多个32 size的短肽，feed到encoder中估计每个扭角正态分布的latent parameter（two per angle：μ，σ），对每个正态分布进行采样，feed到decoder中解码得到每个canvas，最后结合起来得到φ,ψ的von Mises distributions冯·米塞斯分布。

（参考：https://blog.csdn.net/weixin_42576437/article/details/107124117）

与DRAW关键区别在于，这里需要根据蛋白质sequence和MSA得到的condition vector作为条件概率调整encoder和decoder。

 

![img](https://i.loli.net/2021/06/26/rhAkTQoIVugbxcK.jpg) 

## **2.4.** ***\*模型评估\****

虽然我们上面已经得到了，氨基酸残基之间的距离分布和扭角分布，但是由于这只是一个概率分布，而不是一个确定的值，所以Alphafold 定义了一个protein-specific potential，也就是一个蛋白质的势能，势能越低，而对应蛋白质的三维构像越稳定

 

定义：x=G(φ,ψ), X为氨基酸残基βC的坐标，dij=||Xi-Xj||

 

工具：Rosetta  Rosetta是基于蒙特卡罗模拟退火为算法核心的高分子建模软件库，通过Rosetta内置的打分函数以及各类方法来采样生物过程，评估和优化这些高分子结构

Reference：https://zhuanlan.zhihu.com/p/65248904

 

l Distance potential

 预测的距离概率分布的势能

![img](https://i.loli.net/2021/06/26/W4vqFyCMfe3AdRV.jpg) 

背景概率分布的势能

(也就是统计上两个氨基酸距离为length，氨基酸种类分别为AB时，距离为dij的概率)

![img](https://i.loli.net/2021/06/26/BncpVrwvtj3xa9h.jpg) 

l geometric potential

![img](https://i.loli.net/2021/06/26/3gmfeAiqWZv2GIX.jpg) 

l smooth potential

由于前面的仅仅预测链骨架的结构，缺乏对侧链的考虑，例如图中在肽链旋转的过程中，侧链会发生立体碰撞，使得原子间距离过进，范德华力突变为排斥力，使得分子势能增加，这种构象是不稳定的，所以我们引入a van der Waals term 

![img](https://i.loli.net/2021/06/26/3xdYuf7nozs26BJ.png)![img](https://i.loli.net/2021/06/26/4FogWbhA9lSntKc.png) 

Reference：https://proteopedia.org/wiki/index.php/Dihedral/Index

  

总的势能函数可以写为

![img](https://i.loli.net/2021/06/26/q6hLo9jvFD2wWHX.jpg) 

（这三个势能的通过交叉验证，在相同weight的情况下有最好的效果）

然后利用Rosetta对分子构象进行打分

所以整个Vtotal势能可看做φ,ψ的函数，通过随机采样torsion初始化φ,ψ，利用L-BFGS梯度下降minimize Vtotal。

由于整个优化的结构非常以依赖初始化的情况，所以在structure pool里保留20个具有最低势能的结构，使得初始化的φ,ψ 90%来自 structure poool，10%来自分布采样，再加上30°的噪声

# **3.** ***\*总结\****

## **3.1.** ***\*Alphafold\****

所以整个模型算法流程可以看做如下

![img](https://i.loli.net/2021/06/26/RbrlPpVBezc26Fq.jpg) 

将带预测的蛋白质序列数据与蛋白质数据库进行多序列比对，提取对应的MSA特征并与序列数据一起输入深度残差网络，预测扭角(2*L)分布和成对残基(L^2)之间的距离分布，将扭角，距离，范德华力合并作为待优化的势能，从扭角分布初始化φ,ψ，进行L-BFGS梯度下降，将产生的低势能结构，扭角采样，noise加权平均，更新初始化的φ,ψ，反复迭代直至结构收敛。

![img](https://i.loli.net/2021/06/26/QbaVq8SCT1km6zF.jpg) 

从上图我们可以看到在反复迭代过程中，rmsd不断下降，TM(template model)score不断上升，整个蛋白质构象也趋于稳定。

 

但是尽管在CASP13 中夺冠的 AlphaFold 虽然在预测蛋白质结构中取得一定程度的效果，但在Alphafold的github中有这么一段话。![img](https://i.loli.net/2021/06/26/4t7rgSlKLVAP8qv.jpg)

Alphafold的出现，在一定程度上缓解了人工的压力，但是却对硬件有极高的要求，适用性不强，所以在一定程度上 AlphaFold 的胜利也可以说是 DeepMind 硬件的胜利，并没有从根本上找到「机器学习」解。

 

## **3.2.** ***\*Alphafold2\****

但在2020年CASP14中，不少工业界的队伍也参与进来，除了Deepmind 的alphafold2，比如微软的BrainFold，微软亚洲研究院的两款算法，TOWER and FoldX 和 NOVA，腾讯的数款算法tFold，但是他们的结果在Alphafold2面前只能用两个字形容：血虐。

 

![img](https://i.loli.net/2021/06/26/QBZ4GjqvAFLtwrW.png) 

而在新的Alphafold2模型中，抛弃了传统的CNN，RNN，选择利用NLP中比较火的transformer和注意力机制。

![img](https://i.loli.net/2021/06/26/lCjknUfu7gAJdmo.jpg) 

下面是Deepmind博客提供的结果，绿色表示实验结果，蓝色表示预测结果，我们可以看到二者几乎相差无几。

![img](https://i.loli.net/2021/06/26/Rr4EAaWLFZNpygu.png) 

正如[Nature的一篇评论](https://link.zhihu.com/?target=https://www.nature.com/articles/d41586-020-03348-4)，人们可以花更多的时间思考，花更少的时间拿移液枪了，如果说Alphafold没有从根本问题上找到机器学习的解，那么Alphafold2已经找到了机器学习的近似解，Alphafold2的出现对于生物学具有里程碑式的意义。

 

Reference：

 [AlphaFold2。 - 知乎](https://zhuanlan.zhihu.com/p/315497173)

 [An Introduction to AlphaFold and Protein Modeling • brett koonce](https://brettkoonce.com/talks/an-introduction-to-alphafold-and-protein-modeling/)

 [How AlphaFold Works in Predicting 3D Protein Shapes [3\] | Bio-Med | PI IP LAW](https://piip.co.kr/en/blog/how-AlphaFold-works)

 [DeepMind's AlphaFold 2 Explained! AI Breakthrough in Protein Folding! What we know (& what we don't) - YouTube](https://www.youtube.com/watch?v=B9PL__gVxLI)

 
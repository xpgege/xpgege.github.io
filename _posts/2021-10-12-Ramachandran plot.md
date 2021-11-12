---
layout:     post
title:      Ramachandran Plot
subtitle:   Ramachandran Plot
date:       2021-10-12
author:     xpgege
header-img: img/post-web.jpg
catalog: true
tags:
	- 知识记录
---

Ramachandran 图是肽中包含的残基（氨基酸）的扭转角 - phi (φ) 和 psi (ψ) - 的图。φ 是 N(i-1),C(i),Ca(i),N(i) 扭转角，ψ 是 C(i),Ca(i),N(i),C (i+1) 扭转角。

<img src="https://pic3.zhimg.com/80/v2-7fad8886893037900e9509e52eb1af3d_720w.jpg?source=1940ef5c" alt="img" style="zoom:80%;" />



由 G. N. Ramachandran 等人于 1963 年开发。通过在 x 轴上绘制 φ 值，在 y 轴上绘制 ψ 值，

![Fig. 23](https://ars.els-cdn.com/content/image/3-s2.0-B9780128096338204846-f20484-23-9780128114148.jpg)

以这种方式绘制扭转角以图形方式显示哪些角度组合是可能的。肽中每个残基的扭转角通过相对于两个相邻的平面肽键定位其平面肽键来定义其与其两个相邻残基的连接的几何形状，从而扭转角决定了残基和肽的构象。由于空间位阻，许多角度组合以及因此残基的构象是不可能的。通过绘制Ramachandran图，蛋白质结构科学家可以确定允许哪些扭转角，并可以深入了解肽的结构。

## 1.二级结构绘图

肽的二级结构是具有有序重复结构的肽段，重复结构是由于残基的重复构象，最终是 φ 和 ψ 的重复值。

不同的二级结构可以通过它们的 φ 和 ψ 值范围来区分，不同二级结构的值映射到Ramachandran图的不同区域。

### 1.1 α-helix

​	α-helix的Ramachandran图具有聚集在 φ= -57°和 ψ= -47° 值附近的点，它们是 α-螺旋的平均值。将另外两个螺旋段的点（不同颜色）加入，而来自所有三个螺旋段的数据都出现在一个大集群中，说明无法通过它们的 φ 和 ψ 值的差异来区分螺旋段

![image-20211012150317030](https://i.loli.net/2021/10/12/Vbj6p2Eqzk9Wext.png)

### 1.2 β-sheets

​	β-sheets的Ramachandran图具有聚集在 φ= -130°和 ψ= +140° 值附近的点，



## 2.空间位阻绘图

大多数 φ 和 ψ 的组合在空间上是被禁止的，如三肽 Glu-Ser-Ala 所示，大多数 φ 和 ψ 的组合在空间上是被禁止的。 

天然肽图中，数据点将在未发生空间位阻的几个区域中形成簇。

核心区域（图中蓝色）包含最有利的 φ 和 ψ 组合，并且包含最多的点。一些蛋白质的图在右上象限包含一个小的第三个核心区域。

允许的区域（图中绿色）可以位于核心区域周围，也可以与核心区域无关，但它们包含的数据点少于核心区域。

剩余的区域被认为是不允许的。

![Image:Ramaplot.png](https://proteopedia.org/wiki/images/5/57/Ramaplot.png)

### 2.1 Glycine, Proline, etc.

现代Ramachandran 标准对具有不同局部空间位阻特性的氨基酸子集使用单独的函数

由于 Gly 只有一个氢作为侧链，因此当 φ 和 ψ 旋转一系列值时，空间位阻不太可能发生。出于这个原因，Gly 会经常在一般情况下的 Ramachandran 图的不允许区域中绘图。上图中禁止区域中的几乎所有数据点都是 Gly 点。

脯氨酸的侧链与前面的主链 N 共价连接，比一般情况下的残基受到更严格的限制。 Pro（称为“prePro”）之前的残基与脯氨酸环的空间相互作用有一些限制。

 Ile 和 Val 的分支 β 碳也使它们具有独特的禁止拉马钱德兰区域的形状。

其他 16 种氨基酸类型对非常有利区域的偏好各不相同，但它们将允许区域与异常区域分开的外部轮廓都非常一致，因此它们都被归为“一般情况”分布

![img](https://www.bonvinlab.org/education/molmod_online/ramachandran.png)



### 2.2 Functionally relevant residues

​	功能相关的残基比其他残基更有可能具有绘制到拉马钱德兰图的允许但不受欢迎的区域的扭转角。这些功能相关残基的特定几何形状虽然在能量上有些不利，但可能对蛋白质的催化或其他功能很重要。这种构象需要通过蛋白质使用 H 键、空间堆积或其他方式来稳定，并且对于高度暴露于溶剂的残基应该很少发生。

最后说明，Ramachandran Plot只是氨基酸残基扭转角的一种可视化表示，没有考虑能量问题，所以拉氏图也只是一个同源建模最基本的检测分析，作为一个参考。

## 3.绘制Ramachandran 图

**PROCHECK**



**PROCHECK**是最常用的蛋白质几何形状评估工具之一，主要用于评价蛋白质结构的立体化学参数，它不仅适用于实验过程中已经结晶的蛋白，对同源建模所得模型也同样适用。它忽略蛋白质系统的能量，主要研究氨基酸残基碳骨架的二面角在蛋白质或多肽链中的分布情况，并生成Ramachandran统计图。该图被不同颜色的曲线划分为 4 部分，分别是最合理区（the most favored regions），其它合理区（the additional allowed regions），一般合理区（the generous allowed regions），不合理区（the disallowed regions）。从Ramachandran 图中我们可以得到各个氨基酸残基之间的二面角φ 、ψ，当不少于90%的二面角分布在最佳合理区时，则蛋白质的结构是合理的。



[http://www.csb.yale.edu/poststructure/procheck/procheck.html](https://link.zhihu.com/?target=http%3A//www.csb.yale.edu/poststructure/procheck/procheck.html)



![img](https://pic3.zhimg.com/80/v2-d3285b58ee4e2299d65ed65500935a66_720w.jpg)

**VERIFY_3D**



**VERIFY_3D**对同源建模所得模型的质量进行了可视化分析，主要是用来判断模型与其本身氨基酸序列之间的兼容性，对氨基酸残基数大于100的蛋白质，VERIFY_3D的评估结果更准确。一般当不低于80%的氨基酸残基得分大 0.2时，即可认为目的蛋白建模所得模型属于高质量的模型结构。





[https://servicesn.mbi.ucla.edu/Verify3d/](https://link.zhihu.com/?target=https%3A//servicesn.mbi.ucla.edu/Verify3d/)



![img](https://pic3.zhimg.com/80/v2-d3285b58ee4e2299d65ed65500935a66_720w.jpg)

**PROSA**





**PROSA**的评价结果反映的是目的蛋白三维结构氨基酸残基之间的相互作用能，由公式计算得到：



![img](https://pic2.zhimg.com/80/v2-0d75a792a38fff7c63ff4f98be271f31_720w.jpg)



其中 ES,C表示序列S间的残基在C构象空间中的相互作用能，而 ÊS,C指序列S间的残基在C构象空间的相互作用能的平均值，下面的σS为相对应的标准偏差。将含有目的蛋 白 空 间 结 构 的PDB格 式 的 文 件 上 传 到 在 线 服 务 器（[https://www.came.sbg.ac.at/prosa.php](https://link.zhihu.com/?target=https%3A//www.came.sbg.ac.at/prosa.php)）便可得到PROSA的评价结果，计算所得数值ZS,C即为 Z-score，该值通常为负数。在 PROSA 结果图中还包括所有已经实验确定的现存于 RCSB PDB数据库中与目的蛋白大小相近的蛋白质的Z-score值，当目标蛋白的Z-score 值分布在这些已知蛋白的Z-score值绘制的图形范围内时，则证明目标蛋白结构的能量是合理的，反之，则认为结构是不合理的。



[https://prosa.services.came.sbg.ac.at/prosa.php](https://link.zhihu.com/?target=https%3A//prosa.services.came.sbg.ac.at/prosa.php)





reference：

https://proteopedia.org/wiki/index.php/Ramachandran_Plots

https://zhuanlan.zhihu.com/p/399348442




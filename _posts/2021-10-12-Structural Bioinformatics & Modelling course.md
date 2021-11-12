---
layout:     post
title:      Structural Bioinformatics & Modelling
subtitle:   Structural Bioinformatics & Modelling
date:       2021-10-12
author:     xpgege
header-img: img/post-web.jpg
catalog: true
tags:
	- 知识记录

---

# 1.Homology modeling 

原文教程通过swiss model进行同源建模分析，[Homology modelling](https://www.bonvinlab.org/education/molmod_online/modelling)。

随着alphafold的出现，其高准确率和精度足以替换绝大多数同源建模方法

- 可以进入注册登陆谷歌[colab](https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb)，直接网页端进行预测（需要翻墙），但这是一个简化版本，但耗时更低。
- 若需要高精度预测，请根据其[github](https://github.com/deepmind/alphafold)源码进行配置,也可参考我之前的安装采坑。



课程研究的蛋白质为MDM2 mouse protein，它的N端区域会结合p53反式激活domain，序列为：

```
>sp|P23804|MDM2_MOUSE E3 ubiquitin-protein ligase Mdm2 OS=Mus musculus OX=10090 GN=Mdm2 PE=1 SV=3
MCNTNMSVSTEGAASTSQIPASEQETLVRPKPLLLKLLKSVGAQNDTYTMKEIIFYIGQY
IMTKRLYDEKQQHIVYCSNDLLGDVFGVPSFSVKEHRKIYAMIYRNLVAVSQQDSGTSLS
ESRRQPEGGSDLKDPLQAPPEEKPSSSDLISRLSTSSRRRSISETEENTDELPGERHRKR
RRSLSFDPSLGLCELREMCSGGSSSSSSSSSESTETPSHQDLDDGVSEHSGDCLDQDSVS
DQFSVEFEVESLDSEDYSLSDEGHELSDEDDEVYRVTVYQTGESDTDSFEGDPEISLADY
WKCTSCNEMNPPLPSHCKRCWTLRENWLPDDKGKDKVEISEKAKLENSAQAEEGLDVPDG
KKLTENDAKEPCAEEDSEEKAEQTPLSQESDDYSQPSTSSSIVYSSQESVKELKEETQDK
DESVESSFSLNAIEPCVICQGRPKNGCIVHGKTGHLMSCFTCAKKLKKRNKPCPVCRQPI
QMIVLTYFN
```

具体蛋白质信息参考uniprot：https://www.uniprot.org/uniprot/P23804



uniprot中存在两种结构，一种是Alphafold预测出的结构（https://alphafold.ebi.ac.uk/entry/P23804），储存在[AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)中

![image-20211012163820051](https://i.loli.net/2021/10/12/QDw3MSCbihNIYRG.png)

但个人感觉这个结构看起来就很不靠谱



另一种是同源建模产生的一个二聚体（https://swissmodel.expasy.org/repository/uniprot/P23804?csm=4ABF489A82038DF4）

![image-20211012164213981](https://i.loli.net/2021/10/12/emqRsAa1tfn85Fi.png)

看起来结构尚可，通过这个链接下载pdb结构，https://swissmodel.expasy.org/repository/6160d883d3f30106a7792464.pdb

# 2.Molecular Dynamics Simulation of the p53 N-terminal peptide

### Selecting an initial structure

第一步显然是选择起始结构。本教程的目的是模拟 p53 反式激活域 N 端序列的肽段。下面以 FASTA 格式给出了该肽的序列

```text
>P53_MOUSE
SQETFSGLWKLLPPE
```

肽通常是非常灵活的分子，具有短暂的二级结构元素。有些甚至可以采用不同的结构，具体取决于它们与哪种蛋白质伴侣相互作用，如果在溶液中游离，则保持无序状态。因此，使用高级方法（例如对该肽段进行同源性建模）的努力很可能是没有根据的。相反，生成三种理想构象（螺旋、片状和聚脯氨酸-2）的肽结构是可能的，而且是合理的，这些构象已被证明代表了 PDB 中沉积的大多数肽。生成这些结构是一个简单的操作骨架二面角的问题。 Pymol 有一个实用程序脚本来执行此操作，由 Robert Campbell 编写，如有必要，可在此处获取http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/

```
build_seq peptide_helix, SQETFSGLWKLLPPE, ss=helix
save p53_helix.pdb, peptide_helix
```

### Preparing the initial structure

对于从 RCSB PDB 下载的结构，重要的是确保没有丢失原子，并检查是否存在非标准氨基酸和其他小配体。

力场通常包含天然氨基酸和核苷酸、一些翻译后修饰、水和离子的参数。药物和辅助因子等外来分子通常必须手动参数化，这本身就是一门科学。始终判断这些外来物种的存在是否必要。在某些情况下，可以安全地忽略配体并将其从结构中移除。至于缺失的残基和原子，除了氢，在开始模拟之前重建它们是绝对必要的。 MODELLER 是用于此目的的出色程序。此外，一些晶体以足够好的分辨率衍射以区分密度网格中的水分子。除了这些水域是研究对象的非常特殊的情况外，最好的策略是将它们从结构中完全移除。幸运的是，大多数这些“有问题”的分子在 PDB 文件中显示为杂原子 (HETATM)，因此可以使用简单的 sed 命令轻松删除：

```
sed -e ‘/^HETATM/d’ 1XYZ.pdb > 1XYZ_clean.pdb
```

由于 p53 肽的初始结构是使用 Pymol 和理想几何形状生成的，因此无需进行此类检查

### Structure Conversion and Topology Generation

一个分子不仅由其原子的三维坐标定义，而且还由这些原子如何连接以及它们如何相互作用的描述来定义。在上一步中生成或下载的 PDB 文件仅包含前者。

系统在原子类型、电荷、键等方面的描述包含在拓扑中，它特定于模拟中使用的力场。力场的选择不能掉以轻心。对于生物分子系统，几乎没有主要的力场——例如CHARMM, AMBER, GROMOS, OPLS – 已被参数化以重现生物分子（即蛋白质）的特性。自分子动力学模拟的第一天以来，这一直是并且继续是一个活跃的研究领域。 Pubmed 中有几篇文献综述可以评估每个力场及其多个版本的质量和适当性。有些以其artifacts,而闻名，例如对 alpha 螺旋构象的偏向倾向。

在这里，在本教程中，我们使用 AMBER99SB-ILDN 力场，它广泛用于采样和折叠模拟，并已被证明可以很好地再现实验数据（来源）。这个选择背后的另一个更实际的原因是 GROMACS 中这个力场的可用性

由于模拟发生在溶剂化环境中，即一盒水分子，我们还必须选择合适的溶剂模型。该模型只是对包含水分子特性的力场的补充，参数化以重现特定特性，例如密度、冻结和汽化温度。因此，特定的水模型往往与特定的力场相关联。由于在计算上再现水的特性很困难——是的，即使是这样一个简单的分子！ - 一些模型代表具有超过 3 个原子的水，使用额外的伪粒子来改善其静电分布等特性。建议与 AMBER 力场一起使用的水模型，我们将在本次模拟中使用的是 TIP3P 模型（用于 3 点的可转移相互作用势），它实际上是由 OPLS 力场的作者开发的。真的可以转让！

```
gmx pdb2gmx -f peptide.pdb -o peptide.gro -p peptide.top -ignh -ter
```

GROMACS 程序 pdb2gmx 采用初始结构并返回拓扑文件 (peptide.top) 和遵守力场原子命名约定的新结构 (peptide.gro)。为了转换结构和构建拓扑，pdb2gmx 将分子分成几个块，例如氨基酸，并使用此类构建块的力场特定库进行必要的转换。通常，与库的匹配是通过 PDB 文件中每个 ATOM/HETATM 行上的残基/原子名称来完成的。如果残基（或原子）未被识别，程序将停止并返回错误。

不同的力场定义不同的原子类型和/或为相同的原子类型赋予不同的名称。虽然大多数重原子（即非氢）在大多数力场中具有相同的命名，但氢却没有。因此，标志 -ignh 表示 GROMACS 在读取结构时应忽略这些原子，并使用力场中定义的理想几何参数（重新）生成它们的坐标。

此外，该程序允许用户通过 -ter 标志定义分子末端的状态。 Termini 可以带电（例如 NH3+ 和 COO-）、不带电（例如 NH2 和 COOH），或者被额外的化学基团封端（例如 N 端乙酰基和 C 端酰胺）。这非常重要，因为让末端带电（默认）会导致人工电荷-电荷相互作用，特别是在小分子中。如果肽是更大结构的一部分，那么将末端加帽以中和它们的电荷是有意义的，因为这在现实中会发生。通读 pdb2gmx 的输出并检查程序对组氨酸质子化状态的选择以及由此产生的肽电荷。

```
gmx editconf -f peptide.gro -o peptide.pdb

#What are the differences between both files. What did GROMACS add/remove to the structure?
```

新生成的拓扑文件也值得关注。它包含所有残基及其相应原子的列表，详细说明了原子类型、质量和电荷。此外，它还包含分子中所有键、角度和二面角的列表。请注意，拓扑文件不包含有关其化学的任何信息。此信息存储在拓扑文件最顶部定义的内部参数库中。

```text
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
; residue   1 GLN rtp NGLN q +1.0
     1         N3      1    GLN      N      1     0.1493      14.01   ; qtot 0.1493
     2          H      1    GLN     H1      2     0.1996      1.008   ; qtot 0.3489
     3          H      1    GLN     H2      3     0.1996      1.008   ; qtot 0.5485
     4          H      1    GLN     H3      4     0.1996      1.008   ; qtot 0.7481
     5         CT      1    GLN     CA      5     0.0536      12.01   ; qtot 0.8017
     6         HP      1    GLN     HA      6     0.1015      1.008   ; qtot 0.9032
     7         CT      1    GLN     CB      7     0.0651      12.01   ; qtot 0.9683
     8         HC      1    GLN    HB1      8      0.005      1.008   ; qtot 0.9733
     9         HC      1    GLN    HB2      9      0.005      1.008   ; qtot 0.9783
    10         CT      1    GLN     CG     10    -0.0903      12.01   ; qtot 0.888
    11         HC      1    GLN    HG1     11     0.0331      1.008   ; qtot 0.9211
    12         HC      1    GLN    HG2     12     0.0331      1.008   ; qtot 0.9542
    13          C      1    GLN     CD     13     0.7354      12.01   ; qtot 1.69
    14          O      1    GLN    OE1     14    -0.6133         16   ; qtot 1.076
    15          N      1    GLN    NE2     15    -1.0031      14.01   ; qtot 0.0732
    16          H      1    GLN   HE21     16     0.4429      1.008   ; qtot 0.5161
    17          H      1    GLN   HE22     17     0.4429      1.008   ; qtot 0.959
    18          C      1    GLN      C     18     0.6123      12.01   ; qtot 1.571
    19          O      1    GLN      O     19    -0.5713         16   ; qtot 1
```

### Periodic Boundary Conditions

这种转换后的结构包括几个原子，即氢，它们仅根据理想的几何参数添加。如果使用 Pymol 生成，它也具有理想的骨架几何形状。如果它是从 RCSB PDB 下载的，则结构也可能包含某些化学方面（键长、角度、原子间距离），这些方面被力场认为是不理想的。

事实上，仅仅改变力场也会导致理想的定义发生变化。准备系统的第一步是尽可能地消除这些“缺陷”，这通常是通过系统的能量最小化来实现的。这种优化方法本质上迫使一组原子尽可能地遵守力场的定义。系统中原子的数量越多，就越难让所有原子完美地符合所有定义。例如，将两个原子靠得更近以减少违反范德华力强加的定义的应变可能会导致静电项的应变增加。

在最小化系统之前，必须选择模拟设置的总体布局。换句话说，肽必须放置在某个地方才能实现这种最小化。大多数现代蛋白质和肽模拟都定义了周期性边界条件 (PBC)，它设置了一个可以无限堆叠的单个晶胞。因此，定义了一个无限的周期系统，避免了分子可以真正碰到的硬边界（壁）的问题。当蛋白质穿过左侧的壁时，其右侧的周期性图像进入当前晶胞，在每个晶胞中保持恒定数量的原子。将 PBC 合理化的一种更简单的方法是将它们与旧诺基亚手机中的蛇游戏进行比较。当蛇的头部越过屏幕边界时，它会重新出现在截然相反的边缘。

晶胞形状的选择也很重要，因为这将定义模拟分子的体积。分子动力学模拟在计算上要求很高。系统中的分子越多，每一步需要计算的力就越多。因此，虽然立方体可以无限地完美堆叠，但从体积的角度来看，它并不是最有效的形状（请记住，模拟通常发生在溶剂化环境中！）。将形状近似为球体是理想的，但球体不能堆叠。因此，只有少数通用形状支持周期性边界条件的设置。其中之一是菱形十二面体，它对应于球体的最佳堆积，因此是自由旋转分子（如肽或蛋白质）的最佳选择。

设置 PBC 时要记住的另一件事是单元格的大小。继续用蛇类比，让蛇的头看到自己的尾巴是不合适的。换句话说，单元格必须足够大以允许分子跨越边界，并且仍然与下一个图像保持足够的距离，以便在它们之间不进行力计算。在 GROMACS 中，此设置定义为从分子到晶胞壁的距离。这个距离也不应该任意大，否则盒子太大，模拟计算效率低下。以用于计算力场中非键相互作用（长距离）的截止值作为经验法则。到墙壁的距离必须大于这个值。

```
#使用肽和晶胞壁之间 1.4 nm 的最小距离设置周期性边界条件
gmx editconf -f peptide.gro -o peptide-PBC.gro -bt dodecahedron -c -d 1.4
```

与 pdb2gmx 一样，GROMACS 程序 editconf 生成一个相当大的输出，其中包含例如它刚刚创建的单位单元的体积和尺寸。维度使用三斜矩阵表示，其中前三个数字指定对角线元素（xx、yy、zz），后六个数字指定非对角线元素（xy、xz、yx、yz、zx、zy）。

### Energy minimization of the structure in vacuum

定义了可以进行模拟的物理空间后，现在可以最小化分子的能量。 GROMACS 对任何涉及分子和力场的计算使用两步法。首先，用户必须将结构和拓扑数据以及仿真参数组合在一个控制文件中。该文件包含有关系统的所有信息，并确保模拟的再现性，前提是机器上具有相同的力场。

拥有这样一个独立文件的另一个优点是准备工作可以在一台机器上进行，而计算则在另一台机器上运行。同样，模拟在计算上要求很高。虽然该系统可以在笔记本电脑上轻松准备，借助 Pymol、支持 GUI 的文本编辑器以及具有屏幕的所有其他优势，但计算通常在具有数百个处理核心的专用集群上运行，这些核心仅提供一个命令——线路接口接入。这在运行生产模拟时很重要。准备系统的中间计算非常小，可以在笔记本电脑上运行。

模拟参数包含在单独的文件中，通常带有 .mdp 扩展名。为简单起见，我们在我们的 GitHub 存储库中提供这些文件），如果您正在使用它，也已经在我们的虚拟映像中（请参阅 $MOLMOD_DATA/mdp/）。例如，这些参数指定用于计算非键相互作用的截止值、用于计算每个原子的邻居的算法、周期性边界条件的类型（例如三维、二维）以及用于计算的算法非键相互作用。他们还指定了模拟的类型，例如能量最小化或分子动力学，以及它的长度和时间步长（如果合适）。最后，他们还描述了 GROMACS 将坐标和能量值写入磁盘的频率。根据模拟的目的，可以增加此写入频率以具有更高的时间分辨率，但会牺牲一些计算效率（写入需要时间）。 MDP 文件支持数百个参数设置，所有这些在 GROMACS 手册中都有详细说明 [GROMACS manual](https://manual.gromacs.org/2019-current/user-guide/mdp-options.html).

```
gmx grompp -v -f $MOLMOD_DATA/mdp/01_em_vac_PME.mdp -c peptide-PBC.gro -p peptide.top -o peptide-EM-vacuum.tpr -maxwarn 1
gmx mdrun -v -deffnm peptide-EM-vacuum
```

尽管 GROMACS 由几个实用程序组成，但它的核心是 mdrun 程序。正是这段代码运行了所有的模拟。 -deffnm 标志是一个非常方便的选项，它为所有文件选项（包括输入和输出）设置默认文件名，避免多个单独的定义。 -v  verbose，在这种情况下，打印系统的势能和最小化每一步的最大力。

```
Steepest Descents:
   Tolerance (Fmax)   =  1.00000e+01
   Number of steps    =         5000
Step=    0, Dmax= 1.0e-02 nm, Epot=  4.80138e+03 Fmax= 1.83867e+04, atom= 57
Step=    1, Dmax= 1.0e-02 nm, Epot=  2.69755e+03 Fmax= 6.83824e+03, atom= 56
Step=    2, Dmax= 1.2e-02 nm, Epot=  1.46780e+03 Fmax= 4.07714e+03, atom= 57
Step=    3, Dmax= 1.4e-02 nm, Epot=  1.69036e+03 Fmax= 1.47891e+04, atom= 56
Step=    4, Dmax= 7.2e-03 nm, Epot=  1.19448e+03 Fmax= 6.00160e+03, atom= 56
Step=    5, Dmax= 8.6e-03 nm, Epot=  1.03838e+03 Fmax= 4.85784e+03, atom= 56
Step=    6, Dmax= 1.0e-02 nm, Epot=  1.11613e+03 Fmax= 1.01399e+04, atom= 56
Step=    7, Dmax= 5.2e-03 nm, Epot=  8.98891e+02 Fmax= 2.73856e+03, atom= 56
Step=    8, Dmax= 6.2e-03 nm, Epot=  8.39895e+02 Fmax= 4.76931e+03, atom= 56
Step=    9, Dmax= 7.5e-03 nm, Epot=  8.05094e+02 Fmax= 5.81049e+03, atom= 56
Step=   10, Dmax= 9.0e-03 nm, Epot=  7.77891e+02 Fmax= 5.97918e+03, atom= 56
```

此最小化中使用的最速下降算法计算系统每一步的能量梯度，并提取将系统推向能量最小值的力。因此，势能必须减少。对于分子动力学和其他最小化算法，情况并非如此。当满足以下两个条件之一时，最小化结束：要么最大力小于提供的阈值 (10 kJ.mol-1)，并且最小化收敛，要么算法达到参数文件中定义的最大步数（ 5000）。理想情况下，最小化应该运行直到收敛，但除了非常特定的场景，如正常模式分析，这不是一个严格的要求。

### Solvating the simulation box

下一步是将溶剂添加到模拟框。第一次蛋白质分子动力学模拟是在真空中完成的，但研究人员很快意识到这是一个主要限制。水分子与蛋白质相互作用，介导残基之间的相互作用。此外，作为溶剂的水对长程相互作用（例如静电）具有屏蔽作用。在真空中，即使在很远的距离内，也没有什么可以阻止两个相反电荷的原子相互感知，只要它们在用于模拟的截止范围内即可。随着水的加入，这种相互作用显着减弱。在选择盒子的大小时，水介导的相互作用的影响也很重要。溶质肽的存在会诱导其附近水分子的特定排序。这可能会产生涟漪效应，传播溶质的影响并导致远远超出理论非键合截止值的伪影。

在运行 pdb2gmx 时，您应该已经选择了合适的水模型 - TIP3P。求解不需要拓扑文件。本质上，这个操作只是将预先计算好的水分子块放入盒子中，并去除那些与蛋白质原子重叠的水分子。不涉及化学。但是，必须更新拓扑以反映溶剂的添加。

```
gmx solvate -cp peptide-EM-vacuum.gro -cs spc216.gro -o peptide-water.gro -p peptide.top
```

GROMACS 在更新之前备份先前的拓扑文件。通常，GROMACS 从不覆盖文件，而是复制前一个文件并用 # 符号重命名它。在新拓扑文件的末尾，有一个额外的条目列出了现在结构中的水分子数量。它还添加了一个加载水模型参数的定义。

```
GROMACS 在更新之前备份先前的拓扑文件。通常，GROMACS 从不覆盖文件，而是复制前一个文件并用 # 符号重命名它。在新拓扑文件的末尾，有一个额外的条目列出了现在结构中的水分子数量。它还添加了一个加载水模型参数的定义。
```

### Addition of ions: counter charge and physiological concentration

除了水之外，细胞环境还包含许多离子，这些离子可以保持系统的某种化学中性。将其中一些添加到模拟框也增加了模拟的真实性。 GROMACS 程序 genion 执行此任务，但需要 .tpr 文件作为输入。离子的添加是通过替换模拟框中已经存在的某些原子来完成的。由于不太需要去除肽的原子，请注意您在运行 genion 时选择的组。 -neutral 标志表示允许过量的一种离子种类来中和系统的电荷（如果有的话）。

```
gmx grompp -v -f $MOLMOD_DATA/mdp/02_em_sol_PME.mdp -c peptide-water.gro -p peptide.top -o peptide-water.tpr -maxwarn 1
gmx genion -s peptide-water.tpr -o peptide-solvated.gro -conc 0.15 -neutral -pname NA+ -nname CL-
```

### Energy minimization of the solvated system

添加离子是为模拟设置系统（化学）的最后一步。从这里开始，所需要做的就是以受控方式放松系统。添加溶剂和离子可能会导致一些不利的相互作用，例如重叠的原子和放置得太近的相等电荷。

```
mx grompp -v -f $MOLMOD_DATA/mdp/02_em_sol_PME.mdp -c peptide-solvated.gro -p peptide.top -o peptide-EM-solvated.tpr
gmx mdrun -v -deffnm peptide-EM-solvated
```

### Restrained MD – relaxation of solvent and hydrogen atoms

尽管耗散了系统中的大部分应变，但能量最小化不考虑温度，因此不考虑速度和动能。当第一次运行分子动力学时，算法会为原子分配速度，这再次给系统带来压力，并可能导致模拟变得不稳定。为了避免可能的不稳定性，这里描述的准备设置包括分子动力学的几个阶段，这些阶段逐步消除对系统的限制，因此，让它慢慢适应生产模拟运行的条件。

此模拟的 .mdp 文件与用于最小化运行的文件大不相同。首先，积分器现在是 md，它指示 mdrun 实际运行分子动力学。然后，有几个与此算法特别相关的新选项：dt、t_coupl、ref_t 和 gen_vel。在文件的顶部，有一个预处理选项，它定义了一个特定的标志 -DPOSRES。在拓扑文件中，有一个特定的语句只有在设置了这个标志时才被激活，它与 pdb2gmx 创建的文件有关——posre.itp。该文件包含系统某些原子的位置限制，防止它们在模拟过程中自由移动。

```
cp $MOLMOD_DATA/mdp/03_nvt_pr1000_PME.mdp ~/
gmx grompp -v -f ~/03_nvt_pr1000_PME.mdp -c peptide-EM-solvated.gro -r peptide-EM-solvated.gro -p peptide.top -o peptide-NVT-PR1000.tpr
gmx mdrun -v -deffnm peptide-NVT-PR1000
```

在这个系统中包含速度导致粒子和系统获得动能。此信息以扩展名为 .edr 的二进制文件格式存储，可以使用 GROMACS 公用事业能源读取。该实用程序将能量文件中的信息提取到表格文件中，然后可以将其转换为图表。通过依次键入它们的数字然后按 Enter 来选择感兴趣的术语。要退出，请键入 0 并按 Enter。使用 xvg_plot.py 实用程序绘制生成的 .xvg 文件，传递 -i 标志以打开交互式会话。如果要更改绘图的颜色，请使用 -h 标志运行脚本并参考此页面以获取可用的颜色图。https://github.com/JoaoRodrigues/gmx-tools

```
gmx energy -f peptide-NVT-PR1000.edr -o thermodynamics-NVT-PR1000.xvg
$MOLMOD_BIN/xvg_plot.py -i thermodynamics-NVT-PR1000.xvg
```

### Coupling the barostat – simulating in NPT conditions

平衡通常分两个阶段进行：首先，系统在经典系综 (NVT) 下进行模拟，其中分子数量、体积和温度保持恒定。目标是让系统达到并稳定在所需的温度。第二步是将恒压器耦合到模拟中并保持恒定压力，这更类似于实验条件。通过调节粒子的速度来控制温度，而通过改变模拟箱的体积（PV=NRTPV=NRT）来保持压力恒定。

```
gmx grompp -v -f $MOLMOD_DATA/mdp/04_npt_pr_PME.mdp -c peptide-NVT-PR1000.gro -r peptide-NVT-PR1000.gro -p peptide.top -o peptide-NPT-PR1000.tpr
gmx mdrun -v -deffnm peptide-NPT-PR1000
gmx energy -f peptide-NPT-PR1000.edr -o thermodynamics-NPT-PR1000.xvg
$MOLMOD_BIN/xvg_plot.py -i thermodynamics-NPT-PR1000.xvg
```

### Releasing the position restraints

到现在为止，系统有时间调整注入速度以及引入温度和压力。然而，肽的重原子仍被限制在它们的初始位置。模拟设置的下一步和最后一步逐步释放这些限制，直到系统完全不受限制并在所需的温度和压力下完全平衡，从而为生产模拟做好准备

约束的强度在由 pdb2gmx 创建的 posre.itp 文件中定义。力常数的值定义了原子受约束的严格程度。因此，解除限制就像修改文件上的数字一样简单。

```
[ position_restraints ]
; atom  type	  fx	  fy	  fz
     1     1    1000  1000  1000
     4     1    1000  1000  1000
     6     1    1000  1000  1000
     9     1    1000  1000  1000
    12     1    1000  1000  1000
    13     1    1000  1000  1000
    14     1    1000  1000  1000
    17     1    1000  1000  1000
```

```
#Decrease the strength of the force constant of the position restraints and re-run the system under NPT.

cp posre.itp posrest.itp.1000 # Make a backup of the original file
sed -i -e 's/1000  1000  1000/ 100   100   100/g' posre.itp
gmx grompp -v -f $MOLMOD_DATA/mdp/04_npt_pr_PME.mdp -c peptide-NPT-PR1000.gro -r peptide-NPT-PR1000.gro -p peptide.top -o peptide-NPT-PR100.tpr
gmx mdrun -v -deffnm peptide-NPT-PR100

cp posre.itp posrest.itp.100
sed -i -e 's/100   100   100/ 10    10    10/g' posre.itp
gmx grompp -v -f $MOLMOD_DATA/mdp/04_npt_pr_PME.mdp -c peptide-NPT-PR100.gro -r peptide-NPT-PR100.gro -p peptide.top -o peptide-NPT-PR10.tpr
gmx mdrun -v -deffnm peptide-NPT-PR10
```

最后的平衡步骤是完全消除位置限制。这是通过删除 .mdp 文件开头的 -DPOSRES 定义来完成的，同时保留所有其他参数。为简单起见，我们提供了另一个没有此定义的 .mdp 文件。

```
gmx grompp -v -f $MOLMOD_DATA/mdp/05_npt_NOpr_PME.mdp -c peptide-NPT-PR10.gro -p peptide.top -o peptide-NPT-noPR.tpr
gmx mdrun -v -deffnm peptide-NPT-noPR
```

### Production Simulation

尽管做出了所有这些努力，但系统不太可能已经处于平衡状态。模拟的前几纳秒，取决于系统，实际上是一个平衡期，在对感兴趣的属性进行任何分析时应该放弃。要为生产设置模拟，只需生成一个包含所需参数的新 .tpr 文件，即定义模拟长度的步骤数。在这个阶段，有很多问题需要解决，这些问题对计算性能有不同程度的影响：

所研究的过程发生在什么时间尺度？模拟应该运行多长时间？

回答研究问题所需的时间分辨率是多少？

是否需要经常存储速度和能量？

应该多久将模拟信息写入日志文件？

模拟将运行 50 纳秒，这足以让我们深入了解这种小肽的构象动力学。请记住，对整个景观进行全面和详尽采样的适当模拟应该持续更长时间，并且可能会使用更先进的分子动力学协议，例如副本交换。在这种情况下，由于预计有几个学生使用不同的随机种子并从不同的初始构象开始研究相同的肽，因此我们假设 50 纳秒的单个模拟足以提供信息。

```
cp $MOLMOD_DATA/mdp/06_md_PME.mdp ~/06_md_PME.mdp
gmx grompp -v -f ~/06_md_PME.mdp -c peptide-NPT-noPR.gro -p peptide.top -o p53_helix_CAH.tpr
```

如果您希望检查 .tpr 文件的内容，请使用 GROMACS 的转储实用程序，顾名思义，它会将文件的全部内容输出到屏幕上。将命令的输出通过管道传输到文本处理器，例如 less 或 more（Linux 笑话）以对输出进行分页。按 q 退出程序。

gmx dump -s p53_helix_CAH.tpr | more

### Analysis of the Molecular Dynamics Simulation

生产运行只是分子动力学模拟背后真正工作的开始。模拟分析可以分为几个部分，并且根据模拟的目标和所提出的研究问题而有很大差异。通常，分析的第一部分是从整体上评估仿真的质量和稳定性。如果这些表明模拟存在任何问题，即周期性图像相互作用、不稳定的温度或压力，或溶质的不受控制的动力学（即蛋白质的意外展开），则可能必须重复模拟。如果模拟是稳定的，则分析会继续提取可能有助于回答研究问题的数据。

生产模拟产生许多文件，每个文件包含不同的信息。根据提供给 mdrun 的选项，名称可能会有所不同。但是，扩展名保持不变。对于大多数分析，唯一的要求是压缩轨迹 (.xtc) 和能量 (.edr) 文件。

topol.tpr：运行输入文件，包含模拟开始时系统的完整描述。

confout.gro：结构文件，包含模拟最后一步的坐标和速度。

traj.trr：全精度轨迹，包含随时间变化的位置、速度和力。

traj.xtc：压缩轨迹，只包含坐标（低精度：0.001 nm）

ener.edr：能量文件，包含能量、温度、压力等随时间变化的相关参数

md.log：包含有关模拟的信息的日志文件，即性能、警告和错误。

### Quality Assurance

首先，必须确保模拟正确完成。许多变量会导致模拟崩溃，尤其是与力场相关的问题（如果您使用自定义参数）或系统平衡不足或不足。

```
gmx check -f p53_helix_CAH.xtc
```

关于模拟及其成功结论的另一个重要信息来源是日志文件。该文件的大部分内容包含有关模拟每个步骤的能量的信息。最后，有几个表格包含有关模拟性能的详细信息。

### Visually inspecting the simulation

尽管大多数分析归结为提取数据并绘制它们，但分子动力学首先是关于动力学的。因此，可以从轨迹中提取帧并将它们组合成电影。仅此一项就可以在整个模拟过程中充分了解肽的完整性。以下 Pymol 命令以sausage-like表示形式显示肽，从 N 端到 C 端按顺序着色。要操作轨迹文件，请使用 trjconv，GROMACS swiss-knife 实用程序。当要求选择要输出的组时，只选择蛋白质，否则你最终会得到一盒泥泞的水分子，掩盖了真正的动作！

```
gmx trjconv -f p53_helix_CAH.xtc -s p53_helix_CAH.tpr -o p53_helix_CAH-nojump.pdb -pbc nojump -dt 500
```

```
#pymol
cartoon tube
set cartoon_tube_radius, 1.5
as cartoon
spectrum count, rainbow, byres=1
smooth # Optional, for less jerky movie
unset movie_loop
mplay
```

肽在盒子周围移动，在它通过水分子扩散并超出盒子边界时摆动。当电影结束时，使用intra_fit命令对齐所有帧，以便更好地观察肽段运动。然后重播轨迹。

```
intra_fit name ca+c+n+o
zoom vis
mplay
```

随意使用 Pymol。放大特定区域，例如肽段最刚性或最柔韧的区域，并检查侧链构象（显示棒状图）。随意浪费一些（CPU）时间来制作漂亮的图像，使用 ray 和 png。请注意，过于复杂的场景可能会导致 Pymol 的内置光线追踪器崩溃，因此在这种情况下，您只能直接使用 png 在屏幕上获取图像。查看 Pymol Gallery 获取灵感，或向您的讲师寻求建议。如果你真的有很多时间可以浪费，你也可以制作轨迹的电影，尽管出于性能原因，这可能最好在课程的虚拟机之外完成。您可能需要从模拟中提取更多帧来制作相当大的电影，具体取决于您选择的帧速率。

```pymol
viewport 640, 480 # No HD, unless you really want to waste time!
set ray_trace_frames, 1
set ray_opaque_background, 0
mpng frame_.png
```

```
convert -delay 1 -loop 0 -dispose Background frame_*.png dynamics.gif
```

### Quantitative Quality Assurance

在对轨迹进行第一次目视检查后，假设模拟顺利进行，现在是对模拟质量进行额外和更彻底检查的时候了。该分析涉及测试热力学参数的收敛性，例如温度、压力以及势能和动能。有时，还根据每个帧的原子坐标相对于初始结构和/或平均结构的均方根偏差 (RMSD) 来检查模拟的收敛性。由于这个模拟是一个非常小且灵活的肽，预计它不会收敛，尽管可能会有惊喜！最后，还必须检查周期性图像之间相互作用的发生，因为如果确实发生了这些相互作用，它们可能会导致模拟中出现伪影。

### Convergence of the thermodynamical parameters

首先从能量文件中提取热力学参数，如前所述。感兴趣的是温度、压力、势能、动能、晶胞体积、密度和盒子尺寸。模拟的能量文件包含几十项。一些能量术语被分成组。这些组在 .mdp 文件中定义，可用于隔离系统的特定部分以供将来分析，例如，查看特定残基之间的相互作用。

```
gmx energy -f p53_helix_CAH.edr -o md_temperature.xvg
$MOLMOD_BIN/xvg_plot.py -i md_temperature.xvg
```

查看该图，了解温度如何在指定值 (310 K) 附近波动。系统的热容量也可以从这些波动中计算出来。必须从 .edr 能量文件中提取系统温度以及焓（对于 NPT）或 Etot（对于 NVT）值。此外，我们必须使用 -nmol 选项明确说明系统中有多少分子（您可以参考拓扑文件的末尾以获取系统中的分子总数）。这将允许 gmx energy 自动计算热容量并在其输出结束时显示。查看 GROMACS 手册了解更多详细信息。

```
gmx energy -f p53_helix_CAH.edr -fluct_props -nmol XXXX
```

某些项的平衡需要比其他项更长的时间。特别是，温度迅速收敛到其平衡值，而例如系统不同部分之间的相互作用可能需要更长的时间。

### Calculation of the minimum distance between periodic images

使用周期性边界条件的任何分子动力学模拟分析的一个关键点是检查相邻图像之间是否存在任何直接相互作用。由于周期性图像只是避免硬边界的一种技巧，因此此类交互是非物理的自交互，并使模拟结果无效。

```
gmx mindist -f p53_helix_CAH.xtc -s p53_helix_CAH.tpr -od md_mindist.xvg -pi
$MOLMOD_BIN/xvg_plot.py -i md_mindist.xvg
```

如果周期性图像瞄准非常短暂且不频繁，则可以忽略它的发生。如果它确实在模拟的一段时间内频繁或持续发生，则需要返回并重新进行整个设置。此外，不仅直接交互受到关注。如前所述，溶质周围的水与本体水具有不同的结构。为了安全起见，在计算允许的最小距离时增加一个额外的纳米。

### Conformational dynamics and stability I – Radius of Gyration

在分析任何结构参数之前，必须对轨迹进行按摩以避免由于周期性边界条件而造成的伪影。此外，如果轨迹仅包含必要的（蛋白质）原子及其信息，则所有分析工具的工作速度都会更快。

```
gmx trjconv -f p53_helix_CAH.xtc -s p53_helix_CAH.tpr -o p53_helix_CAH_reduced.xtc -pbc nojump
```

也许与这个特定的模拟并不完全相关，因为目标是对许多构象进行采样，但模拟质量保证的另一部分是检查结构本身的收敛性。这可以通过计算每帧原子坐标相对于初始结构或平均结构的均方根偏差 (RMSD) 来完成，也可以通过计算结构在轨迹上的回转半径来完成。回转半径表示分子的形状，并与实验获得的流体动力学半径进行比较。

```
gmx gyrate -f p53_helix_CAH_reduced.xtc -s p53_helix_CAH.tpr -o md_radius-of-gyration.xvg
$MOLMOD_BIN/xvg_plot.py -i md_radius-of-gyration.xvg
```

### Conformational dynamics and stability II – Root Mean Square Fluctuation (RMSF)

肽的结构在整个模拟过程中都会发生变化，但并不完全相同。一些区域比其他区域更灵活，通常是由于氨基酸序列的差异。均方根波动捕捉每个原子关于其平均位置的波动，并且通常对应于晶体温度（或 b）因素。将此实验测量与 RMSF 配置文件进行比较可以作为模拟的额外质量检查。温度因子越高，原子的移动性越强。此分析的一个有趣的附属品是平均结构的计算，可用于未来的分析

```
gmx rmsf -f p53_helix_CAH_reduced.xtc -s p53_helix_CAH.tpr -o md_rmsf.xvg -ox md_average.pdb -oq md_temperature-factors-residue.pdb -res
$MOLMOD_BIN/xvg_plot.py -i md_rmsf.xvg
```

### Conformational dynamics and stability III – Root Mean Square Deviation (RMSD)

由于 RMSF 的计算也产生了一个平均结构，现在可以计算整个轨迹的均方根偏差。该度量通常用作结构向平衡状态收敛的指标。 RMSD 是一种距离度量，因此对于低值最有意义。与平均结构相差 10Å 的两个框架很可能是完全不同的构象。 GROMACS 工具 rms 允许进行此类计算，特别是仅选择分子的特定原子组，例如骨架。

```
gmx rms -f p53_helix_CAH_reduced.xtc -s p53_helix_CAH.tpr -o md_bbrmsd_from_start.xvg
gmx rms -f p53_helix_CAH_reduced.xtc -s p53_helix_CAH.tpr -o md_aarmsd_from_start.xvg
gmx rms -f p53_helix_CAH_reduced.xtc -s md_average.pdb -o md_bbrmsd_from_average.xvg
gmx rms -f p53_helix_CAH_reduced.xtc -s md_average.pdb -o md_aarmsd_from_average.xvg
```

虽然相对于初始结构的 RMSD 是相关的，但如果它在相对较高的值处稳定，则不会告知构象的稳定性。如上所述，10Å 处的两种结构可能非常不同。出于这个原因，相对于平均结构的 RMSD 可能会提供一个更好的视角来了解整个模拟过程中结构变化的演变。

## Structural Analysis

When asked for a selection choose “Protein” if no selection is specifically stated or does not follow logically from the text.

Having assured that the simulation has converged to an equilibrium state, and that its results are likely to be valid, it is time for some real analysis that provides answers to a research question. Analysis of simulation data can be divided in several categories. One comprises the interpretation of single conformations according to some functions to obtain a value, or a number of values, for each time point. Example of these are the previously calculated RMSD and radius of gyration metrics. Next to that, the analysis can be done in the time domain, e.g. through averaging, such as (auto)correlations or fluctuations. In the next section, several different types of analyses will be performed, each providing a different but complementary view into the trajectory. Some might not be strictly necessary for this simulation, but are included as an example of what can be done elsewhere.

When running the GROMACS programs to perform the analyses, pay attention to their output as well as the plots they generate.

### Hydrogen Bonds

Secondary structures (of proteins) are maintained by specific hydrogen bonding networks. Thus, the number of hydrogen bonds, both internal and between the peptide and the solvent. The presence or absence of a hydrogen bond is inferred from the distance between a donor/acceptor pair and the H-donor-acceptor angle. OH and NH groups are regarded as donors, while O is always classified as an acceptor. N is an acceptor by default as well, unless specifically disabled. GROMACS can calculate hydrogen bonds over full trajectories with the `hbond` program. The program output informs on the number of hydrogen bonds, distance and angle distributions, and an existence matrix of all internal hydrogen bonds over all frames. The number of hydrogen bonds alone is a proxy for the existence of secondary structures.

Calculate the number of internal and protein-solvent hydrogen bonds over the trajectory. Note that for determining hydrogen bonds to the solvent the reduced trajectory cannot be used.gmx hbond -f p53_helix_CAH_reduced.xtc -s p53_helix_CAH.tpr -num md_hbond_internal.xvg
gmx hbond -f p53_helix_CAH.xtc -s p53_helix_CAH.tpr -num md_hbond_solvent.xvg
How does the number of internal hydrogen bonds correlate with the radius of gyration?Comment on the relation between the internal and the solvent hydrogen bond populations.

In addition to global analyses, many GROMACS programs support index files, which are created with the `make_ndx` program. These index files allow the creation of user-specified groups, such as single residues or stretches of residues. For example, it is possible to evaluate the creation of β-hairpins by checking the existence of hydrogen bonds between the two halves of the peptide. Assume you are working on a 14-residue long peptide. The syntax within `make_ndx` to create an index file to check for hydrogen bonds between the two halves is as follows:

```
r 1-7
name 19 half_1
r 8-14
name 20 half_2
q
```

Create an index file to assess the existence of hydrogen bonds that might justify a β-hairpin structure.gmx make_ndx -f p53_helix_CAH.tpr -o my_index.ndx
gmx hbond -f p53_helix_CAH_reduced.xtc -s p53_helix_CAH.tpr -n my_index.ndx -num beta_hairpin_hbond.xvg

On the basis of this analysis, is your peptide adopting a β-hairpin structure during the simulation?

### Secondary Structure

Among the most common parameters to analyse protein structure is the assignment of secondary structure elements, such as α-helices and β-sheets. One of the most popular tools for this purpose is the `dssp` software. Although not part of the GROMACS distribution, `dssp` can be freely obtained online at the [CMBI website](https://swift.cmbi.ru.nl/gv/dssp/), and integrated in many of its analysis tools. Specifically, the `do_dssp` tool produces a plot of the different secondary structure elements of each residue in the peptide as a function of time. This matrix, in *.xpm* format, can be converted into a Postscript file using the `gmx xpm2ps` tool, and then into a PDF file using `ps2pdf`. The `xpm2ps` utility allows a scaling flag, `-by`, that is useful for very short sequences, as well as a `-rainbow` flag that controls the coloring of the output.

Perform an analysis of the secondary structure of the peptide throughout the trajectory using *dssp*.gmx do_dssp -f p53_helix_CAH_reduced.xtc -s p53_helix_CAH.tpr -o md_secondary-structure.xpm
gmx xpm2ps -f md_secondary-structure.xpm -o md_secondary-structure.eps -by 20 -rainbow blue
ps2pdf md_secondary-structure.eps md_secondary-structure.pdfDiscuss the changes in secondary structure, if any.Compare and discuss the stability of the different secondary structures in the different conformations.

------

## Analysis of time-averaged properties

This simulation considers only one conformation. To obtain proper sampling of the peptide conformational landscape, 50 nanoseconds do not suffice. However, trajectories starting from different initial structures or starting from the same structure with a different initial random seed explore different regions of the conformational landscape. It is then desirable to combine different trajectories together and therefore obtain a much larger body of data.

Obtain different (full) trajectories from 2 of your colleagues. If possible, try to be as diverse as possible regarding initial structures.

### Preparation of a concatenated trajectory

The first step is to trim the trajectories in order to remove the first 10 nanoseconds, which can be conservatively considered as equilibration. This operation is possible through `trjconv` and its `-b` flag, which allows the user to specify an offset previous to which the trajectory data is ignored. To be able to extract only the peptide atoms, `trjconv` requires an *dummy* index file.

Trim the first 10 nanoseconds of each trajectory to discard the equilibration stage.gmx make_ndx -f p53_helix_CAH.tpr -o p53_helix_CAH.ndx
gmx trjconv -f p53_helix_CAH.xtc -s p53_helix_CAH.tpr -n p53_helix_CAH.ndx -pbc nojump -dt 50 -b 10000 -o p53_helix_CAH_reduced_10-50ns.xtc

Why doesn’t it matter which topology file is used to process the different trajectory files?

After all three trajectories are trimmed, they can be concatenated using the GROMACS program `trjcat`. Make sure to note down the order in which the trajectories are provided to `trjcat`. The concatenation requires two particular flags to be provided as input to the program: `-cat`, which avoids discarding double time frames, and `-settime`, which changes the starting time of the different trajectories interactively. Effectively, the second trajectory will start at 40 ns and the third at 80 ns. The program will prompt for an action during the concatenation: press `c`, which tells `trjcat` to append the next trajectory right after the last frame of the previous one.

Concatenate all three trajectories into a single one for further processing.gmx trjcat -f p53_helix_CAH_reduced_10-50ns.xtc p53_sheet_CAH_reduced_10-50ns.xtc p53_polypro_CAH_reduced_10-50ns.xtc -o p53_concatenated.xtc -cat -settime

### Root Mean Square Deviations – Part II

Although the root mean square deviation (RMSD) was already calculated to check for the convergence of the simulation, it can be used for a more advanced and in-depth analysis of conformational diversity. After all, the RMSD is a metric that compares structures. By performing an all-vs-all comparison with all frames in the concatenated trajectory, it is possible to identify groups of frames that share similar structures. This also quantifies the conformational diversity of a particular trajectory (or trajectories). The matrix allows also to detect and quantify the number of transitions between different conformations during the simulations. It is as relevant to have 10 different conformations or 2 that interconvert quickly. Since an all-vs-all RMSD matrix entails a very large number of pairwise comparisons, and the peptide conformations are different enough, use only backbone atoms to fit and calculate the RMSD.

Calculate and plot an all-vs-all RMSD matrix for the concatenated trajectory.gmx rms -f p53_concatenated.xtc -f2 p53_concatenated.xtc -s p53_helix_CAH.tpr -m p53_concatenated_RMSD-matrix.xpm
gmx xpm2ps -f p53_concatenated_RMSD-matrix.xpm -o p53_concatenated_RMSD-matrix.eps -rainbow blue
ps2pdf p53_concatenated_RMSD-matrix.eps p53_concatenated_RMSD-matrix.pdf
How many groups of similar structures do you see in the RMSD matrix?Are there overlapping regions of the conformational landscape in the different trajectories?

### Cluster Analysis

Using the all-vs-all RMSD matrix calculated in the previous step, it is possible to quantitatively establish the number of groups of similar structures that a trajectory (or concatenated trajectories) samples. Using an unsupervised classification algorithm, *clustering*, structures that are similar to each other within a certain RMSD threshold are grouped together. The size of a cluster, the number of structures that belong to it, is also an indication of how favourable that particular region of the conformational landscape is in terms of free energy. GROMACS implements several clustering algorithms in the `cluster` program. Here, we will use the `gromos` clustering algorithm with a cutoff of 2 Å. Briefly, the algorithm first calculates how many frames are within 2 Å of each particular frame, based on the RMSD matrix, and then selects the frame with the largest number of neighbors to form the first cluster. These structures are *removed* from the pool of available frames, and the calculation proceeds iteratively, until the next largest group is smaller than a pre-defined number. The `cluster` program produces a very large number of output files that inform on several different properties of the clusters. Importantly, it also produces a PDB file with the centroids, or representatives, of each cluster.

Cluster the RMSD matrix using the GROMOS method to quantitatively extract representative structures of the simulation. Choose peptide backbone for fitting and all-atoms of peptide as output. This is important, since we have will use the output structures for docking.

gmx cluster -f p53_concatenated.xtc -s p53_helix_CAH.tpr -dm p53_concatenated_RMSD-matrix.xpm -dist p53_concatenated_rmsd-distribution.xvg -o p53_concatenated_clusters.xpm -sz p53_concatenated_cluster-sizes.xvg -tr p53_concatenated_transitions.xpm -ntr p53_concatenated_transitions.xvg -clid p53_concatenated_cluster-id-over-time.xvg -cl p53_concatenated_clusters.pdb -cutoff 0.2 -method gromos

How many clusters did the algorithm find? Tune the cutoff to obtain a reasonable number of clusters (e.g. 10-15).

What is the clustering cutoff that allows the definition of that number of clusters? Do you think these clusters are meaningful, i.e. contain only similar structures?

Open the resulting PDB file in Pymol and compare the centroids of each cluster with the others.

disable all
intra_fit name ca+n+c+o
split_states p53_concatenated_clusters
delete p53_concatenated_clusters
dssp all, [PATH TO DSSP e.g. /opt/bin/dssp]
as cartoon

Are there any meaningful differences between the largest clusters?

------

## Picking representatives of the simulation

The aim of this simulation exercise was the sample the conformational landscape of the p53 N-terminal transactivation peptide, in order to extract representatives that could be used to generate models of its interaction with the MDM2 protein. The last step of clustering provides an unbiased method to select structures that were sampled throughout most of the trajectory (large clusters) and are likely good candidates for seeding the docking calculations.

Select 5 representatives of the clusters you obtained in the previous step and create individual PDB files using Pymol.



# 3.HADDOCKing of the p53 N-terminal peptide to MDM2

## A bite of theory

蛋白质-蛋白质相互作用介导细胞中的大多数细胞过程，例如分化、增殖、信号转导和细胞死亡。然而，它们的结构表征并不总是微不足道的，即使 X 射线晶体学和核磁共振光谱学不断发展。罪魁祸首千差万别，从可能使它们难以纯化或结晶的复合物的天然环境到系统的规模太大，当前的方法无法掌握。更重要的是，体内平衡通常取决于非常严格调控的瞬态相互作用，这是这些复合物纯化和表征的另一个障碍。尽管如此，结构生物学三十年的惊人进展表明，蛋白质折叠空间是有限的，折叠游戏的规则也很明确。这促使开发了几种计算方法，旨在补充实验技术，以寻求蛋白质组的完整 3D 视图。用于预测蛋白质-蛋白质复合物的广泛使用的计算方法是分子对接，其目的是从其天然成分的结构（或模型）开始生成这种复合物的结构

本教程将介绍 HADDOCK（High Ambiguity Driven DOCKing）作为一种方法来预测蛋白质 - 蛋白质复合物的三维结构，使用各种信息源来指导对接过程并对预测模型进行评分。 HADDOCK 是一组源自 ARIA 的 Python 脚本，它们利用 CNS（Crystallography and NMR System）的功能来计算分子复合物的结构。 HADDOCK 与其他对接软件的不同之处在于它继承自 CNS 的能力，能够将实验数据作为约束并使用这些数据来指导对接过程以及传统的能量学和形状互补。此外，与 CNS 的密切耦合使 HADDOCK 能够实际生成足够质量的模型，以便在蛋白质数据库中存档。

HADDOCK 的一个核心方面是ambiguous interaction restraints或 AIR 的定义。这些允许将原始数据（如 NMR 化学位移扰动或诱变实验）转换为距离限制，这些限制被纳入计算中使用的能量函数中。空气是通过分为两类的残留物清单来定义的：主动和被动。通常，活性残基是相互作用的核心残基，例如敲除消除相互作用的残基或化学位移扰动较高的残基。在整个模拟过程中，如果可能，这些活性残基被限制为界面的一部分，否则会导致得分惩罚。被动残基是那些有助于相互作用但被认为不太重要的残基。如果这样的残基不属于界面，则没有得分惩罚。因此，仔细选择哪些残基是主动的，哪些是被动的，对于对接的成功至关重要。

HADDOCK 的对接协议被设计成使分子经历不同程度的灵活性和不同的化学环境，它可以分为三个不同的阶段，每个阶段都有一个明确的目标和特征：

**1. Randomization of orientations and rigid-body minimization (it0)**

在这个初始阶段，相互作用的伙伴被视为刚体，这意味着所有几何参数，如键长、键角和二面角都被冻结。partner在空间中分开并围绕它们的质心随机旋转。接下来是刚体能量最小化步骤，其中允许合作伙伴旋转和平移以优化交互。 AIR 在此阶段的作用尤为重要。由于它们包含在最小化的能量函数中，因此产生的复合物将偏向于它们。例如，定义一组非常严格的 AIR 会导致构象空间的采样非常狭窄，这意味着生成的姿势将非常相似。相反，非常稀疏的约束（例如伴侣的整个表面）将导致非常不同的解决方案，在绑定区域显示更大的可变性。

![img](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-protein-protein-basic/haddock_mini.gif)

**2. Semi-flexible simulated annealing in torsion angle space (it1)**

对接协议的第二阶段通过基于分子动力学的三步改进为交互伙伴引入了灵活性，以优化界面包装。值得注意的是，扭转角空间的灵活性意味着键长和角度仍处于冻结状态（二面角变化）。相互作用的伙伴首先保持刚性，仅优化它们的方向。然后在界面中引入了灵活性，该界面是根据对 5Å 截止值内的分子间接触分析自动定义的。这允许来自 it0 的不同绑定姿势定义不同的灵活区域。然后允许属于该界面区域的残基在第二个细化步骤中移动它们的侧链。最后，灵活接口的主链和侧链都被赋予了自由。 AIR 在此阶段再次发挥重要作用，因为它们可能会驱动构象变化。

![img](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-protein-protein-basic/haddock_sa.gif)

**3. Refinement in Cartesian space with explicit solvent (water)**

注意：这个阶段是标准 HADDOCK 协议的一部分，直到（包括）v2.2。从 v2.4 开始，它不再默认执行，但用户仍然可以选择启用它。取而代之的是执行短暂的能量最小化。对接协议的最后阶段将复合物浸入溶剂壳中，以改善相互作用的能量学。 HADDOCK 目前支持水（TIP3P 模型）和 DMSO 环境。后者可用作膜模拟物。在这个简短的显式溶剂改进中，模型在 300K 下进行了简短的分子动力学模拟，对非界面重原子进行了位置限制。这些限制后来被放宽以允许优化所有侧链。

![img](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-protein-protein-basic/haddock_water.gif)

该协议的性能当然取决于每一步生成的模型数量。很少有模型不太可能捕捉到正确的绑定姿势，而夸大的数字在计算上会变得不合理。标准的 HADDOCK 协议在刚体最小化阶段生成 1000 个模型，然后在 it1 和水中优化最好的 200 个——关于能量函数。但是请注意，虽然默认情况下在 it0 中生成了 1000 个模型，但它们是五次最小化试验的结果，并且对于每个模型，还对 180º 对称解进行采样。实际上，写入磁盘的 1000 个模型是 10.000 个对接解决方案的抽样结果。

最终模型基于特定的相似性度量自动聚类 - 位置界面配体 RMSD (iL-RMSD) - 通过拟合受体（第一个分子）的界面并计算界面上的 RMSD 来捕获界面的构象变化较小的合作伙伴。此计算中使用的接口是根据对所有模型中所有接触的分析自动定义的。

## Predicting the interface of p53 on Mdm2

鉴于有某种信息可以指导对接，HADDOCK 擅长预测蛋白质复合物的结构。在缺乏实验信息的情况下，可以使用序列保守性和表面残基的生物物理特性等特征来推断蛋白质表面上的假定界面。由于同源建模模块创建了小鼠 MDM2 的同源物列表，因此可以评估哪些残基更保守。

首先我们需要再次找到序列同源物。这次我们将使用 UniProt 运行 BLAST 搜索。我们可以回到同源建模部分的条目，在那里我们在 Uniprot 中查找鼠标 MDM2。

应该会出现一个熟悉的页面，其中包含所有先前描述的信息。直接转到“序列”部分。您看到的序列是canonical序列，这意味着它要么是最普遍、与其他物种中的直系同源序列最相似的序列，要么是在没有任何信息的情况下最长的序列。在序列的右侧，可以运行 BLAST（基本局部比对搜索工具）搜索。

接下来，将打开一个带有 BLAST 搜索的新窗口。可以输入蛋白质或核苷酸序列或 UniProt 标识符。

选择所有序列，然后单击“对齐”部分中的“对齐”。运行完成后，下载 FASTA 格式的压缩对齐。

![image-20211013163920916](https://i.loli.net/2021/10/13/qFtZCAEKoJp4Sms.png)

为了可视化比对，以及哪些位置更保守，最简单的方法是生成一个序列标志。对于序列中的每个位置，标识标识最常出现的残基并根据保护分数缩放其单字母代码。我们将使用 [weblogo](http://weblogo.berkeley.edu/logo.cgi) 服务器，以便为 BLAST 生成的比对生成序列logo。

除了序列保守性，其他特征可用于预测蛋白质结构上可能的界面。例如，某些残基往往在蛋白质-蛋白质界面处过多。该信息与进化守恒和表面聚类算法相结合，该算法可以找到满足上述两个标准的表面残基组，从而产生相当准确的预测。CPORT 是一种预测蛋白质-蛋白质界面残基的算法。它将六种界面预测方法组合成一个共识预测器。预测可以用作 HADDOCK 中的主动和被动残基。 [CPORT](https://alcazar.science.uu.nl/services/CPORT/)旨在为 HADDOCK 提供预测。服务器还返回一个原始结构的 PDB 文件，其中加载了温度因子列中的预测。这对于可视化 Pymol 中的预测非常有帮助。

```
Predicted residues (active residues in HADDOCK):
14, 16, 33, 40, 41, 55, 58, 61, 62, 63, 64, 65, 66, 67, 68,
69, 70, 71, 72, 73, 74, 76, 78, 80, 83, 89, 90, 92, 93, 108,
111,

Surrounding residues (passive residues in HADDOCK):
11, 12, 13, 15, 17, 18, 19, 20, 21, 27, 29, 30, 31, 32, 35,
36, 39, 42, 51, 52, 54, 59, 79, 81, 84, 86, 87, 88, 94, 95,
96, 98, 104, 105, 106, 107, 109, 110,
```

![image-20211013164542930](https://i.loli.net/2021/10/13/KXZlMa21he58uEc.png

## Preparing the structures for the docking calculation

为了使用 HADDOCK 执行对接计算，MDM2 和 p53 的初始结构必须满足一些要求。首先，PDB 文件的最后一行必须有一个 END 语句。文件也不能包含多个占用的原子。然而，可以提交一个结构集合，这对于从分子动力学模拟中提取的 p53 的代表很有用。提交这样的集合时，所有成员必须包含完全相同的原子。幸运的是，SWISS-MODEL 和 GROMACS 生成的文件都遵守 END 语句和多占用要求，因此这里不需要任何操作。

源自分子动力学模拟的 p53 肽的结构可以作为一个整体提交给 HADDOCK。 pdb-tools 集的 pdb_mkensemble 实用程序提供了一种从隔离的 PDB 文件构建这种集成结构的快速方法。它还向 PDB 文件添加了正确的 END 语句。最后，它具有对集成完整性的内置检查，即所有成员都具有完全相同的原子构成。

```
pdb_mkensemble p53_cluster_1.pdb p53_cluster_2.pdb p53_cluster_3.pdb > p53_ensemble.pdb
```

## Setting up the docking calculation using the HADDOCK web server

要获得 HADDOCK 帐户，请访问 https://wenmr.science.uu.nl/haddock2.4/ 并单击注册。
准备好初始结构并构建了一个假定的界面残基列表后，是时候使用 HADDOCK 2.4 Web 服务器界面提交对接计算了。
在这里，您可以获得有用的信息，例如指向新 HADDOCK 最佳实践指南的链接，其中包含不同对接场景的设置。在服务器信息下，您可以找到网络服务器的默认设置，这对于了解如何处理限制很重要，以及支持的修饰氨基酸列表和当前 HADDOCK 版本。

要开始提交作业，请单击提交新作业。



### Submission and validation of structures

For this we will make us of the [HADDOCK 2.4 interface](https://wenmr.science.uu.nl/haddock2.4/submit/1) of the HADDOCK web server.

步骤 1：在“作业名称”字段中定义对接运行的名称，例如MDM2-p53。
步骤 2：选择要停靠的分子数，本例中默认为 2。
步骤 3：输入第一个蛋白质PDB文件。对于此展开分子 1 - 如果尚未展开，则输入。

第一个分子：提供的结构在哪里？ -> “我正在提交”

使用哪条链？ -> 全部（对于这种特殊情况）

PDB结构提交->浏览选择modelx.pdb（你之前准备的同源模型）

注意：将所有其他选项保留为默认值。

步骤 4：输入第二个蛋白质 PDB 文件。这是三种肽构象的集合。对于此展开分子 2 - 如果尚未展开，则输入。

第一个分子：提供的结构在哪里？ -> “我正在提交”
使用哪条链？ -> 全部（对于这种特殊情况）
要提交的 PDB 结构 -> 浏览并选择 p53_ensemble.pdb（包含您之前准备的肽模型集合的 PDB 文件）

由于我们的同源模型和肽不对应于完整序列，因此最好具有不带电荷的末端。

步骤 4：点击界面左下角的下一步按钮。这会将结构上传到 HADDOCK 网络服务器，在那里它们将被处理和验证（检查格式错误）。服务器利用 Molprobity 检查侧链构象，最终交换它们（例如天冬酰胺）并定义组氨酸残基的质子化状态。

### Definition of restraints

如果一切顺利，界面窗口应该已经更新，现在应该显示分子 1 和 2 的残基列表。我们将利用每个分子残基序列下方的文本框来指定活性残基列表用于对接运行。约束的定义确实需要一些思考。 HADDOCK 中的活性残基是那些需要在界面上的残基。另一方面，被动残基是那些可能在界面上的残基。不明确的交互约束或 AIRs 是在一个伙伴的每个主动残基与另一个伙伴的主动和被动残基的组合之间创建的。不在界面上的活性残基会导致能量损失，而被动残基则不是这种情况。对于 MDM2 和 p53 的对接，MDM2 上的活性残基取自 CPORT 预测，而肽仅定义为被动残基。这遵循在我们的 Structure 2013 论文中发布的配方。以这种方式，蛋白质的活性残基将吸引肽，而肽残基本身并不具有所有接触。

步骤 6：指定第一个分子的活性残基。为此展开分子 1 - 参数（如果尚未展开）。在这个阶段，我们将利用 CPORT 返回的活性残基用于 MDM2

活性残基（直接参与相互作用）-> 在此处输入 CPORT 为 MDM2 返回的活性残基列表

在主动残基周围自动定义被动残基 -> 取消选中（只有在为第二个分子定义了主动残基时才应定义被动）



步骤 7：指定第二个分子的活性残基。为此展开分子 2 - 参数（如果它尚未展开）

活性残基（直接参与相互作用）-> 留空（在这种情况下肽没有活性）
在主动残基周围自动定义被动残基 -> 取消选中（默认选中）
被动残基（周围表面残基）-> 在此处输入肽的所有残基作为逗号分隔的列表

### Definition of fully flexible segments

第 8 步：由于肽具有高度的灵活性，我们将为肽提供更大的灵活性以允许更大的构象变化。为此展开分子 2 的完全柔性段选项卡并输入：

Fully flexible segments -> Enter here all residues of the peptide as a comma-separated list

在这里您也可以简单地再次选择整个肽序列。这将导致 HADDOCK 认为肽残基在模拟退火精炼阶段的所有阶段都是完全灵活的，因此增加了采样。

### Increased sampling

- **Step 9:** This can be done in one simple step by choosing the bioinformatics predictions settings described [here](https://wenmr.science.uu.nl/haddock2.4/settings#bioinfo).

Optimize run for bioinformatics predictions -> **check**

 

- **Step 10:** Click on the **Next** button on the bottom of the page.

检查生物信息学预测设置更改自动采样到这些参数：

Number of structures for rigid body docking -> 10000

Number of structures for semi-flexible refinement -> 400

Number of structures for the final refinement -> 400

Number of trials for rigid body minimisation -> 1

Number of structures to analyze -> 400

使用一组结构转化为在刚体阶段每个构象的采样更差。每个起始构象将被采样有限的次数，定义为在刚体对接阶段采样的模型总数除以集合中的模型数。使用 1000 个刚体模型的默认值，10 个起始构象的集合转换为每个集合成员生成的 100 个模型。每个模型采样的减少可能会降低对接计算的准确性，特别是如果限制是模糊的，就像使用生物信息学预测时的情况一样。出于这个原因，建议增加在对接协议的各个阶段生成的结构数量。根据经验，合奏中的每个成员有 1000 个刚体模型是一个不错的数字。选择到它1和水的模型数量可以简单地增加一倍。这些细化阶段的计算成本不允许按比例增加。这些数字可以在采样参数选项卡中一一编辑，但在检查生物信息学预测设置时会自动完成。

注意：由于在初始结构集合的情况下每个模型的采样减少，建议限制初始集合中的构象数量。如果每个分子对接都使用构象集合，这一点就更加重要。例如，在两个分子对接的情况下，每个分子 10 个模型将生成 100 个组合。刚体对接阶段采样10000次，每个组合只会采样100次。请注意，服务器将 it0 模型的数量限制为最大 10000。

### Clustering parameters

- **Step 11:** For this unfold the **Parameters for clustering menu**.

HADDOCK 提供两种不同的聚类算法。有关详细信息，请参阅在线手册。对于肽和小分子，我们建议使用 RMSD 聚类。还必须调整聚类算法以适应小尺寸的肽。 7.5Å（界面-配体RMSD）的默认截止值针对蛋白质-蛋白质对接进行了优化，并且在蛋白质-肽复合物的情况下很可能太大。使用此值进行聚类很可能会生成非常大且多样化的聚类。因此，我们应该减少聚类截止：

Clustering method (RMSD or Fraction of Common Contacts (FCC)) -> Select RMSD

RMSD Cutoff for clustering (Recommended: 7.5A for RMSD, 0.75 for FCC) -> 5.0



### Automatic restraining of secondary structure elements



- **Step 12:** For this unfold the **Dihedral and hydrogen bonds restraints menu**.

HADDOCK 提供了一个选项，可以根据输入结构自动定义二面角约束。这可以应用于整个序列，或仅应用于 alpha 螺旋段或 alpha 和 beta 段。这些是根据测量的二面角自动检测的。对于柔性肽，由于我们将它们视为完全柔性的，因此建议打开此选项。

Automatically define backbone dihedral restraints from structure? -> Select Only for alpha helix



### Advanced sampling parameters

- **Step 13:** Adjust the number of flexible refinement steps to increase the sampling of peptide conformations in **Advanced sampling parameters menu**.

number of MD steps for rigid body high temperature TAD -> 1000

number of MD steps during first rigid body cooling stage -> 1000

number of MD steps during second cooling stage with flexible side-chains at interface -> 2000

number of MD steps during third cooling stage with fully flexible interface –> 2000



### Job submission

这个接口允许我们修改许多控制 HADDOCK 行为的参数，但在我们的例子中，默认值都是合适的。它还允许我们下载对接运行的输入结构（以 tgz 存档的形式）和一个 haddockparameter 文件，其中包含我们运行的所有设置和输入结构（以 json 格式）。我们强烈建议下载此文件，因为它可以让您在上传到 HADDOCK 网络服务器的文件上传界面后重复运行。它可以作为运行的输入参考。例如，还可以编辑此文件以更改一些参数

- **Step 14:** Click on the **Submit** button at the bottom left of the interface.

第二个链接指向模拟的结果页面。由于常规对接模拟平均持续几个小时，因此此页面显示其当前状态，但未完成为 PROCESSING、QUEUED 或 RUNNING。如果严重错误阻止模拟继续，无论是由于输入数据的问题，还是模拟本身的问题，网页都会显示错误消息。大多数这些状态更改都伴随着一封电子邮件，该电子邮件发送到链接到用户帐户的地址。如果出现错误，此电子邮件还会提供有关原因的其他详细信息。对于学生，由于所有帐户都已预先配置，因此电子邮件通知已关闭。



## Analyzing the docking calculation results

模拟完成后，将生成结果页面并向用户发送通知电子邮件。此结果页面包括前十个集群的概述，按其四个最佳结构的平均 HADDOCK 得分排名，包括每个集群的能量术语和其他结构度量的统计数据。这允许快速评估生成模型的质量。此外，在页面底部，模型的不同质量控制措施的几幅图也提供了一眼检查对接模拟质量的机会。



## Visual inspection of the cluster representatives

任何分子模拟，包括对接，都缺乏产生单一良好模型的准确性。但是，通过足够的尝试，合理的模型很可能会填充结果。特别是 HADDOCK，鉴于其数据驱动的特性，如果数据质量足够好，可以生成质量更高的模型。在模拟结束时，所有模型都聚集在一起，以过滤掉与模型池中很少有其他模型相似的孤立结构。然后在集群基础上分析模型，并且最佳集群中的最佳模型被认为是模拟的最佳代表模型。尽管如此，仔细的目视检查至关重要，并且应该始终是任何模拟分析的第一步。为了更好地可视化每个集群之间的差异，将所有模型叠加到最大伙伴的主干上，在这种情况下是 MDM2 分子。为此，我们必须选择一个参考模型，我们将其任意定义为 cluster1_1 的模型。



```
select refe, cluster1_1 and chain A
align cluster2_1, refe
align cluster3_1, refe
…
```


## 蛋白质相关模拟分析


#### 周期性校正

主体模拟的轨迹得到之后，有时候会发现体系中的蛋白或者配体有跨盒子的现象，这时候就需要对轨迹的周期性边界条件（PBC）进行校正。一来可以避免分析过程中出现意外情况（比如RMSD曲线有较大的突跃和突降等），二来也是可以获得较为美观的轨迹，方便我们直接观察。

1. 蛋白配体PBC校正一般过程

gmx trjconv的参数多种多样，可以有多种校正轨迹的组合，自己使用的时候，为了获得好的校正效果，可以在理解参数的基础上自由组合尝试。对蛋白配体体系，我个人一般使用的流程为：

- 选取某一个原子，以其对体系居中
- 保证体系中分子的完整性
- 消除平动和转动

下面简单叙述我校正周期性的具体步骤。

1.1 居中

首先生成一个体系的索引文件，并在里面添加 center 组，组内只有一个原子，之后会以这个原子为盒子中心对体系进行居中。

```bash
# 把蛋白配体组合到一个组
gmx make_ndx -f md.tpr -o prolig.ndx
# 手动添加一个 center 组
echo -e "\n[ center ]\n500\n" >> prolig.ndx
```

这里例子里选择的是序号为500的原子。

比如说，对于蛋白而言，可以先计算一下蛋白的几何中心，然后选择一个距离蛋白质几何中心最近的原子作为中心原子对体系进行居中。

这里我们以npt.gro文件为例，寻找一个可以对蛋白进行居中的中心原子。可以使用DuIvyTools的`find_center`命令来进行中心原子的寻找，这个命令需要两个输入文件，一个是模拟的gro文件，一个是体系的索引文件。

```bash
dit find_center -f npt.gro index.ndx
```

输入这个命令之后需要选组，比如说要对蛋白质居中的话，就选蛋白质就行了，这个命令可以找到蛋白质内的最靠近中心的原子的序号。

最后输出的信息，就是离中心最近的原子了，在第三列写着的就是这个原子的序号。把这个序号写到 center 组里面就可以了。

执行下面的命令对体系进行居中：

```bash
gmx trjconv -s md.tpr -f md.xtc -o center.xtc -n prolig.ndx -pbc atom -center
# 第一次输入要居中的组，选 center 组
# 第二次输入要输出的组，选你自己定义的蛋白配体组
```

1.2 分子完整

对居中之后的轨迹做一遍保证分子完整性的校正。当然可以用whole选项保证分子完整性，我这里用下面这条命令：

```bash
gmx trjconv -s md.tpr -f center.xtc -o mol.xtc -n prolig.ndx -pbc mol -ur compact
# 选择蛋白配体组即可
```

1.3 fit一下

如果幸运的话，做完上述两步应该就可以了，我个人习惯还对蛋白去除一下平动和转动，方便查看配体和蛋白质之间的相互运动：

```bash
gmx trjconv -s md.tpr -f mol.xtc -o fit.xtc -n prolig.ndx -fit rot+trans
# 第一次输入，选择对蛋白质进行fit
# 第二次输入，选择蛋白配体组进行轨迹输出
```

这里得到的fit.xtc就是周期性校正之后的蛋白质（或者蛋白质配体复合物）的轨迹了，可以用于后续分析。有些命令执行的时候还需要tpr文件和xtc文件里面的原子数目一致，因而也需要对tpr文件做一下处理使之和xtc文件的原子数目一致。可以使用`gmx convert-tpr`命令完成这个操作：

```bash
gmx convert-tpr -s md.tpr -o fit.tpr -n index.ndx
```

选择和前面生成xtc的同样的组就行了，如此就可以保证tpr和xtc有同样的原子数目。


不管校正到哪一步了，都可以把轨迹可视化出来，自己visual check一下。我一般习惯把xtc转成pdb，用pymol查看。

```bash
gmx trjconv -s md.tpr -f fit.xtc -o fit.pdb -dt 1000 -n prolig.ndx 
```

加个-dt 1000，选择合适的时间间隔，防止输出的pdb文件太大。


#### 动态互相关矩阵(dynamics cross-correlation matrix, DCCM)

`gmx covar`命令用于计算并对角化（质量加权）的协方差矩阵，所有结构都会叠合到结构文件中的参考结构。计算得到的本征向量会写到轨迹文件中。而计算得到的协方差矩阵会通过-xpm和-xpma写到xpm文件中。-xpma得到的xpm文件包含了每个原子对三个坐标的协方差的总和，因而是$N*N$的，N为原子数目。而-xpm得到的xpm则是包含每一维坐标的，因而是$3N*3N$的。要做后处理，当然或许可以直接读取xpm文件的数据然后计算互相关指数得到DCCM；也可以通过covar命令的-ascii参数输出整个的协方差矩阵到一个文本文件($3N∗3N$)，如此得到的数据或许会稍微精准一些。

首先我们需要利用covar命令得到协方差矩阵：

```bash
gmx covar -s md.tpr -f md.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm -ascii covar.dat
```

按照需要选择对齐的组和计算的组，之后协方差矩阵会被保存到xpm文件中，协方差矩阵的数据信息会被保存到dat文件中。

之后可以使用DuIvyTools将covar.dat处理成dccm.xpm：

```bash
dit dccm_ascii -f covar.dat -o dccm.xpm
```

之后使用DuIvyTools对dccm.xpm进行可视化：

```bash
dit xpm_show -f dccm.xpm -o dccm.png -zmin -1 -zmax 1 -cmap bwr 
```

#### 残基距离接触矩阵(residue distance contact matrix)

使用`gmx mdmat`可以得到模拟轨迹的平均残基距离接触矩阵，也可以通过`-b`和`-e`等来设置分析时间：

```bash
gmx mdmat -f md.xtc -s md.tpr -mean rdcm.xpm
```

之后使用DuIvyTools将rdcm.xpm可视化：

```bash
dit xpm_show -f rdcm.xpm -o rdcm.png
```

#### 主成分分析(principal component analysis, PCA)

要获得主成分，需要对轨迹进行主成分分析，gmx提供了`covar`命令和`anaeig`命令帮助我们进行相关分析。请一定注意要在PCA之前消除轨迹的平动和转动，以免分子的整体运动影响分子内部运动的分析。

`covar`命令的作用是对轨迹进行协方差矩阵和本征向量的计算。

```bash
gmx covar -s md.tpr -f md.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covapic.xpm 
```

- eigenvalues.xvg里面记录了分析得出的多个本征值的序号和大小
- eigenvectors.trr即是投影到本征向量之后的轨迹
- covapic.xpm即是轨迹的协方差矩阵

命令执行之后会让你选用哪些原子做align，选C-alpha或者backbone，之后会让你选对哪些原子做PCA，也选C-alpha或者backbone，可以两次选的不一样，按自己需求来就行。当然你也可以选择其它的索引组来进行align和PCA，这里这样子只是因为这俩原子数稍微少一些，分析会快一些。

eigenvalues.xvg里面存储的按本征值大小排序的多个本征值。我们通常希望前两个本征值的大小的和可以越大越好，这意味着前两个主成分就可以代表蛋白的大部分运动信息。因而分析PCA的时候可以将eigenvalues.xvg的第二列的总和算出来，再算一下前三个本征值的大小占总和的比，给出前三个主成分的各自占比以及占比总和。

covapic.xpm存储的是协方差矩阵，同样是xpm文件，可以用`dit xpm_show`命令查看。

之后我们需要利用anaeig命令将轨迹投影到前两个主成分上，也即生成pc1.xvg和pc2.xvg。

```bash
gmx anaeig -s pro.tpr -f pro20.xtc -v eigenvectors.trr -first 1 -last 1 -proj pc1.xvg
gmx anaeig -s pro.tpr -f pro20.xtc -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg
```

#### 自由能形貌图(free energy landscape, FEL)

**方法一： 利用RMSD和Gyrate绘制FEL**

首先获得蛋白质的RMSD(对齐到backbone)和蛋白质的回旋半径数据并存储到rmsd.xvg和gyrate.xvg中：

```bash
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
# 选择backbone进行对齐，选择Protein计算和输出
gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
# 选择Protein进行计算
```

要使用`gmx sham`计算FEL,还需要将数据整理一下，第一列是时间、第二列是rmsd、第三列是回旋半径（gyrate.xvg里面包含了四列回旋半径的数据，分别是总体的，XYZ三个方向上各一列，这里需要的是总的回旋半径列）。

可以使用DuIvyTools工具进行rmsd和gyrate数据的组合：

```bash
dit xvg_combine -f rmsd.xvg gyrate.xvg -c 0,1 1 -o sham.xvg
```

得到组合文件之后，就可以利用gmx的sham命令来生成自由能形貌图了:

```bash
gmx sham -tsham 310 -nlevels 100 -f output.xvg -ls gibbs.xpm -g gibbs.log -lsh enthalpy.xpm -lss entropy.xpm
```

简单介绍下这些参数：

-tsham ： 设定温度
-nlevels： 设定FEL的层次数量
-f : 读入组合文件
-ls : 输出自由能形貌图（Gibbs自由能）
-g ：输出log文件
-lsh : 焓的形貌图
-lss : 熵的形貌图

在上面的命令执行完成之后，会输出好几个文件，除了上文提到的，还包括一个index文件(bindex.ndx) 和一个ener.xvg文件。

gibbs.xpm即是最关心的是Gibbs自由能的形貌图。可以使用DuIvyTools对其进行绘图：

```bash
dit xpm_show -f gibbs.xpm -m countour -cmap jet -o fel.png
```

在完成了自由能形貌图的绘制之后，可能需要寻找FEL中最低能量对应的蛋白质构象。

这里需要用到gibbs.log，bindex.ndx这两个文件。gibbs.log中记录了索引组与能量的关系，bindex.ndx则记录了索引组与帧数的关系。如下所示：

gibbs.log部分截取如下：

```gibbs.log
Minimum 0 at index 26 energy      7.589
Minimum 1 at index 46 energy      7.589
Minimum 2 at index 50 energy      7.589
Minimum 3 at index 56 energy      7.589
Minimum 4 at index 141 energy      5.803
```

bindex.ndx部分截取如下：

```bindex.ndx
[ 26 ] # 这是索引
1274   # 这是索引对应的时间帧
[ 46 ]
2
[ 141 ]
4
1282
```

那么我们如何找最低能量的构象呢？假如说Minimum 4 是能量最低点，它的索引组是141，那么我们就到bindex.ndx中找到索引[ 141 ]，查看到这个索引对应的时间帧是第4帧和第1282帧，然后我们可以根据这个帧数去计算对应的模拟的时刻，从而使用`gmx trjconv`命令提取对应时刻的蛋白质构象。假如说用于分析的数据是从0ns开始，每一帧的时间间隔为10ps,那么第4帧对应的时刻即为30ps。`gmx sham`对于帧的计数是从1开始的，也即这里的0ns即为第1帧。需要注意的是，假如说获取用于绘制FEL的数据或者执行`gmx sham`的时候对数据的时间做了设置，比如说获得的RMSD和gyrate是从10ns开始的，或者说执行`gmx sham`的时候使用了`-b -e`等参数，那么计算时间帧的时候也需要仔细考虑到这些因素，算出正确的时间帧。

数据的开始时刻和时间间隔都可以从`sham.xvg`里面看到，所以实际上第几帧也就是`sham.xvg`的数据的第几行，那么对应的从第几行就可以获得对应的时间帧，这也是一种简洁的时间帧获得方法。

假如说要提取时刻为1010ps的构象：

```bash
gmx trjconv -f md.xtc -s md.tpr -b 1010 -e 1010 -o protein.pdb
```

如此就可以获得了对应时刻的构象了。


**方法二： 利用PCA绘制FEL**

依照上述主成分分析方法获得前两个主成分的xvg文件之后，使用DuIvyTools组合两个文件以获得sham.xvg。

```bash
dit xvg_combine -f pc1.xvg pc2.xvg -c 0,1 1 -o sham.xvg
```

之后按照方法一中的描述利用`gmx sham`命令得到gibbs.xpm并用DuIvyTools可视化即可。











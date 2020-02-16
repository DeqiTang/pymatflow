Pymatflow使用指南
===================
Pymatflow适合哪些人使用？

pymatflow的开发目的是为了实现大规模的高通量计算以及提高计算准备的效率。但这一切的前提是使用者要对背后的计算软件如quantum espresso、SIESTA、CP2K等有较为熟练的掌握，以保证合理使用本工具。

目前Pymatflow仅支持Linux平台，暂无支持Windows平台的打算。

## 简单示例

这里以一个简单的静态计算作为示例，展示如何使用Pymatflow来加速计算的准备工作:

```
~$ qe-scf.py -f CO2.xyz
```

即可自动生成用于QE的静态SCF输入文件并运行。

这里CO2.xyz有一点需要注意的是其中第二行需要以一定格式来定义cell以被Pymatflow读取:

```
3
cell: 10 0 0 | 0 10 0 | 0 0 10
C 0.0000000 0.000000 0.000000
O -1.160000 0.000000 0.000000
O 1.1500000 0.000000 0.000000
```

## Pymatflow支持的计算程序

Pymatflow可以自动生成许多常见的第一性原理计算程序，比如Quantum Espresso、CP2K、SIESTA、Abinit。未来将支持更多的程序。

## Pymatflow的安装

```
pip install git+https://gitlab.com/deqitang/pymatflow
```

## Pymatflow的定位

为了打造我自己的DFT模拟工作流，在本科毕业的那个夏天，我开始为几个DFT程序的常见任务编写成通用脚本来生成输入文件，该项目被命名为`emuhelper`。但是随着项目的推进我发现它有非常大的扩展空间，于是在研一开学后我利用课余时间对其进行开发和完善，并将其重新命名为`Pymatflow`。

由于我爱好开源项目，所以对商业的VASP我并没有在Pymatflow中进行支持，不过考虑到以后可能会使用VASP工作，便另外启动了一个叫做vaspstudio的项目来针对VASP构建我的工作流。

## 使用Pymatflow命令行的通式

```
xxx-xxx.py xxx
```





## Pymatflow之qe计算器

Pymatflow的qe计算助手开发比较完善，目前支持静态scf、nscf、能带、态密度、结构优化、过渡态搜索、pp.x各类计算等。

## Pymatflow之cp2k计算器

Pymatflow开发之初的cp2k输入文件生成核心非常小，但是逐渐作者发现cp2k的输入文件的分块、分级的结构可以用类的包含关系来反应，于是改写了输入文件生成核心的代码，这样的好处是使用cp2k的生成器会非常简单有逻辑，但是坏处是类非常多，导致内存占用较多。

为了平衡内存占用过多的缺点，在使用时，对于有需要处理大量不同结构的需求时，建议不要以为每一个结构定义一个`pymatflow.cp2k`对象的方式来批处理，而是仅定义一个`cp2k`类，然后通过`get_xyz`函数来迭代获取结构进行计算。

## Pymatflow之SIESTA计算器

Pymatflow对SIESTA的支持目前仅停留于siesta上，对transiesta、tbtrans及对应的输运计算模块的支持还在开发之中。

## Pymatflow之Abinit计算器

Abinit输入文件自带的多数据(multidata)模式没有得到Pymatflow的支持，因为作者对Abinit的使用并不太多对应的Pymatflow方面的支持还不是特别完善。

## Pymatflow之elk计算器

开发中

## Pymatflow之dalton计算器

开发中

## Pymatflow之orca计算器

开发中

## Pymatflow之lammps计算器

开发中

## Pymatflow之gromacs计算器

开发中

## 特定类型计算流程控制核心`pymatflow.flow`

`pymatflow.flow`中着力对一些特定类型的计算进行开发，比如势能面搜索等。

* `flow-pes.py`静态势能面搜索。

## 实用结构相关工具

结构转换工具:

* `cif-to-xyz-modified.py`依赖于第三方工具`cif2cell`
* `cif-to-pdb.py`
* `pdb-to-cif.py`
* `xyz-modified-to-cif.py`依赖于第三方工具`vasp2cif` [p.s. `xyz-modified-to-cif.py`位于vaspstudio项目中]
* `xyz-modified-to-cif.py`
* `xyz-modified-to-crystal.py`

结构信息分析工具:

* `xyzinfo`利用seekpath库来获取k点路径，以及晶体结构信息。

超胞构建工具:

* `xyz-build-supercell.py`

其它:

* `qe-fix-atom.py`且qe生成固定原子的xyz文件

## 服务器交互实用工具

* `thq.py`
* `thpull.py`
* `thcancel.py`
* `thcmd.py`
* `threport.py`
* `sz-cmd.py`
* `sz-del.py`
* `sz-q.py`
* `sz-pull.py`

## `CP2K`助手使用Wiki

除了命令脚本`cp2k-*.py`的使用，也可以利用助手来实现批量计算或者说高通量计算.  下面是一个简单示例。

```
from pymatflow.cp2k.opt import opt_run
task = opt_run()
task.get_xyz("xxx.xyz")
task.set_geo_opt()
task.geo_opt(directory="xxx")
```



## 各类工具使用Wiki

**`pes-flow-static.py`**

`pes-flow.py`的工作方式是，读取xyz文件，指定最后n个原子作为一个集合进行运动，对每个运动的镜像进行静态scf计算获取能量，然后作图。

注意`--xrange`指定的是运动原子的x坐标变化的**相对**变化范围，`--yrange`同理。该工具还会自动生成用于生成xyz轨道文件的bash脚本`get_traj.sh`(运行较慢)。xyz格式的轨道文件可以通过`xcrystal --xyz trajectory.xyz`可视化。

**`pes-flow-relax.py`**

目前`pes-flow-static.py`工具通过静态的scf计算来获取每个坐标点对应的能量。实际上还可以考虑在每个坐标点进行relax计算，relax的时候保持x和y不变，只允许z方向移动，这样得到的势能面对于判断吸附位点更有意义，因为在每个xy点上对z进行了优化，得到了对应xy上的能量最小的z。这个时候可以保证势能面上xy点对应的是最小能量状态。

那么该如何实用`pes-flow-relax.py`以及`pes-flow-static.py`呢? 很简单，我们需要准备一个包含体系结构的xyz文件，文件中包含有计算中所涉及到的所有原子，但是要注意，用于扫描势能面的要移动的原子或原子团簇需要放在xyz文件的最后，然后通过`--last-n-move`传递给`pes-flow-relax|static.py`脚本。然后还是老规矩xyz文件第二行要定义晶胞参数来工`Pymatflow`使用。

当然为了实现定制化的势能面计算，我们还需要给程序额外的参数。与几何优化relax有关的参数这里就不作阐述，我们来看看与势能面有关的参数。首先是前面有讲过的`--last-n-move`参数用于指定进行移动扫描的原子或团簇。然后就是设定扫描空间范围的`--xrange`、`--yrange`、`--zshift`、参数。xyz文件中的移动的原子的初始位置将作为定义扫描空间范围的参考位置，xrange确定了移动原子的x的相对变化范围和步长，比如`--xrange -2 2 0.5`表示所有移动的原子的坐标将在原始值+(-2)到原始值+(2)这个范围内进行改变，改变 的步长是`0.5 angstrom`。



**`qe-get-matdyn-qpoints-from-bands-calc.py`**

通过分析qe能带结构的输出来获取k点用与设置qe的matdyn.x计算中使用的q点。非常不错的方案。

注意再作能带计算的时候，如果pw.x的band类型计算的高对称点的设定是以crystal_b的方式进行，那么要注意在进行pw.x的band类型计算后，再用band.x计算后在band.x输出文件中，你会发现高对称点的坐标不一定一样了，因为band.x输出的k点坐标不是crystal_b类型，你会发现其类型被转换成了tpiba_b，单位是2pi/a。所以你会发现如果是四方晶胞，且a=8.82，c=6.22，时，crystal_b类型的A(0.5, 0.5, 0.5)在bands.x的输出文件中会对应于点(0.5, 0.5, 0.709)，0.5与0.709的比值恰好等于c与a的比值。因为crystal_b类型的kx, ky, kz的单位是对应的2pi/a, 2pi/b, 2pi/c，那么当转换为tpiba_b时，都以2pi/a为单位，则新的坐标肯定只有kx不会变，ky，kz均会变化，只不过这里四方a=b，那么ky坐标变化比例为1，也就不变了。但是c != a，于是kz变了。

然后此脚本会根据能带结构计算中的k点来给出matdyn.x计算用得到的q点，从Quantum Espresso官网给出的Summer school on [Advanced Materials and Molecular Modelling with Quantum ESPRESSO](http://qe2019.ijs.si/), Ljubljana, Slovenia, September 15-20, 2019教程的[Github仓库](https://gitlab.com/QEF/material-for-ljubljana-qe-summer-school)的Day-3/example1b的matdyn.Si.in文件中的q点设置我们可以初步判断matdyn.x的q点应该也是tpiba_b类型(我还不能完全确定!!!!!!!)，那么本脚本在给出q点的时候主要就要统一使用bands.x输出中的tpiba_b类型的k点，就算pw.x的band类型的计算输入文件中以crystal_b类型给出高对称点，我们也只需要利用到其中对应的#label注释，然后对应高对称点转化为tpiba_b后，也是在bands.x的输出中，所以在给出matdyn.x的qpint的时候我们以bands.x中的输出为准。



## 关于各程序中K点的设置

此前已经将由seekpath生成高对称点以及K点路径集成到了能带计算和Phonopy声子谱计算中，但是我发现seekpath给出的k点对应的cell参数与输入的参数相比可能是变化了的，因此可能不值得信赖。

我打算将集成的seekpath生成k点部分改写为手动读取K点文件，然后以后计算时单独提供K点文件，如果要用seekpath，也可以用脚本生成可用的K点文件，然后pymatflow读取来进行计算。目前QE已经实现了手动读取，但是还没有移除掉内部的seekpath部分代码。

对于`Quantum ESPRESSO`有一个技巧就是:

在pw.x的bands类型计算的crystal_b类型的高对称点设置时，如果某个点后的 整数个连接点被设置为0，那么其效果是最后可以从bands.x的输出中查看到两个高对称点的xcoord是一样的。这是合理的，因为两个点之间的链接数为0，那么中间距离肯定就是0了。这样有一个好处就是，我们在设置$\Gamma–X–M–\Gamma–Z–R–A–Z|X–R|M–A$类型的布里渊区K-Path的时候，我们可以设置高对称点$Z|X$、$R|M$中链接的$Z$和$X$以及$R$和$M$的点数为0，那么最后相对应的$Z$与$X$以及$R$与$M$点的xcoord就会是一样的，我们在作能带图的时候就可以用$Z|X$、$R|M$这样的符号来表示，很方便。并且我的qe.post.band已经能够自动处理这种情况，从pw.x的bands计算输入文件中读取高对称点符号，从bands.x输出中获得对应高对称点的xcoord，下一步就是很不错的一步，对以后的高对称点label和xcoord进行refine处理，如果发现某两个相邻点的xcoord一样，那么就将两个合并，丢掉一个xcoord，并且将两个label进行合并为"xxx|xxx"的形式。

在SIESTA的Bandsline模式中，我也想这样实现，但是不幸的是SIESTA中设置0个点来连接相邻高对称点时，计算输出的SystemLabel.bands文件中可以看到，程序直接从0个连接点处断开，后面的点就没输出了，意味着按照QE那种设置0个点链接断开点的技巧$\Gamma–X–M–\Gamma–Z–R–A–Z|X–R|M–A$这样一条路径，最后只会计算输出$\Gamma–X–M–\Gamma–Z–R–A–Z$，SIESTA读取到链接Z到X的点为0个时，就会舍弃后面的！！！！这个方法是行不通的，我又尝试了设置为1个点来链接，但是SystemLabel.bands文件末查到的对应两个高对称点的xcoord就不一样了，这不适合脚本向qe后处理那样自动处理。

**`post-xxx-phonoy.py`**的使用

`--qpath`参数，与能带计算中的``--kpath`参数格式一样，尽管phonoy不需要单独设置每个高对称q点之间链接点的数量，但是为了使用习惯的一致性，`--qpath`参数仍然需要`--qpath '0.0 0.0 0.0 GAMMA 10' ' 0.5 0.0 0.0 X |'`这样来设置。其中每个q点的最后一个值如果是"|"，那么表示其未与后面的q链接。如果不是"|"就表示是链接的。由于该值不会被脚本传递给Phonopy使用，因此除了"|"以外的任何合法字符都会被认为是表明是链接状态。

如果是在其它脚本中`--kpath`的每个k点的最后一个参数应该要么是"|"要么是一个整数，表示该点与后面的点断开或者以该整数个点进行链接。

### K点的单位类型

就像实空间坐标有直角坐标(Cartesian)和晶体坐标(Crystal又叫做分数坐标Fractional)一样，倒易空间的K点也有直角坐标和晶体坐标(有时也叫作reciprocal坐标)。QE的crystal_b、CP2K的B_VECTOR、SIESTA的ReciprocalLatticeVectors类型的BandLinesScale、Abinit的kpt、kptbounds，以及vasp设置rec类型的K点都是晶体坐标或者说reciproca primitive vector形式。这样的kx、ky、kz的单位长度其实就分别是b1、b2、b3。

QE的tpiba_b类型的K点就是直角坐标，kx、ky、kz的单位长度都是2pi/a。

晶体坐标类型的K点大多都是有理数，使用起来比较方便。

## 关于为什么，我的几何结构核心选取xyz结构:
* 由于我的计划会在pymatflow中实现许多对结构的操作，而我目前的几何只是只能支持我在Cartesian这样的直角坐标上进行操作。
* 另外，没有必要支持所有类型的结构文件的直接结构操作，我们以xyz文件(modified version)为核心进行操作，然后提供工具
# 将其转换为其它类型的结构文件就可以了。

## TODO

* VASP计算数据的提取，单独作PDOS、BAND图。
* QE考虑SOC，以及对应PDOS、BAND分析工具对SOC的支持。
* QE的DFPT计算声子谱的后处理，作图考虑"|".
* phonopy声子谱数据的提取，额外作图。
* Abinit、提交脚本的支持。
* Abinit的Phonopy和能带测试。
* SIESTA的Phonopy和能带测试。
  * 目前我在服务器上安装的SIESTA计算退出似乎有异常，在计算完一个disp后无法正常继续计算下一个disp，且仍然占用计算资源，我发现作业脚本应该没问题，或许是SIESTA安装有问题。
* CP2K能带后处理。
* 一个单独的SeekPath模块用于为各类计算器准备K点，目前已经有了脚本，但是我想设计出模块，方便以后高通量计算使用。以前其实每个计算器都耦合了，但是我把它们移除了，打算单独抽出来，统一一个kpath的数据格式，来供各计算器使用。
* 各类计算时间信息的分析。
* 为各类计算编写保障性代码。
* QE的epsilon.x以及turbo_davidson.x、turbo_lanczos.x、turbo_spectrum.x已经写好了，但是需要测试。
* 构建一个qe类，将其它如pwscf、ph.x、epsilon.x、projwfc.x单独抽象为类后包含到qe类中。
* 将auto默认设置为在服务器上运行、
* 将xyz类更进一步分离出来，与计算器的结构部分代码进一步解离。
* 编写测试代码，可以改造一些base_xyz，让他除了能够从xyz文件获取结构，也能直接从一个python数据结构获取结构。

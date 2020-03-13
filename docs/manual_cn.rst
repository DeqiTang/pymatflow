.. _header-n0:

Pymatflow用户手册
=================

``Pymatflow``\ 是一个第一性原理模拟的工作流自动化管理软件，有两种工作模式。其一是终端命令行模式，目的是减少用户进行DFT模拟所需要的输入文件准备和结果后处理时间，让用户在服务器集群上通过一行命令就可以提交任务和进行结果的提取。其二是使用\ ``Pymatflow``\ 提供的\ ``API``\ 进行高通量计算，这部分的使用需要用户具有一定的编程能力。目前支持的计算程序有\ ``Abinit``\ 、\ ``Quantum ESPRESSO``\ 、\ ``CP2K``\ 、\ ``SIESTA``\ 、\ ``VASP``

由于个人精力有限，本手册内容并不完善，目前主要对命令行工具\ ``matflow``\ 、\ ``postflow``\ 和\ ``structflow``\ 的使用方式进行了初步介绍，以方便用户试用。

.. _header-n4:

安装Pymatflow
-------------

.. code:: shell

   ~$ pip install pymatflow

目前推荐的方式是在服务器中使用\ ``pip``\ 安装已发布的pymatflow。

如果自己不是系统管理员，没有root权限:

.. code:: shell

   ~$ pip install pymatflow --user

并保证\ ``~/.local/bin``\ 在环境变量\ ``PATH``\ 之中。

``pymatflow``\ 支持本地与远程服务器交互操作，但是需要做一定配置，后续做补充。

.. _header-n11:

使用须知
~~~~~~~~

首先用户使用本软件需要对结果自行负责，Pymatflow不保证结果的正确性。需要用户具有一定的背景知识，能够鉴别计算任务的合理性。

1. 结构文件

可以通过\ ``--cif``\ 、\ ``--xsd``\ 、\ ``--xsf``\ 或者\ ``--xyz``\ 三个参数之一传递结构文件路径给\ ``matflow``\ 或者\ ``postflow``\ 。用户需要注意的是这四个参数一次只能设置一个。如果是进行过渡态计算，则通过\ ``--images``\ 参数输入结构。

1. 赝势文件

``matflow``\ 提供了\ ``--pot``\ 选项供用户选择赝势文件，有三种设置方式。第一种是指定一个文件夹，该文件夹中需要包含所有的赝势文件，并且建议用户在该文件夹中仅存放赝势文件，因为目前\ ``matflow``\ 会将所有该文件夹下的内容都拷贝到当前路径。其二是默认形式，即不需要指定\ ``--pot``\ 选项，但用户需要保证所有赝势文件都在当前路径下。第三种方式是\ ``--pot auto``\ ，此时\ ``matflow``\ 会自动会你准备好赝势文件。

对于第三种方式，我们需要对服务器进行一下配置，将赝势库都准备在用户家目录下面。各个程序的赝势存放位置分别是:``~/.pot-vasp``\ 、\ ``~/.pot-qe``\ 、\ ``~/.pot-siesta``\ 、\ ``~/.pot-abinit``\ 。需要注意的是目前对\ ``CP2K``\ 的支持不需要指定赝势，默认只会使用\ ``DZVP-MOLOPT-SR-GTH``\ 基组和\ ``GTH-PBE``\ 类型的赝势。

1. 高对称点路径

在对晶体材料的体相进行研究时，常会用到倒空间这个强有力的工具。在使用\ ``Pymatflow``\ 的过程中，所有涉及倒空间高对称点路径的任务都以一个规范的形式进行管理。你可以通过\ ``--kpath-manual``\ 或者\ ``--kpath-file``\ 向\ ``matflow``\ 传入高对称点路径。

使用\ ``--kpath-manual``\ 时，需要在其后按照一定格式手动输入路径:

.. code:: shell

   --kpath-manual '0.000000 0.000000 0.000000 GAMMA 15' '0.500000 0.000000 0.000000 X |' '0.000000 0.500000 0.000000 Y 15' '0.000000 0.000000 0.000000 GAMMA 15' '0.000000 0.000000 0.500000 Z |' '0.500000 0.500000 0.500000 R 15' '0.000000 0.000000 0.000000 GAMMA 15' '0.000000 0.500000 0.500000 T |' '0.500000 0.000000 0.500000 U 15' '0.000000 0.000000 0.000000 GAMMA 15' '0.500000 0.500000 0.000000 V 15'

上面的参数给出的高对称点路径是::math:`\Gamma-X|Y-\Gamma-Z|R-\Gamma-T|U-\Gamma-V`\ 。

Notes:

-  如果一个K点的最后一个符号
   为一个整数，那么它将通过该整数个点连接到下一个K点。

-  如果一个K点的最后一个符号为\ ``|``\ ，那么它与后面的K点处于断开状态。

-  ``kx, ky, kz``\ 是倒空间K点坐标的分数坐标。

此外，用户可以将高对称点路径信息准备在文件中，然后将文件路径通过\ ``--kpath-file``\ 传递给\ ``Pymatflow``\ 。文件格式如下:

.. code:: text

   11
   0.000000 0.000000 0.000000 #GAMMA 15
   0.500000 0.000000 0.000000 #X |
   0.000000 0.500000 0.000000 #Y 15
   0.000000 0.000000 0.000000 #GAMMA 15
   0.000000 0.000000 0.500000 #Z |
   0.500000 0.500000 0.500000 #R 15
   0.000000 0.000000 0.000000 #GAMMA 15
   0.000000 0.500000 0.500000 #T |
   0.500000 0.000000 0.500000 #U 15
   0.000000 0.000000 0.000000 #GAMMA 15
   0.500000 0.500000 0.000000 #V 15

文件的第一行制定了文件后面给出的行数。

每一行对应于\ ``kpath``\ 数据结构中的一个元素:

-  行中前三个元素对应于K点的x、y和z坐标

-  行中第四个元素对应于K点的符号(需要是大写)

-  行中第五个元素对应于K点与其后K点的连接关系，意义和上述\ ``kpath``\ 数据结构中一致

用户往往需要通过\ ``--kpath-file``\ 参数向\ ``Pymatflow``\ 的命令行工具(如\ ``matflow``\ 、\ ``postflow``\ 、\ ``structflow``)指定包含有高对称K点的文件路径。要注意这里同样坐标必须是分数坐标。注意\ ``--kpath-manual``\ 的优先级高于\ ``--kpath-file``\ 。

1. 环境变量

用户需要为\ ``Pymatflow``\ 指定所用程序的环境变量，比如将下列内容加入到你的\ ``~/.bashrc``\ 中:

.. code:: shell

   export PMF_VASP_STD=/path/to/your/vasp_std
   export PMF_VASP_NCL=/path/to/your/vasp_ncl
   export PMF_VASP_GAM=/path/to/your/vasp_gam

然后通过下面的命令使其立即生效:

.. code:: shell

   ~$ source ~/.bashrc

.. _header-n56:

``matflow``\ 简介
-----------------

``matflow``\ 提供了计算任务的自动化生成和提交功能。使用方式归结为一句话就是"主命令+子命令+可选参数"。主命令当然就是\ ``matflow``\ 了，至于子命令，目前一共有5个，即:

-  ``abinit``

-  ``qe``

-  ``cp2k``

-  ``siesta``

-  ``vasp``

使用时也就代表了命令使用到的计算器。用户可以通过\ ``matflow -h``\ 查看使用帮助。如果想要具体了解某一个计算器提供的可选参数和使用帮助需要加上子命令选项，比如命令\ ``matflow vasp -h``\ 可以查看\ ``vasp``\ 对应的帮助文档。下面我们通过一个例子来看如何在服务器上使用\ ``matflow``\ 来辅助第一性原理计算。

另外为了减少输入，用户可以使用\ ``mflow``\ 代替\ ``matflow``\ 。

.. _header-n71:

``matflow``\ 命令工具使用示例
-----------------------------

这里以一个简单的优化\ ``LiH``\ 立方晶胞参数的例子，展示如何使用Pymatflow来加速计算的准备工作。

首先，准备好\ ``LiH``\ 的晶体结构文件\ ``lih.cif``\ 。然后在对应路径运行一下命令:

.. code:: shell

   ~$ matflow vasp -r 2  --cif lih.cif --encut 300 --ibrion 2 --isif 2 --kpoints-mp 3 3 3 0 0 0 -d lih-cubic

即可在指定的\ ``lih-cubic``\ 文件夹中自动生成输入文件，并提交任务到服务器中，要注意的是这里默认提交的是PBS类型的任务。目前支持的服务器只有吕梁天河二号和PBS集群，或者以单机模式运行。具体使用参见\ ``--server``\ 参数的帮助文档。至于结果的提取由\ ``postflow``\ 工具提供。目前\ ``matflow``\ 对vasp支持的计算类型包括:
静态计算(scf、nscf、bands)、结构优化、立方晶胞参数优化、六方晶胞参数优化、四方晶胞参数优化、结合VTST进行过渡态计算、VASP自带声子谱计算、结合Phonopy进行声子谱计算。

.. _header-n76:

``postflow``\ 简介
------------------

``postflow``\ 对Pymatflow的后处理进行了部分封装，目前还未完全成熟。可能有部分选项无法使用。这里展示一下通常的使用流程。比如对上述的\ ``LiH``\ 晶胞参数优化任务的后处理，只需要运行命令:

.. code:: shell

   ~$ postflow vasp -r 2 -d lih-cubic

就可以进行结果的提取，生成的文件统一在\ ``lih-cubic/post-processing``\ 目录下。

目前Pymatflow的后处理功能主要还是由API提供，但是往\ ``postflow``\ 移植的工作正在进行中，后续会完善。

.. _header-n81:

``structflow``\ 简介
--------------------

``structflow``\ 主要提供了常见的结构文件转换，以及进行固定原子结构优化等功能，其中结构文件转换得益于ASE项目的\ ``ase.io``\ 的帮助。这部分内容会在后续的文章中进行介绍

.. _header-n83:

使用Pymatflow命令行的通式
-------------------------

.. _header-n84:

``matflow``\ 之VASP计算器
-------------------------

.. _header-n85:

VASP静态计算
~~~~~~~~~~~~

.. code:: shell

   ~$ matflow vasp -r 0 --cif YOUR_STRUCTURE.cif -d DIRECTORY --encut VALUE --ediff VALUE --kpath-file KPATH_FILE_PATH --kpoints-mp VALUE

参数\ ``-r 0``\ 表示静态计算，其行为是scf、nscf、bands一起进行计算。其中需要注意\ ``-d``\ 参数不能指定为当前路径，否者命令会抛出警告并退出。后处理命令如下:

.. code:: shell

   ~$ postflow vasp -r 0 -d DIRECOTRY --kpath-file KPATH_FILE_PATH --kpoints-mp VALUE

.. _header-n89:

VASP结构优化
~~~~~~~~~~~~

.. code:: shell

   ~$ matflow vasp -r 1 --cif YOUR_STRUCTURE.cif -d DIRECTORY --encut VALUE --ediff VALUE --ediffg VALUE    --ibirion VALUE --isif VALUE --kpoints-mp VALUE

参数\ ``-r 1``\ 表示进行结构优化。其中需要注意\ ``-d``\ 参数不能指定为当前路径，否者命令会抛出警告并退出。后处理命令如下:

.. _header-n92:

VASP晶胞参数优化
~~~~~~~~~~~~~~~~

.. code:: shell

   ~$ matflow vasp -r [2|3|4] --cif YOUR_STRUCTURE.cif -d DIRECTORY --encut VALUE --ediff VALUE --kpoints-mp VALUE

参数\ ``-r 2``\ 是进行cubic晶胞的cell参数的优化，\ ``-r 3``\ 进行hexagonal晶胞的cell参数的优化，\ ``-r 4``\ 进行tetragonal晶胞的cell参数的优化。其中需要注意\ ``-d``\ 参数不能指定为当前路径，否者命令会抛出警告并退出。后处理命令如下:

.. code:: shell

   ~$ postflow vasp -r [2|3|4] -d DIRECOTRY

.. _header-n96:

VASP过渡态计算(VTST)
~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   ~$ matflow vasp -r 5 --images INITIAL.cif FINAL.cif -d DIRECTORY --encut VALUE --ediff VALUE --kpoints-mp VALUE

参数\ ``-r 5``\ 代表进行过渡态计算。其中\ ``--images``\ 用于指定初始和终态结构。

.. _header-n99:

VASP声子谱计算(内置IBRION=[5, 6, 7, 8])
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   ~$ matflow vasp -r 6 --cif YOUR_STRUCTURE.cif -d DIRECTORY --encut VALUE --ediff VALUE --kpoints-mp VALUE --ibirion [5|6|7|8] --supercell-n VALUE

参数\ ``-r  6``\ 是利用vasp内置的算法进行声子谱计算，需要设置\ ``--ibrion``\ 为\ ``5|6|7|8``\ 之一。其中需要注意\ ``-d``\ 参数不能指定为当前路径，否者命令会抛出警告并退出。后处理命令如下:

.. code:: .

   ~$ postflow vasp -r 6 --cif YOUR_STRUCTURE.cif -d DIRECTORY --kpath-file KPATH_FILE_PATH

.. _header-n103:

VASP声子谱计算(Phonopy)
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   ~$ matflow vasp -r 7 --cif YOUR_STRUCTURE.cif -d DIRECTORY --encut VALUE --ediff VALUE --kpoints-mp VALUE --supercell-n VALUE

.. _header-n105:

特定类型计算
------------

``Pymatflow``\ 提供部分特定类型计算任务，文档暂无。

.. _header-n107:

服务器交互实用工具
------------------

-  ``thq.py``

-  ``thpull.py``

-  ``thcancel.py``

-  ``thcmd.py``

-  ``threport.py``

-  ``sz-cmd.py``

-  ``sz-del.py``

-  ``sz-q.py``

-  ``sz-pull.py``

这些工具是本地与远程服务器交互的工具，具体使用需要参数配置，后续补充。

.. _header-n128:

API
---

文档暂无

.. _header-n130:

未来
----

首先欢迎有兴趣的朋友可以一起参与开发。也欢迎用户对程序的使用提出建议，或者提交功能需求。

-  项目地址: https://gitlab.com/DeqiTang/pymatflow

-  文档地址: http://pymatflow.readthedocs.org/

.. _header-n137:

问题反馈
--------

本文档所有权归属\ ``Pymatflow``\ 项目。对项目有任何问题反馈，请发送邮件至

📫 pymatflow@163.com

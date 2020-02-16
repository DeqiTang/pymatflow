Phonopy+Quantum ESPRESSO计算固体声子谱(Pymatflow简化版)
=======================================================
本教程将简单介绍如何使用Pymatflow简单地实现利用Phonopy+Quantum ESPRESSO实现对固体声子谱的计算。

Phonopy通过体系的原子间力常数来计算声子谱，而力常数的计算需要借助于第一性原理程序，目前非常多的DFT程序都得到了支持，比如VASP、QE、Abinit、CP2K等。

Pymatflow工具将使用Phonopy的过程进行了封装，能够方便地利用Phonopy结合其它DFT程序，比如QE、CP2K、Abinit来快速进行声子谱的计算。

使用流程
-------

* 准备好贋势文件和xyz结构文件，注意第二行需要按照一定方式注明cell参数.
* 运行命令提交计算

案例演示
-------

最近待在家里，切身体验柴米油盐酱醋茶的滋味。为此这次的研究对象就选取每天都能尝到的NaCl好了。这里的NaCl为常见的面心立方结构，空间群为225，结构如图所示:

.. image:: /_images/nacl.png

准备xyz类型的结构文件(nacl.xyz)如下::

    8
    cell:   5.6916936   0.0000000   0.0000000  |   0.0000000   5.6916936   0.0000000  |   0.0000000   0.0000000   5.6916936 
    Na      0.0000000   0.0000000   0.0000000
    Na      0.0000000   2.8458468   2.8458468
    Na      2.8458468   0.0000000   2.8458468
    Na      2.8458468   2.8458468   0.0000000
    Cl      0.0000000   0.0000000   2.8458468
    Cl      0.0000000   2.8458468   0.0000000
    Cl      2.8458468   0.0000000   0.0000000
    Cl      2.8458468   2.8458468   2.8458468

然后在服务器上运行以下命令:

``qe-phonopy.py -f nacl.xyz --ecutwfc 60 --kpoints-mp 3 3 3 0 0 0 --supercell-n 1 1 1``

然后你就可以静静等待任务运行结束。然后运行后处理程序。下面解释一下参数的意义: 

* -f指定结构文件
* --kpoints-mp指定Monkhorst-Pack的K点方案
* --supercell-n是通过pymatflow传递给Phonopy的超胞参数，对应于Phonopy的--dim参数

下面便是结果:

声子谱:

.. image:: /Figure_3.png


这次的内容就到这里了，后面我还会继续更新Pymatflow的使用案例，并继续整理文档，欢迎使用和讨论，有意参与开发者或者有任何建议者可以联系我。

更多文章
----------

`Band Structure <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&mid=2247484108&idx=1&sn=8433830398824f147bd98b46893803c6&chksm=fdb9b539cace3c2f1e4b673f9d5f5cc039dbd8e382f874e77a935515b86b1fb2c12baddec5ae&token=1365138185&lang=zh_CN#rd>`_

`Phonopy+QE <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&mid=2247484116&idx=1&sn=a3415739cc964015938009c1d8656988&chksm=fdb9b521cace3c3783d1c747a3a0253ba9716277db906637cbaf39242ca4bcb472d39086aeee&token=1365138185&lang=zh_CN#rd>`_



微信公众号
-----------

.. image:: /_images/qrcode_for_shengcaiyoudao_1.jpg

欢迎关注微信公众号"生材有道", 可以订阅更多资讯。

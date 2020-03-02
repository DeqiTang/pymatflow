中文版教程
=========
目前关于Pymatflow的模块接口使用的教程暂未完成。用户现在如何使用已经编写好的命令行工具
来用于日常的研究当中，了解使用方式的最佳方式是使用案例(Examples)。

命令行脚本使用示例(部分为公众号文章链接)
----------------
.. toctree::
   :maxdepth: 2

   `使用Phonopy+CP2K计算固体声子谱(Pymatflow篇) <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&tempkey=MTA1MF9keDgwYjk5STdpRlZFVGVyZnlSRGlZTXI4LW5IdUN6RTd6a3RmNm9tOHVMc0ZvVmRTQnUwQWJDeENUZWNMOEx1NW1UVFl5aW1EM1owZGpvVDRZWFNTN0FmTWtKRU1mTWQ4UEc0cnJlT202VlJNUU5JUEJzc2VXSVc1NHo4cFBjSGxiaUFVQnJaR2YzYjVOeDNJYm00WWdXdGhTZE92WEs3M0phS3N3fn4%3D&chksm=7db9b52f4ace3c3908f0b9f37b78bf5682ab04a1d432f7ac7dd2a73d5d3e977ce2935d5987ba#rd>`_
   examples/phonopy_qe.rst
   examples/bandstructure_qe.rst


注意事项
--------
在使用Pymatflow命令行工具自动进行第一性原理模拟输入文件的生成、任务提交和结果后处理的过程中，如果涉及
到能带、声子谱计算等涉及到倒空间高对称点路径的模拟，用户可统一通过``--kpath``或者``--kpath-file``提供
高对称点路径信息(推荐使用后者)。``--kpath-file``用于指定文件包含有高对称点路径，文件格式如下::

    5
    0.000000 0.000000 0.000000 #GAMMA 15
    0.500000 0.000000 0.000000 #X |
    0.000000 0.500000 0.000000 #Y 15
    0.000000 0.000000 0.000000 #GAMMA 15
    0.000000 0.000000 0.500000 #Z


其中，第一行是一个整数指定一共有多少个高对称点。后面的行用于指定高对称点的坐标信息符号信息(必须大写)以及
连接情况。如果某一个高对称点与其后的高对称点是相连的，那么每一行的最后应当是一个整数用于指定
两个高对称点连线中间的k点个数。如果一个高对称点与其后的高对称点断开，则在该行末尾用``|``指定。

如图中第一个``GAMMA``点通过15个点连接到了``X``点，然后``X``点与其后的``Y``点处于断开状态。

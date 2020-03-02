.. Pymatflow documentation master file, created by
   sphinx-quickstart on Sun Feb 16 14:35:45 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python Materials WorkFlow
==========================
The Python Materials WorkFlow (Pymatflow) is an workflow simplifier for reserach
on Materials Science by means of Ab initio simulation.
Pymatflow provides module interface to several popular open source Ab initio programs
along with many useful command line tools which can help setting up, post-processing
the simulations.
A basic usage of Pymatflow (via module interaface) is like this::

  >>> from pymatflow.qe.import static_run
  >>> task = static_run()
  >>> task.get_xyz("CO2.xyz")
  >>> params = {
  >>>   "encutwfc":60,
  >>>   "degauss": 0.001
  >>>   }
  >>> task.set_params(params)
  >>> task.scf(runopt="genrun", auto=0)

The above code will read the xyz structure file ``CO2.xyz`` and do a basic scf
calculation on it. Note, to make sure ``task.get_xyz()`` run normally we must specify cell
parameters in the second line of ``CO2.xyz`` manually, in format like this::

  3
  cell: 10 0 0 | 0 10 0 | 0 0 10
  C  0.0000000 0.000000 0.000000
  O -1.1600000 0.000000 0.000000
  O  1.1500000 0.000000 0.000000

Tutorials
-----------------------
.. toctree::
   :maxdepth: 1

   tutorials/tutorials_en.rst
   tutorials/tutorials_cn.rst

For Chinese User
----------------
Pymatflow官方微信公众号"生材有道"目前已经开始试运营。欢迎订阅查看更多关于Pymatflow的最新
信息和使用教程。

.. image:: /_images/qrcode_for_shengcaiyoudao_1.jpg

往期公众号文章
`使用Phonopy+CP2K计算固体声子谱(Pymatflow篇) <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&tempkey=MTA1MF9keDgwYjk5STdpRlZFVGVyZnlSRGlZTXI4LW5IdUN6RTd6a3RmNm9tOHVMc0ZvVmRTQnUwQWJDeENUZWNMOEx1NW1UVFl5aW1EM1owZGpvVDRZWFNTN0FmTWtKRU1mTWQ4UEc0cnJlT202VlJNUU5JUEJzc2VXSVc1NHo4cFBjSGxiaUFVQnJaR2YzYjVOeDNJYm00WWdXdGhTZE92WEs3M0phS3N3fn4%3D&chksm=7db9b52f4ace3c3908f0b9f37b78bf5682ab04a1d432f7ac7dd2a73d5d3e977ce2935d5987ba#rd>`_

`使用Pymatflow辅助Quantum ESPRESSO计算能带 <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&mid=2247484108&idx=1&sn=8433830398824f147bd98b46893803c6&chksm=fdb9b539cace3c2f1e4b673f9d5f5cc039dbd8e382f874e77a935515b86b1fb2c12baddec5ae&token=1365138185&lang=zh_CN#rd>`_

`Phonopy+Quantum ESPRESSO计算固体声子谱(Pymatflow简化版) <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&mid=2247484116&idx=1&sn=a3415739cc964015938009c1d8656988&chksm=fdb9b521cace3c3783d1c747a3a0253ba9716277db906637cbaf39242ca4bcb472d39086aeee&token=1365138185&lang=zh_CN#rd>`_



Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. Pymatflow documentation master file, created by
   sphinx-quickstart on Sun Feb 16 14:35:45 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pymatflow's documentation!
=====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules/modules.rst
   tutorials/tutorials.rst
   wiki/server.rst

About Pymatflow
----------------
see :ref:`about`

Basic usage of Pymatflow
------------------------
use of pymatflow follows code like this::

  from pymatflow.qe.import static_run
  task = static_run()
  task.get_xyz("xxx.xyz")
  params = {}
  params["encutwfc"] = 60
  params["degauss"] = 0.001
  task.set_params(params)
  task.scf(runopt="genrun", auto=0)

The above code will read the xyz structure file ``xxx.xyz`` and do a basic scf
calculation on it. Note, to make pymatflow run normally we must specify cell
parameters in the second line of ``xxx.xyz`` manually, in format like this::

  ```
  3
  cell: 10 0 0 | 0 10 0 | 0 0 10
  C 0.0000000 0.000000 0.000000
  O -1.160000 0.000000 0.000000
  O 1.1500000 0.000000 0.000000
  ```
Tutorials
-----------------------
see :ref:`tutorials/tutorials`

More Articles
--------------

`Band Structure <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&mid=2247484108&idx=1&sn=8433830398824f147bd98b46893803c6&chksm=fdb9b539cace3c2f1e4b673f9d5f5cc039dbd8e382f874e77a935515b86b1fb2c12baddec5ae&token=1365138185&lang=zh_CN#rd>`_

`Phonopy+QE <https://mp.weixin.qq.com/s?__biz=MzU4MjM5NDUyMg==&mid=2247484116&idx=1&sn=a3415739cc964015938009c1d8656988&chksm=fdb9b521cace3c3783d1c747a3a0253ba9716277db906637cbaf39242ca4bcb472d39086aeee&token=1365138185&lang=zh_CN#rd>`_

WeChat Official Account
-----------------------

.. image:: /_images/qrcode_for_shengcaiyoudao_1.jpg

Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

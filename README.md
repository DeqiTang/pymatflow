# Pymatflow

With scripts to provide preparing and/or post-processing for your Emulations(DFT, MD).

## Goals
* powerful and stable calculator supports for popular Ab initio programs
* fully featured results analysis
* user friendly structure builder
* automated workflow for high throughput calculation[including database handling]
  * automatical choice of appropriate parameters
  * automatical correction of error during running
* statistical analysis for properties, basic machine learning module
* Note: I am currently not willing to implement a general representation of calculations
    ```
    because I think a powerful input generation engine for different dft programs is more
    important at present. and an overall representation is just a trick to port task to
    specific calculators, which I am not preparing to make one now.[but probably there will
    be one in the future.]
    ```
* support for system visulization by third-party molecular visualizer like Avogadro, Vipster, paraview, vtk, vmd, xcrysden, vesta, blender
* support for third-party structure generation tools, like mbuild, atomsk, msmbuilder,
* using spgroup to get space group of crystals.
* add Sphix docs and host it on readthedoc.org



## About the th\*.py scripts
in order to use th\*.py scripts and the function of --auto in running scripts(like qe-relax.py), you have to prepare a server.conf file in ~/.emuhelper/server.conf, in format like this:
```
[server]
ip = x.x.x.x
user = xxx
serverdir = /xx/xx/xx # like /THFS/home/xxx/
password = xxx
```

## License
```
MIT License

Copyright (c) 2019 DeqiTang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## About

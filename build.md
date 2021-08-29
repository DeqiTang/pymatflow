# Build and Distribute Pymatflow

# References
* https://realpython.com/python-wheels/#different-types-of-wheels
* https://github.com/pypa/manylinux

# Multi Platform Support
Debian
```
docker pull quay.io/pypa/manylinux_2_24_x86_64
docker run -v /path/to/pymatflow:/root/pymatflow  -it quay.io/pypa/manylinux_2_24_x86_64 bash
# -----------------------
# inside docker container
# -----------------------
apt update
apt install libatlas-dev libblas-dev liblapack-dev  # needed for python3.10 to build scipy
# add pybind11 support, need to add pybind11 to PATH and CPLUS_INCLUDE_PATH
# or pybind11 will not be found by cmake
export PATH=$PATH:${HOME}/.local/bin
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${HOME}/.local/include

for py in cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
/opt/python/${py}/bin/pip install --user scikit-build pybind11[global]==2.7.1
done


for py in cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
# old build might destroy current build, so remove it
rm -rf /root/pymatflow/_skbuild
/opt/python/${py}/bin/pip wheel --use-feature=in-tree-build /root/pymatflow -w /root/pymatflow/dist/ 
done


for whl in /root/pymatflow/dist/pymatflow*-linux-*.whl
do
auditwheel repair ${whl} -w /root/pymatflow/dist
done
```
CentOS 7
```
docker pull quay.io/pypa/manylinux2014_x86_64
docker run -v /path/to/pymatflow:/root/pymatflow  -it quay.io/pypa/manylinux2014_x86_64 bash
# -----------------------
# inside docker container
# -----------------------

yum install atlas blas lapack  # needed for python3.10 to build scipy

# add pybind11 support, need to add pybind11 to PATH and CPLUS_INCLUDE_PATH
# or pybind11 will not be found by cmake
export PATH=$PATH:${HOME}/.local/bin
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${HOME}/.local/include

for py in cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
/opt/python/${py}/bin/pip install --user scikit-build pybind11[global]==2.7.1
done

for py in cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
# old build might destroy current build, so remove it
rm -rf /root/pymatflow/_skbuild
/opt/python/${py}/bin/pip wheel --use-feature=in-tree-build /root/pymatflow -w /root/pymatflow/dist/ 
done


for whl in /root/pymatflow/dist/pymatflow*-linux-*.whl
do
auditwheel repair ${whl} -w /root/pymatflow/dist
done
```

## pybind11 support
in server you cannot use apt and yum to install pybind11, you can use pip --user to install global version of pybind11. however, you need to add pybind11 to the env, by
```
pip3 install --user pybind11[global]==2.7.1
export PATH=$PATH:${HOME}/.local/bin
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${HOME}/.local/include
```
# or pybind11 will not be found by cmake
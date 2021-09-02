
apt update
apt install -y libatlas-dev libblas-dev liblapack-dev  # needed for python3.10 to build scipy
# need working pip3 command
apt install -y python3-pip

for py in cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
/opt/python/${py}/bin/pip3 install --user scikit-build cython
done

cd /root/pymatflow/
for py in cp36-cp36m cp37-cp37m cp38-cp38 cp39-cp39 cp310-cp310
do
# old build might destroy current build, so remove it
rm -rf _skbuild
/opt/python/${py}/bin/python3 setup.py bdist_wheel 
done

for whl in /root/pymatflow/dist/pymatflow*-linux_*.whl
do
auditwheel repair ${whl} -w /root/pymatflow/dist
done

#!/bin/bash

easy_install argparse

rpm -Uvh p7zip-9.20.1-1.el6.rf.x86_64.rpm

tar zxf cmake-2.8.12.2.tar.gz
cd cmake-2.8.12.2
./configure
make
make install
cd


tar zxf Cython-0.20.1.tar.gz
cd Cython-0.20.1
python setup.py build
python setup.py install
cd 

tar zxf numpy-1.7.2.tar.gz
cd numpy-1.7.2
python setup.py build
python setup.py install
cd

tar zxf scipy-0.13.3.tar.gz
cd scipy-0.13.3
python setup.py build
python setup.py install
cd


tar zxf matplotlib-1.3.1.tar.gz
cd matplotlib-1.3.1
python setup.py build
python setup.py install
cd

tar zxf pca_module-1.1.02.tar.gz
cd pca_module-1.1.02
sed -i".bak" 's/add_ext = True/add_ext = False/' setup.py
rm setup.py.bak
sed -i".bak" 's/old_numeric = True/old_numeric = False/' setup.py
rm setup.py.bak
python setup.py build
python setup.py install
cd

unzip python-graph-1.3.1.zip
cd python-graph
python setup.py build
python setup.py install
cd

unzip opencv-2.4.8.zip
cd opencv-2.4.8
cmake .
make
make install
cd

tar xvf CGAL-4.4.tar.gz
cd CGAL-4.4

# add a macro for GMP-4.3 compatibility

sed -i".bak" -n 'H;${x;s/\/\/ GMP before 5.0.*\n/\/\/ GMP-4.3.* is missing in mpn_sqr\n\n#ifndef mpn_sqr\n#define mpn_sqr(dest,a,n) mpn_mul_n(dest,a,a,n)\n#endif\n\n&/;p;}' include/CGAL/Mpzf.h
rm -f include/CGAL/Mpzf.h.bak

#echo "// GMP-4.3.* is missing in mpn_sqr"
#echo "#ifndef mpn_sqr"
#echo "#define mpn_sqr(dest,a,n) mpn_mul_n(dest,a,a,n)"
#echo "#endif"


cmake .
make
make install
cd ..

tar xvf guppy-0.1.10.tar.gz
cd guppy-0.1.10
python setup.py build
python setup.py install
cd




tar xf opaque_0.1.tar

cd opaque_0.1/src/modules
python setup.py build_ext --inplace

cd medialaxis
python setup.py build_ext --inplace
cd ..

cd alphashape
python setup.py build_ext --inplace
cd ..

cd nelmin
python setup.py build_ext --inplace
cd ../..


mkdir results
cd results
tar zxf ~/test_data.tar.gz
cd ..


mkdir testDir
cd testDir
cp -R ../media .
cp ../*.cfg .
cd


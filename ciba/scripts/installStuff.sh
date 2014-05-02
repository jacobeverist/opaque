#!/bin/bash

rm /etc/udev/rules.d/70-persistent-net.rules
touch /etc/udev/rules.d/70-persistent-net.rules

rm /lib/udev/rules.d/75-persistent-net-generator.rules
touch /lib/udev/rules.d/75-persistent-net-generator.rules

sed -i".bak" '/HOSTNAME/d' /etc/sysconfig/network
rm /etc/sysconfig/network.bak

sed -i".bak" '/HWADDR/d' /etc/sysconfig/network-scripts/ifcfg-eth0
rm /etc/sysconfig/network-scripts/ifcfg-eth0.bak

sed -i".bak" 's/NM_CONTROLLED=\"yes\"/NM_CONTROLLED=\"no\"/' /etc/sysconfig/network-scripts/ifcfg-eth0
rm /etc/sysconfig/network-scripts/ifcfg-eth0.bak

adduser ctuser
echo "Password2014" | passwd --stdin ctuser 

# add ctuser to sudoers file
echo "ctuser   ALL=(ALL)    ALL" >> /etc/sudoers

# add opencv libraries to library search path
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib" >> /etc/bashrc

echo "AllowTcpForwarding yes" >> /etc/ssh/sshd_config
echo "X11DisplayOffset 10" >> /etc/ssh/sshd_config
echo "X11UseLocalhost yes" >> /etc/ssh/sshd_config

export http_proxy=http://26711:Thu@www-proxy-west.aero.org:8080
export https_proxy=https://26711:Thu@www-proxy-west.aero.org:8080

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






exit

tar zxf pystache-0.5.3.tar.gz
cd pystache-0.5.3
python setup.py build
python setup.py install
cd
rm -rf pystache-0.5.3

cd /opt
tar jxf ~/openmpi-1.4.5.tar.bz2
cd openmpi-1.4.5
./configure --prefix=/opt/openmpi-1.4.5
make all install
cd

cd /opt
tar zxf ~/cmake-2.8.12.2.tar.gz
mv cmake-2.8.12.2 cmake-2.8.12
cd cmake-2.8.12
./configure
make
cd


mv cloud_barometer_2.0.0.tar.bz2 ~ctuser/
chown ctuser ~ctuser/cloud_barometer_2.0.0.tar.bz2
chgrp ctuser ~ctuser/cloud_barometer_2.0.0.tar.bz2

cd ~ctuser


sudo -u ctuser tar jxvf cloud_barometer_2.0.0.tar.bz2
cd cloud_barometer_2.0.0/CTBS
su ctuser -c "source build_env.sh ; cmake . ; make; exit"

cp doc/cb_banner.txt /etc/issue
cp doc/cb_banner.txt /etc/motd
echo "Banner /etc/issue" >> /etc/sshd/sshd_config

echo "export DYNINSTAPI_RT_LIB=/usr/lib64/dyninst/libdyninstAPI_RT.so" >> ~ctuser/.bashrc
chown ctuser ~ctuser/.bashrc
chgrp ctuser ~ctuser/.bashrc

cd ~ctuser/cloud_barometer_2.0.0/traceapp/rpms/dependencies
rpm -ivh libdwarf-20130729-1.el6.x86_64.rpm
rpm -ivh dyninst-8.1.2-4.el6.x86_64.rpm
rpm -ivh dyninst-devel-8.1.2-4.el6.x86_64.rpm
rpm -ivh papi-4.1.3-3.el6.x86_64.rpm
rpm -ivh papi-devel-4.1.3-3.el6.x86_64.rpm
cd ..
rpm -ivh TraceApp-2.0.0-1.x86_64.rpm





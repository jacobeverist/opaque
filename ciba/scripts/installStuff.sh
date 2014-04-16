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


echo "AllowTcpForwarding yes" >> /etc/ssh/sshd_config
echo "X11DisplayOffset 10" >> /etc/ssh/sshd_config
echo "X11UseLocalhost yes" >> /etc/ssh/sshd_config


tar zxvf pystache-0.5.3.tar.gz
cd pystache-0.5.3
python setup.py build
python setup.py install
cd
rm -rf pystache-0.5.3

cd /opt
tar jxvf ~/openmpi-1.4.5.tar.bz2
cd openmpi-1.4.5
./configure --prefix=/opt/openmpi-1.4.5
make all install
cd

cd /opt
tar zxvf ~/cmake-2.8.12.2.tar.gz
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





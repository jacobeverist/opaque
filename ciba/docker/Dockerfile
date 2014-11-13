# DOCKER-VERSION 0.11.1
FROM	centos


MAINTAINER Jacob Everist

# Bundle the dependencies
ADD . /src

#RUN cd /src ; chmod +x ./setupRepo.sh ; ./setupRepo.sh

RUN yum -y update

RUN yum -y install cmake tar make scipy Cython python-matplotlib git wget openssh-server openssh-clients gcc gcc-c++ gcc-gfortran autoconf automake libtool sysstat numpy pygtk2 python-imaging python-simplejson evince boost-devel python-devel lapack-devel freetype-devel libpng-devel unzip gmp-d-evel mpfr-devel vim-X11 xorg-x11-server-Xvfb freeglut-devel bzip2 libXaw-devel libXrandr-devel zziplib-devel zlib-devel

RUN yum -y --skip-broken groupinstall "X Window System"

RUN echo "Xvfb -shmem -screen 0 1280x1024x24 &> /var/log/xvfb.log &" >> /etc/bashrc
RUN echo "export LD_LIBRARY_PATH=/usr/lib:/usr/local/lib" >> /etc/bashrc
RUN echo "export DISPLAY=:0" >> /etc/bashrc

RUN cd /src; chmod +x ./installStuff.sh ; ./installStuff.sh

#CMD while true; echo "foo"; sleep 1; done

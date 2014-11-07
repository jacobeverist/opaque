# DOCKER-VERSION 0.11.1
FROM	centos


MAINTAINER Jacob Everist

# Bundle the dependencies
ADD . /src

RUN cd /src ; chmod +x ./setupRepo.sh ; ./setupRepo.sh

RUN yum update

#RUN yum install -y xvfb git wget openssh-server openssh-clients gcc gcc-c++ gcc-gfortran autoconf automake sysstat numpy pygtk2 python-imaging python-simplejson evince boost-devel python-devel lapack-devel freetype-devel libpng-devel unzip gmp-d-evel mpfr-devel vim-X11

#RUN yum groupinstall -y "X Window System"


#RUN cd /src; chmod +x ./installStuff.sh ; ./installStuff.sh

#CMD while true; echo "foo"; sleep 1; done

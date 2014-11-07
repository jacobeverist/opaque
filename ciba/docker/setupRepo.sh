#!/bin/bash

sed -i".bak" 's/mirrorlist/#mirrorlist/' /etc/yum.repos.d/CentOS-Base.repo
rm /etc/yum.repos.d/CentOS-Base.repo.bak

sed -i".bak" 's/#baseurl/baseurl/' /etc/yum.repos.d/CentOS-Base.repo
rm /etc/yum.repos.d/CentOS-Base.repo.bak

sed -i".bak" 's/mirror.centos.org/www.chapterhouse.aero.org\/mirrors\/sites/' /etc/yum.repos.d/CentOS-Base.repo
rm /etc/yum.repos.d/CentOS-Base.repo.bak

ls -l /etc/yum.repos.d/
cat /etc/yum.repos.d/CentOS-Base.repo

echo "proxy=http://www-proxy-west.aero.org:8080" >> /etc/yum.conf
echo "proxy_username=26711" >> /etc/yum.conf
echo "proxy_password=Thu" >> /etc/yum.conf


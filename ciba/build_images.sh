#!/bin/bash

WORKING_COPY_PATH=/home/everist/segd-cloud

fileName=cloud_barometer_2.0.0

rm cloud_barometer_2.0.0.tar.bz2
rm -rf cloud_barometer_2.0.0

svn export $WORKING_COPY_PATH/dist cloud_barometer_2.0.0
tar cvf cloud_barometer_2.0.0.tar cloud_barometer_2.0.0
bzip2 cloud_barometer_2.0.0.tar

rm -rf output_centos_cloud_barometer

packer/packer build scripts/kvm.json

qemu-img convert -O vmdk -o subformat=monolithicSparse output_centos_cloud_barometer/$fileName.raw output_centos_cloud_barometer/$fileName.vmdk.tmp

vboxmanage clonehd --format=VMDK --variant=Stream output_centos_cloud_barometer/$fileName.vmdk.tmp output_centos_cloud_barometer/$fileName.vmdk
 
rm output_centos_cloud_barometer/$fileName.vmdk.tmp

# BYTE_SIZE = bytes of vmdk image
# VMDK_NAME = name of vmdk image: cloud_barometer_rc1.vmdk
# VIRT_NAME = name of virtual machine: "ctbs-base-ovf"
# VIRT_ID = ID of virtual machine: "ctbs-base-ovf"


rm -rf output_centos_cloud_barometer/$fileName\_ovf
mkdir output_centos_cloud_barometer/$fileName\_ovf

mv output_centos_cloud_barometer/$fileName.vmdk output_centos_cloud_barometer/$fileName\_ovf/$fileName.vmdk
#cp output_centos_cloud_barometer/$fileName.vmdk output_centos_cloud_barometer/$fileName\_ovf/$fileName.vmdk

#capSize="10485760000"
capSize=`ls -l output_centos_cloud_barometer/$fileName.raw | awk '{print $5}'`
fileSize=`ls -l output_centos_cloud_barometer/$fileName\_ovf/$fileName.vmdk | awk '{print $5}'`

sed s/\$BYTE_SIZE/$fileSize/ scripts/ctbs-base-ovf.ovf > output_centos_cloud_barometer/$fileName\_ovf/$fileName.ovf
sed -i".bak" s/\$VMDK_NAME/$fileName.vmdk/ output_centos_cloud_barometer/$fileName\_ovf/$fileName.ovf
sed -i".bak" s/\$VIRT_NAME/$fileName/ output_centos_cloud_barometer/$fileName\_ovf/$fileName.ovf
sed -i".bak" s/\$VIRT_ID/$fileName/ output_centos_cloud_barometer/$fileName\_ovf/$fileName.ovf
sed -i".bak" s/\$CAP_SIZE/$capSize/ output_centos_cloud_barometer/$fileName\_ovf/$fileName.ovf
rm output_centos_cloud_barometer/$fileName\_ovf/$fileName.ovf.bak


# create manifest file
vmdkName="SHA1($fileName.vmdk) ="
vmdkSum=`sha1sum output_centos_cloud_barometer/$fileName\_ovf/$fileName.vmdk | awk '{print $1}'`
vmdkCheck="$vmdkName $vmdkSum"
echo $vmdkCheck > output_centos_cloud_barometer/$fileName\_ovf/$fileName.mf

ovfName="SHA1($fileName.ovf) ="
ovfSum=`sha1sum output_centos_cloud_barometer/$fileName\_ovf/$fileName.ovf | awk '{print $1}'`
ovfCheck="$ovfName $ovfSum"
echo $ovfCheck >> output_centos_cloud_barometer/$fileName\_ovf/$fileName.mf

#exit

rm -rf release-v2.0.0
svn export $WORKING_COPY_PATH/release/v2.0.0 release-v2.0.0
rm -rf release-v2.0.0/labels

cp cloud_barometer_2.0.0.tar.bz2 release-v2.0.0/sw-disk1
md5sum cloud_barometer_2.0.0.tar.bz2 > release-v2.0.0/sw-disk1/cloud_src_checksum

cp output_centos_cloud_barometer/$fileName.raw release-v2.0.0/kvm-disk2/
cd release-v2.0.0/kvm-disk2
#tar cvf kvm.tar $fileName.raw
bzip2 $fileName.raw
md5sum $fileName.raw.bz2 > cloud_raw_checksum
#mv kvm.tar.bz2 release-v2.0.0/kvm-disk2

cd ../..

cp -R output_centos_cloud_barometer/$fileName\_ovf release-v2.0.0/vmware-disk3/
cd release-v2.0.0/vmware-disk3
tar cvf $fileName\_ovf.tar $fileName\_ovf
#bzip2 $fileName\_ovf.tar
#md5sum $fileName\_ovf.tar.bz2 > cloud_ovf_checksum
md5sum $fileName\_ovf.tar > cloud_ovf_checksum
#mv vmware.tar.bz2 release-v2.0.0/vmware-disk3
rm -rf $fileName\_ovf

cd ../..

rm -rf isos
mkdir isos
mkisofs -JR -o isos/barometer_src_2.0.0.iso release-v2.0.0/sw-disk1/
mkisofs -JR -o isos/barometer_kvm_2.0.0.iso release-v2.0.0/kvm-disk2/
mkisofs -JR -o isos/barometer_ovf_2.0.0.iso release-v2.0.0/vmware-disk3/
mkisofs -JR -o isos/barometer_all_2.0.0.iso release-v2.0.0/



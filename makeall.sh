#! /bin/sh

for folder in `ls -d v.*`
do
	echo $folder
	cd $folder
	make MODULE_TOPDIR=/usr/local/grass-7.3.svn
	cd ..
done

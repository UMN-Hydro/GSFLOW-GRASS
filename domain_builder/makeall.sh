#! /bin/sh

for folder in `ls -d */`
do
	echo 
	echo $folder
	echo --------------------------------------
	cd $folder
	make MODULE_TOPDIR=/usr/local/grass-7.3.svn
	cd ..
done


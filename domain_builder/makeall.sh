#! /bin/sh

for folder in `ls -d */`
do
	echo 
	echo $folder
	echo --------------------------------------
	cd $folder
	make MODULE_TOPDIR=/usr/local/grass-7.4.svn
	cd ..
done


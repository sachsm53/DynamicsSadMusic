#!/bin/sh

#The source folder has already done nonlinear registration from highres to standard. 
#The destination folder has already done registration from example_func to highres, but has not gone any further
#The purpose of this script is to avoid repeating registration from highres to standard for each functional run for a given subject

if (( $# < 1 ))
then
	echo "usage: copyreg.sh <source feat folder> <destination feat folder>"
	exit
fi


sourcefolder=$1/reg
destfolder=$2/reg

echo Copying from $sourcefolder to $destfolder

cp $sourcefolder/standard.nii* $destfolder/
cp $sourcefolder/standard_head* $destfolder/
cp $sourcefolder/highres_head* $destfolder/
cp $sourcefolder/standard_mask* $destfolder/
cp $sourcefolder/highres2standard.mat $destfolder/
cp $sourcefolder/standard2highres.mat $destfolder/
cp $sourcefolder/highres2highres_jac* $destfolder/
cp $sourcefolder/highres2standard_warp* $destfolder/
cp $sourcefolder/highres2standard_linear* $destfolder/


convert_xfm -omat $destfolder/example_func2standard.mat -concat $destfolder/highres2standard.mat $destfolder/example_func2highres.mat
convertwarp --ref=$destfolder/standard --premat=$destfolder/example_func2highres.mat --warp1=$destfolder/highres2standard_warp --out=$destfolder/example_func2standard_warp
applywarp --ref=$destfolder/standard --in=$destfolder/example_func --out=$destfolder/example_func2standard --warp=$destfolder/example_func2standard_warp
convert_xfm -inverse -omat $destfolder/standard2example_func.mat $destfolder/example_func2standard.mat
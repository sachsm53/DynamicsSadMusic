#!/usr/bin/python

# Slice timing correction files
# JUNE 2016

import sys,os
from subprocess import call, check_output
import argparse
from datetime import datetime
import re
import shutil
import numpy as np 
import glob

def listdir_nohidden(path):
    return glob.glob(os.path.join(path, '*'))

#command line options
if (len(sys.argv) < 2):
	print "\n\tusage: %s <input subjectfolder>\n" % sys.argv[0]
	sys.exit()
else:
	subjectfolder = sys.argv[1]	

#logging colors
sectionColor = "\033[94m"
sectionColor2 = "\033[96m"
groupColor = "\033[90m"
mainColor = "\033[92m"

pink = '\033[95m'
yellow = '\033[93m'
red = '\033[91m'

ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'

runs = ["rest","happy","sadln","sadsh"]
runs_tr =[360,183,530,271]

#set paths
pathbase = "/Volumes/MusicProject/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-1/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-2/NaturalisticNetwork"

dicompath = pathbase + "/fmri_data"
configfile = dicompath + "/dicomdir/config.json"
analysispath = pathbase + "/fmri_analysis"
designpath = analysispath + "/scripts/designs" 
scriptspath = analysispath + "/scripts" 	


#rename export folder to DICOMr
subject = subjectfolder[-7:-1]
subjectnifti = "%s/%s/"  %(dicompath,subject)
subjectdicom = dicompath + '/dicomdir/' + subject
restfolder = subjectfolder + 'rest/'
funcfolder = subjectfolder + 'music/'


for c,run in enumerate(runs): 
	run_tr = runs_tr[c]
	if c == 0: 
		slicetimefile = restfolder + "%s_slicetime.txt" %(run)
	else: 
		slicetimefile = funcfolder + "%s_slicetime.txt" %(run)

	if os.path.exists(slicetimefile):
		if os.stat(slicetimefile).st_size == 0: 
			red + 'Slice timing for %s for %s is empty. Deleting' %(subject, run)
			os.remove(slicetimefile)
		else: 
			print sectionColor2 + "Slice timing file for %s for %s already exists, moving on %s"  %(subject, run, mainColor)
			test = np.loadtxt(slicetimefile)
			if len(test) != 48:
				print red + 'Slice timing for %s for %s is not the right length' %(subject, run)
				sys.exit()

	if not os.path.exists(slicetimefile):
		txtfile = open(slicetimefile, 'w')
		folder = os.path.join(os.path.abspath(subjectdicom),listdir_nohidden(subjectdicom)[0])
		dicoms = listdir_nohidden(folder)
		for d in dicoms: 
			if len(listdir_nohidden(os.path.join(folder,d))) == run_tr: 
				print yellow + "Creating slicetime file for %s...%s" %(run,mainColor) 
				print yellow + run,run_tr,len(listdir_nohidden(os.path.join(folder,d)))
				dicomone = listdir_nohidden(os.path.join(folder,d))[0]
				dicomonepath = os.path.join(folder,d,dicomone)

				#get slice timing files from afni: 
				#command = "dicom_hdr -slice_times %s > %s" %(dicomonepath,slicetimefile)
				command = "dicom_hdr -slice_times %s" %(dicomonepath)
				out = check_output(command, shell = True).split()
				#print out
				for i,item in enumerate(out):
					if i > 4:
						num = float(item)/(1000)
						txtfile.write("%f\n" %num)







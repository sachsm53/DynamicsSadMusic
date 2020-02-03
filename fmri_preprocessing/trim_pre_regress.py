#!/usr/bin/env python

#Preprocessing script for music naturalistic network study (June 2018)

import sys,os
from subprocess import call, check_output
import argparse
from datetime import datetime
import re
import shutil
import numpy as np 

starttime = datetime.now()


##CHANGE THIS HERE
songs = ["sadln"]


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


#Handle command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--subjects",help="process listed subjects",nargs='+',action="store")
parser.add_argument("--all",help="process all subjects", action="store_true")
parser.add_argument("--half",help="process half subjects", action="store_true")
parser.add_argument("--unzip",help="unzip file and move to ISC", action="store_true")
parser.add_argument("--aroma",help="Do same for aroma file", action="store_true")
parser.add_argument("--cut20",help="perform the same on the cut20 files in ", action="store_true")
args = parser.parse_args()

subjects = args.subjects

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
logfilename = analysispath + "/logs/analysis_log.txt"

#os.environ["FSLOUTPUTTYPE"] = "NIFTI"


#develop list of subjects
subjects = args.subjects
removelist = ['sub-07','sub-28','sub-30']

if args.all:
	#Get list of subjects
	subjects = os.listdir(dicompath + '/dicomdir/')
	subjects = [elem for elem in subjects if "sub" in elem]
	subjects = [elem for elem in subjects if "18" not in elem]
	subjects = [elem for elem in subjects if "07" not in elem]
	subjects = [elem for elem in subjects if "28" not in elem]
	subjects = [elem for elem in subjects if "30" not in elem]
	subjects = [elem for elem in subjects if "pil" not in elem]
	subjects.sort()

if args.half:
	#Get list of subjects
	subjects = os.listdir(dicompath + '/dicomdir/')
	subjects = [elem for elem in subjects if "sub" in elem]
	subjects = [elem for elem in subjects if "07" not in elem]
	subjects = [elem for elem in subjects if "pil" not in elem]
	subjects.sort()
	subjects = subjects[0:19]


if subjects:
	print subjects
else:
	print "Subjects must be specified. Use --all for all subjects or --subjects to list specific subjects."
	sys.exit()

##DEFINE THE MODEL HERE:
#model = 'hpf128'
model = '200hpf'


for subject in subjects:

	#print subject
	#rename export folder to DICOM
	subjectfolder = "%s/%s/"  %(analysispath,subject)
	mpragebrain = subjectfolder + "mprage_brain.nii.gz"
	mprage = subjectfolder + "/mprage.nii.gz"
	funcfolder = subjectfolder + 'music/'
	restfolder = subjectfolder + 'rest/'


	for song in songs:
		featfile = "%s%s_model1_pre_%s.feat/" % (funcfolder,song,model)
		restfile = "%sresting_pre.feat/" %(restfolder)
		standard = "%sreg/standard.nii.gz" % restfile
		premat = "%sreg/example_func2highres.mat" % (featfile)
		warpfile = "%sreg/highres2standard_warp.nii.gz" % (restfile)
		if not os.path.exists(warpfile):
			print red + 'Could not find %s. Moving on...%s' %(warpfile, mainColor)
			continue
		#Make 200sec temporal filtering folder within feat_model1
		#newfolder = "%slpf_200/" % (featfile)
		
		if not args.cut20:
			#Full trim file (not cutting the first 20)
			if not args.aroma:

				tffile = "%sfiltered_func_data_trim.nii.gz" % (featfile)
				inputdata = "%sfiltered_func_data_%s_residuals.nii.gz" %(featfile,model)
				outputdata = "%sfiltered_func_data_%s_residuals_standard.nii.gz" % (featfile,model)

				if not os.path.exists(tffile):
					print red + 'Could not find %s. Moving on...%s' %(tffile, mainColor)
					continue

				dst = '%s/ISC/data/%s_%s_filtered_func_%s_standard.nii.gz' %(pathbase,subject,song,model)
				dstunzip = '%s/ISC/data/%s_%s_filtered_func_%s_standard.nii' %(pathbase,subject,song,model)
				if os.path.exists(dstunzip):
					print yellow + 'Warping and unzipping already completed for %s. Moving on%s' %(dstunzip, mainColor)
					continue


				designfile = "%sresid_design_trim_%s.txt" % (featfile,model)
				if not os.path.exists(designfile):
					print red + 'Could not find %s. Moving on...%s' %(designfile, mainColor)
					continue

			#ADD AROMA HERE
			elif args.aroma:
				#Full trim file (not cutting the first 20)
				tffile = "%sfiltered_func_data_aroma_trim.nii.gz" % (featfile)
				inputdata = "%sfiltered_func_data_%s_aroma_residuals.nii.gz" %(featfile,model)
				outputdata = "%sfiltered_func_data_%s_aroma_residuals_standard.nii.gz" % (featfile,model)

				if not os.path.exists(tffile):
					print red + 'Could not find %s. Moving on...%s' %(tffile, mainColor)
					continue

				designfile = "%sresid_design_trim_%s_aroma.txt" % (featfile,model)
				if not os.path.exists(designfile):
					print red + 'Could not find %s. Moving on...%s' %(designfile, mainColor)
					continue

				dst = '%s/ISC/data/%s_%s_filtered_func_%s_standard_aroma.nii.gz' %(pathbase,subject,song,model)
				dstunzip = '%s/ISC/data/%s_%s_filtered_func_%s_standard_aroma.nii' %(pathbase,subject,song,model)
				if os.path.exists(dstunzip):
					print yellow + 'Warping and unzipping already completed for %s. Moving on%s' %(dstunzip, mainColor)
					continue

		elif args.cut20:
			#Trim file (not cutting the first 20)
			tffile = "%sfiltered_func_data_cut20.nii.gz" % (featfile)
			inputdata = "%sfiltered_func_data_%s_cut20_residuals.nii.gz" %(featfile,model)
			outputdata = "%sfiltered_func_data_%s_cut20_residuals_standard.nii.gz" % (featfile,model)
			if not os.path.exists(tffile):
				print red + 'Could not find %s. Moving on...%s' %(tffile, mainColor)
				continue

			dst = '%s/ISC/data/%s_%s_filtered_func_%s_cut20_standard.nii.gz' %(pathbase,subject,song,model)
			dstunzip = '%s/ISC/data/%s_%s_filtered_func_%s_cut20_standard.nii' %(pathbase,subject,song,model)
			if os.path.exists(dstunzip):
				print yellow + 'Warping and unzipping already completed for %s. Moving on%s' %(dstunzip, mainColor)
				continue

			designfile = "%sresid_design_%s_cut20.txt" % (featfile,model)
			if not os.path.exists(designfile):
				print red + 'Could not find %s. Moving on...%s' %(designfile, mainColor)
				continue



		#FIrst regress out confounds
		if not os.path.exists(inputdata):
			command = "fsl_glm -i %s --demean -m %smask.nii.gz -d %s --out_res=%s" % (tffile,featfile,designfile,inputdata)
			print sectionColor2 + 'Regressing out confounds for %s%s' %(tffile, mainColor)
			call(command,shell=True)
		else:
			print yellow + 'Regressing out confounds already completed for %s%s' %(subject, mainColor)
		

		#Next convert to standard space
		if not os.path.exists(outputdata):
			command = "applywarp --ref=%s --in=%s --out=%s --warp=%s --premat=%s" % (standard,inputdata,outputdata,warpfile,premat)
			print sectionColor2 + 'Warping to standard space. Output is %s%s' %(outputdata, mainColor)
			call(command,shell=True)
		else:
			print yellow + 'Already converted to standard for %s%s' %(subject, mainColor)
			print outputdata
		
		print 

		#Unzip here and put in one file for ISC
		if args.unzip:
			if os.path.exists(outputdata):
				scr = outputdata
				if not os.path.exists(dst):
					print 'Found file. Copying file from %s to %s' %(scr,dst)
					shutil.copy2(scr,dst)

				if not os.path.exists(dstunzip):
					print 'Unzipping %s' %(dst)
					command = 'gunzip %s' %(dst)
					call(command,shell=True)
					continue
				else:
					print dstunzip
					continue

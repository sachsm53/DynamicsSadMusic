#!/usr/bin/env python


import sys,os
from subprocess import call, check_output
import argparse
from datetime import datetime
import re
import shutil
from commando import commando
import numpy as np 
import pandas as pd

starttime = datetime.now()

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

#EX: ./pairwise_cor_afni.py --aroma --stand --subgroup

#command line optionss
parser = argparse.ArgumentParser()
parser.add_argument("--aroma",help="process aroma filtered data", action="store_true")
parser.add_argument("--stand",help="process standard space images", action="store_true")
parser.add_argument("--song",help="select which song you want to run",nargs='+',action="store")
args = parser.parse_args()

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

#Make network cohesion directory 
outputdir = os.path.join(pathbase,'group_compare/pairwise_mats')


def checkImageLength(imagename):
	command = 'fslinfo %s' % imagename
	results = check_output(command,shell=True)
	TR = results.split()[9]
	return int(TR)


#develop list of subjects
removelist = ['sub-07','sub-18','sub-28','sub-30']

#Get list of subjects
subjects = os.listdir(dicompath + '/dicomdir/')
subjects = [elem for elem in subjects if "sub" in elem]
subjects = [elem for elem in subjects if "pil" not in elem]
subjects = [elem for elem in subjects if elem not in removelist]
subjects.sort()

if args.stand: 
	textfile = '%s/group_compare/pairwise_data_%dss_aroma_stand.txt' %(pathbase,len(subjects))
else:
	textfile = '%s/group_compare/pairwise_data_%dss_aroma.txt' %(pathbase,len(subjects))

#songs = ["sadln","happy"]
song = "sadln"
song_timepoints =[183,530,271]

count = 0
textout = []
#Loop through subjects
for subject in subjects:

	subjectfolder = "%s/%s/"  %(analysispath,subject)
	funcfolder = subjectfolder + "music/%s_model1_pre_200hpf.feat" %(song)
	if args.aroma: 
		funcimage=  '%s/filtered_func_data_200hpf_aroma_residuals.nii.gz' %funcfolder
		funcimagestand = '%s/filtered_func_data_200hpf_aroma_residuals_standard.nii.gz' %funcfolder

	elif not args.aroma:
		funcimage=  '%s/filtered_func_data_200hpf_cut20_residuals.nii.gz' %funcfolder
		funcimagestand = '%s/filtered_func_data_200hpf_cut20_residuals_standard.nii.gz' %funcfolder
		textfile = '%s/group_compare/pairwise_data_40ss_noaroma.txt' %(pathbase)

	for sub2 in subjects: 
		if subject != sub2: 
			count = count + 1
			subjectfolder2 = "%s/%s/"  %(analysispath,sub2)
			funcfolder2 = subjectfolder2 + "/music/%s_model1_pre_200hpf.feat" %(song)
			if args.aroma: 
				funcimage2 = '%s/filtered_func_data_200hpf_aroma_residuals.nii.gz' %funcfolder2
				funcimagestand2 = '%s/filtered_func_data_200hpf_aroma_residuals_standard.nii.gz' %funcfolder2
				outpre = '%s/%s_%s_aroma' %(outputdir,subject, sub2)
			elif not args.aroma:
				funcimage2 =  '%s/filtered_func_data_200hpf_cut20_residuals.nii.gz' %funcfolder2
				funcimagestand2 = '%s/filtered_func_data_200hpf_cut20_residuals_standard.nii.gz' %funcfolder2
				outpre = '%s/%s_%s' %(outputdir,subject, sub2)
			
			#Run pairwise correlation
			if not args.stand:
				checkfile = '%s+orig.HEAD' %outpre
				if not os.path.exists(checkfile):
					print sectionColor2 + '%d: Correlation %s with %s and saving as %s%s' %(count,subject, sub2, outpre,mainColor)
					command = '3dTcorrelate -prefix %s %s %s' %(outpre,funcimage,funcimage2)
					print command
					call(command, shell = True)
				else:
					print yellow + '%d: Correlation completed for %s with %s. Moving on...%s' %(count,subject, sub2,mainColor)
			elif args.stand:
				checkfile = '%s+tlrc.HEAD' %outpre
				if not os.path.exists(checkfile):
					print sectionColor2 + '%d: Correlation %s with %s and saving as %s%s' %(count,subject, sub2, outpre,mainColor)
					command = '3dTcorrelate -prefix %s %s %s' %(outpre,funcimagestand,funcimagestand2)
					print command
					call(command, shell = True)
				else:
					print yellow + '%d: Correlation completed for %s with %s. Moving on...%s' %(count,subject, sub2,mainColor)

			#Create row
			textrow = [subject,sub2,checkfile]
			if count == 1:
				textout = textrow
			else:
				textout = np.vstack((textout,textrow))
			print textrow

#Write to text file
print textfile,textout 

if not os.path.exists(textfile):
	df = pd.DataFrame(textout)
	df.to_csv(textfile,sep='\t', header=['Subj', 'Subj2', 'InputFile'],index=False)
	print sectionColor2 + 'Saving as %s' %(textfile)












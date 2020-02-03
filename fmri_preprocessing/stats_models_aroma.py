#!/usr/bin/python

#Preprocessing script for music naturalistic network study (June 2018)

import sys,os
from subprocess import call, check_output
import argparse
from datetime import datetime
import re
import shutil
import numpy as np 
import distutils 


#logging colors
sectionColor = "\033[94m"
sectionColor2 = "\033[96m"
groupColor = "\033[90m"
mainColor = "\033[92m"

pink = '\033[95m'
yellow = '\033[93m'
red = '\033[91m'


#Handle command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--subjects",help="process listed subjects",nargs='+',action="store")
parser.add_argument("--all",help="process all subjects", action="store_true")
parser.add_argument("--half",help="process half subjects", action="store_true")
parser.add_argument("--model",help="process listed subjects",nargs='+',action="store") #added a model argument on 10/19
args = parser.parse_args()

subjects = args.subjects

#analysislist = ["generic_resting_pre","generic_pre_music"]
runs = [2]
songs = ["happy","sadln","sadsh"]
#songs = ["SL"]
song_timepoints =[183,495,271] #changed to reflect 20 TRs cut from beginning

#set paths
pathbase = "/Volumes/MusicProject/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-1/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-2/NaturalisticNetwork"

analysispath = pathbase + "/fmri_analysis"
designpath = analysispath + "/scripts/designs" 
scriptspath = analysispath + "/scripts"

#develop list of subjects
subjects = args.subjects
if args.all:
	#Get list of subjects
	subjects = os.listdir(analysispath)
	subjects = [elem for elem in subjects if "sub" in elem]
	subjects = [elem for elem in subjects if "07" not in elem]
	subjects = [elem for elem in subjects if "28" not in elem]
	subjects = [elem for elem in subjects if "pil" not in elem]
	subjects.sort()

if args.half:
	#Get list of subjects
	subjects = os.listdir(analysispath)
	subjects = [elem for elem in subjects if "sub" in elem]
	subjects = [elem for elem in subjects if "07" not in elem]
	subjects = [elem for elem in subjects if "pil" not in elem]
	subjects.sort()
	subjects = subjects[1:19]


if subjects:
	print subjects
else:
	print "Subjects must be specified. Use --all for all subjects or --subjects to list specific subjects."
	sys.exit()

model = args.model
if model:
	print 'Model specified: %s' %(model[0])
	model = int((model[0]))
else:
	print "No Model specified. Running both '1' or '2'."
	model = 3 #DO BOTH
	#sys.exit()


def checkImageLength(imagename):
	command = 'fslinfo %s' % imagename
	results = check_output(command,shell=True)
	TR = results.split()[9]
	return int(TR)

def yes_or_no(question):
    reply = str(raw_input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        return True
    if reply[0] == 'n':
        return False
    else:
        return yes_or_no("Uhhhh... please enter ")

for subject in subjects:

	#rename export folder to DICOM
	subjectfolder = "%s/%s/"  %(analysispath,subject)
	restfolder = subjectfolder + 'rest/'
	restingreg = '%sresting_pre.feat/reg' %(restfolder)
	funcfolder = subjectfolder + 'music/'

	for run in runs:
		song = songs[run-1]
		#filename = "%s%s_model1_pre.feat/filtered_func_data_.nii.gz" % (funcfolder,song)
		prefilename = "%s%s_model1_pre_200hpf.feat/filtered_func_data_200hpf_aroma_residuals_cut20_add.nii.gz" % (funcfolder,song)
		if not os.path.exists(prefilename):
			print red + '%s does not exist. Exiting script' %(prefilename)
			continue

		timepoints = song_timepoints[run-1]
		run = str(run)
		timepoints = str(timepoints)
		print sectionColor + 'Working on %s %s %s %s' %(subject,run,song,timepoints) 


		length = checkImageLength(prefilename)
		#print yellow + "%s length: %d%s" % (song,length,mainColor)
		if str(length) != timepoints:
			print red + "Song %s has incorrect file length. Check files.%s" %(song,mainColor)
			continue

		if model == 1 or model == 3:

			checkfolder1 = funcfolder + '%s_stats_aroma_cut20_allratings_200hpf_emo_wrt_enj.feat' %(song)
			checkfolder2 = funcfolder + '%s_stats_aroma_cut20_allratings_200hpf_enj_wrt_emo.feat' %(song)
			#checkfile = funcfolder + '%s_model1_stats_nopnm_cut20_200hpf_aroma.feat/stats/cope1.nii.gz' %(song)
			genericfile1 = designpath + "/" + 'generic_stats_design_aroma_enj_wrt_emo.fsf'
			genericfile2 =designpath + "/"  + 'generic_stats_design_aroma_emo_wrt_enj.fsf'
			outputfile1 = funcfolder + 'design_music_%s_model1_stats_design_aroma_emo_wrt_enj.fsf' %(song)
			outputfile2 = funcfolder + 'design_music_%s_model2_stats_design_aroma_enj_wrt_emo.fsf' % (song)
			checkregfolder1 = '%s/reg' %(checkfolder1)
			checkregfolder2 ='%s/reg' %(checkfolder2)
			#checkregfile = '%s/reg/standard.nii.gz' %(checkfolder)


			# if os.path.exists(checkfolder):
			# 	print red + 'Deleting old stats level feat folder %s because of temporal filtering errors...%s' %(checkfolder,mainColor)
			# 	shutil.rmtree(checkfolder)


			print sectionColor2 + 'Starting Model 1 stats feat with temporal derivates for %s...%s' %(song,mainColor)
			command = "sed -e 's/DEFINEBASE/%s/g' -e 's/DEFINESONG/%s/g' -e 's/DEFINESUBJECT/%s/g' -e 's/DEFINETIMEPOINTS/%s/g' %s > %s" % (re.escape(analysispath),song,subject,timepoints,genericfile1,outputfile1)
			call(command,shell=True)

			command = "feat " + outputfile1
			print command
			call(command,shell=True)
			prereg = "%s%s_model1_pre_200hpf.feat/reg" % (funcfolder,song)
			print sectionColor2 + 'Copying preprocessed feat reg folder from %s%s' %(prereg,mainColor)
			shutil.copytree(prereg,checkregfolder1)

			print sectionColor2 + 'Starting Model 2 stats feat with temporal derivates for %s...%s' %(song,mainColor)
			command = "sed -e 's/DEFINEBASE/%s/g' -e 's/DEFINESONG/%s/g' -e 's/DEFINESUBJECT/%s/g' -e 's/DEFINETIMEPOINTS/%s/g' %s > %s" % (re.escape(analysispath),song,subject,timepoints,genericfile2,outputfile2)
			call(command,shell=True)

			command = "feat " + outputfile2
			print command
			call(command,shell=True)
			print sectionColor2 + 'Copying preprocessed feat reg folder from %s%s' %(prereg,mainColor)
			shutil.copytree(prereg,checkregfolder2)



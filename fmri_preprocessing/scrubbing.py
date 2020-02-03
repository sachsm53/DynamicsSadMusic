		
#!/usr/bin/python

#Preprocessing script for music naturalistic network study (June 2018)

import sys,os
from subprocess import call, check_output
import argparse
from datetime import datetime
import re
import shutil
from commando import commando
import numpy as np 

starttime = datetime.now()

#analysislist = ["generic_resting_pre","generic_pre_music"]
runs = [1]
songs = ["happy","sadln","sadsh"]
#songs = ["SL"]
song_timepoints =[183,530,271]

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

#command line options
parser = argparse.ArgumentParser()

parser.add_argument("--subjects",help="process listed subjects",nargs='+',action="store")
parser.add_argument("--all",help="process all subjects", action="store_true")
parser.add_argument("--norest",help="skip resting analysis", action="store_true")
parser.add_argument("--norestreg",help="dont do registration to standard space for resting analysis", action="store_true")
parser.add_argument("--nopre",help="skip all preprocessing steps", action="store_true")
parser.add_argument("--nomusic",help="skip all analysis steps", action="store_true")
parser.add_argument("--nodcm",help="skip dicom conversion", action="store_true")
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
skullstripscript = analysispath + "/scripts/skullstrip.py"
slicetimescript = analysispath + "/scripts/slicetimefiles_snl.py"


subjects = os.listdir(dicompath + '/dicomdir/')
subjects = [elem for elem in subjects if ".DS" not in elem]
subjects = [elem for elem in subjects if "pil" not in elem]
subjects = [elem for elem in subjects if "07" not in elem]
subjects.sort()


#Preprocessing steps
for subject in subjects:

	#rename export folder to DICOM
	subjectfolder = "%s/%s/"  %(analysispath,subject)
	mpragebrain = subjectfolder + "mprage_brain.nii.gz"
	mprage = subjectfolder + "/mprage.nii.gz"
	funcfolder = subjectfolder + 'music/'


	for run in runs: 
		song = songs[run]
		origdatafile = funcfolder + '%s.nii.gz' %(song) 
		scrubout = funcfolder + 'scrub_fd_%s' %(song)
		metric_values_text = funcfolder + "metric_values_fd_%s" %(song)
		metric_values_plot = funcfolder + "metric_plot_fd_%s" %(song)

		if os.path.exists(scrubout):
			print yellow + "FSL Motion Outliers already completed for %s. Moving on\n%s"  % (subject,mainColor)

		if not os.path.exists(metric_values_text):
			print sectionColor2 + " Scrubbing for %s to determine mFD. Output is %s\n%s"  % (subject,scrubout,mainColor)
			command = "fsl_motion_outliers -i %s -o %s --fd --thresh=%s -s %s -p %s -v" % (origdatafile, scrubout, "0.5", metric_values_text, metric_values_plot)
			call(command, shell = True)

		else: 
			print yellow + "FSL Motion Outliers already completed for %s. Moving on\n%s"  % (subject,mainColor)





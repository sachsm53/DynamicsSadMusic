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
parser.add_argument("--half",help="process half subjects", action="store_true")
parser.add_argument("--norest",help="skip resting analysis", action="store_true")
parser.add_argument("--norestreg",help="dont do registration to standard space for resting analysis", action="store_true")
parser.add_argument("--nomusic",help="skip all analysis steps", action="store_true")
parser.add_argument("--nodcm",help="skip dicom conversion", action="store_true")
parser.add_argument("--nofieldmap",help="skip fieldmap prep", action="store_true")
parser.add_argument("--noskullstrip",help="skip skullstripping", action="store_true")
parser.add_argument("--nofeat",help="skip feat analysis", action="store_true")
parser.add_argument("--noaroma",help="skip ICA-AROMA", action="store_true")
parser.add_argument("--noseg",help="skip segmentation", action="store_true")
parser.add_argument("--noconfound",help="skip regressing out confound", action="store_true")
parser.add_argument("--notrim",help="skip trimming the data", action="store_true")
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
icascript = scriptspath + "/ICA-AROMA/ICA_AROMA.py"

#develop list of subjects
subjects = args.subjects

if args.all:
	#Get list of subjects
	subjects = os.listdir(dicompath + '/dicomdir/')
	subjects = [elem for elem in subjects if "sub" in elem]
	subjects = [elem for elem in subjects if "07" not in elem]
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


runs = [1] #CHANGE THIS IF YOU WANT TO RUN ALL PIECES 
songs = ["happy","sadln","sadsh"]
#songs = ["SL"]
song_timepoints =[183,530,271]
song = songs[runs[0]-1]

if subjects:
	print subjects
else:
	print "Subjects must be specified. Use --all for all subjects or --subjects to list specific subjects."
	sys.exit()

if len(runs) == 1:
	print songs[runs[0] - 1]
else:
	print "No song selected. Looping through all songs"
	#sys.exit()


#######################################################################
## DEFINE THE FUNCTIONS
#######################################################################

def checkfeat(featfolder):
	testfile = featfolder + "/filtered_func_data.nii.gz" 
	if not os.path.exists(testfile):
		print red + "WARNING: ANALYSIS DID NOT COMPLETE FOR %s%s" %(featfolder, mainColor)
		logfile.write("%s: WARNING: ANALYSIS DID NOT COMPLETE FOR %s\n" % (datetime.now().strftime('%I:%M:%S%p'),featfolder))
		sys.exit()


def doresting(designsuffix):
	
	#move run
	restfile = "%sresting.nii.gz" %(restfolder)
	dicomrestfile = "%s/%s/func/%s_task-rest_bold.nii.gz" % (dicompath,subject,subject)

	if not os.path.exists(restfile):
		shutil.copy2(dicomrestfile,restfile)

	timepoints = str(checkImageLength(restfile))

	if (designsuffix):
		genericfile = designpath + "/" + designsuffix + ".fsf"
		outputfile = restfolder + designsuffix + ".fsf"
		testdir = "%sresting_pre.feat" % (restfolder)
		testfile = testdir + "/filtered_func_data.nii.gz"
		if not os.path.exists(testfile):	
			command = "sed -e 's/DEFINESUBJECT/%s/g' -e 's/DEFINETIMEPOINTS/%s/g' %s > %s" % (subject,timepoints,genericfile,outputfile)
			call(command,shell=True)

			command = "feat " + outputfile
			print command
			call(command,shell=True)
		else: 
			print sectionColor2 + 'Resting state preprocessing feat directory already exists for %s. Moving on.%s' %(subject,mainColor)

	checkfeat(testdir)


def dolowerlevels(designsuffix):
	for run in runs:
		song = songs[run-1]
		timepoints = song_timepoints[run-1]
		run = str(run)
		timepoints = str(timepoints)
		#print run,song,timepoints

		#move run
		filename = "%s%s.nii.gz" % (funcfolder,song)
		niftifile = "%s/%s/func/%s_task-%s_bold.nii.gz" % (dicompath,subject,subject,song)

		if not os.path.exists(filename):
			shutil.copy2(niftifile,filename)

		length = checkImageLength(filename)
		print yellow + "%s length: %d%s" % (song,length,mainColor)
		if str(length) != timepoints:
			print red + "Song %s has incorrect file length. Check files.%s" %(song,mainColor)
			continue

		if (designsuffix):
			genericfile = designpath + "/" + designsuffix + ".fsf"
			outputfile = funcfolder + 'design_music_%s_200hpf.fsf' %song
			testdir = "%s%s_model1_pre_200hpf.feat" % (funcfolder,song)
			testfile = testdir + "/filtered_func_data.nii.gz"
			if not os.path.exists(testfile):
				command = "sed -e 's/DEFINEBASE/%s/g' -e 's/DEFINESONG/%s/g' -e 's/DEFINESUBJECT/%s/g' -e 's/DEFINETIMEPOINTS/%s/g' %s > %s" % (re.escape(analysispath),song,subject,timepoints,genericfile,outputfile)
				call(command,shell=True)

				command = "feat " + outputfile
				print command
				call(command,shell=True)

			else: 
				print sectionColor2 + 'Preprocessing feat directory already exists for %s. Moving on.%s' %(song,mainColor)
		
		#Copying registration from resting reg folder
		checkfeat(testdir)
		restingreg = '%sresting_pre.feat/reg' %(restfolder)
		checkregfile = '%s%s_model1_pre_200hpf.feat/reg/standard.nii.gz' %(funcfolder,song)

		if os.path.exists(restingreg):
			if not os.path.exists(checkregfile): 
				print sectionColor + 'Copying reg folder from resting state %s%s' %(subject,mainColor)
				command = "%s/copyreg.sh %sresting_pre.feat %s%s_model1_pre_200hpf.feat" % (scriptspath,restfolder,funcfolder,song)
				call(command,shell=True)

		else:
			print red + 'Cannot copy reg directory because resting state analysis not completed for %s%s' %(subject,mainColor)


def doseg(mpragebrain,subjectfolder):

	print sectionColor + "----- Segmentation for %s%s" %(subject,mainColor)
	segfile = '%smprage_brain_seg_0.nii.gz' %(subjectfolder)


	if not os.path.exists(segfile):
		command = "fast -g -o %s %s" % (mpragebrain,mpragebrain)
		print command
		call(command,shell=True)
		command="slicer %s %smprage_brain_seg_0 -a %sCSF.png" % (mpragebrain,subjectfolder,subjectfolder)
		call(command,shell=True)
		command="slicer %s %smprage_brain_seg_1 -a %sGM.png" % (mpragebrain,subjectfolder,subjectfolder)
		call(command,shell=True)
		command="slicer %s %smprage_brain_seg_2 -a %sWM.png" % (mpragebrain,subjectfolder,subjectfolder)
		call(command,shell=True)
		#writeToLog("<h2>Segmentation</h2><br>CSF:<br><img src=CSF.png><br><br>White matter:<br><img src=WM.png><br><br>Gray matter:<br><img src=GM.png><br><br><hr>",reportfile)
	else: 
		print sectionColor2 + "Segmentation already completed for %s. Moving on\n%s"  % (subject,mainColor) 


def doconfound(task):

	# CSF is mprage_brain_seg_0 is and WM is mprage_brain_seg_2
	print sectionColor + "----- Extracting confound timeseries" + mainColor
	
	for run in runs:
		song = songs[run-1]
		timepoints = song_timepoints[run-1]
		run = str(run)
		timepoints = str(timepoints)

		featfile = "%s%s_model1_pre_200hpf.feat/" % (funcfolder,song)
		
		if not args.noaroma: 
			inputfile = "%sfiltered_func_data_aroma.nii.gz" % (featfile)
			csffile = featfile + "csf_200hpf_aroma.txt"
			wmfile = featfile + "wm_200hpf_aroma.txt"

		if args.noaroma:
			inputfile = "%sfiltered_func_data.nii.gz" % (featfile)
			csffile = featfile + "csf_200hpf.txt"
			wmfile = featfile + "wm_200hpf.txt"

		if not os.path.exists(csffile):
			print yellow + "Exracting time series for confounds for %s%s"  % (song,mainColor)			
			command = "flirt -in %smprage_brain_seg_0 -ref %s/example_func.nii.gz -applyxfm -init %s/reg/highres2example_func.mat -interp nearestneighbour -o %sCSF_200hpf.nii.gz" % (subjectfolder,featfile,featfile,featfile)
			call(command,shell = True)
			command = "flirt -in %smprage_brain_seg_2 -ref %s/example_func.nii.gz -applyxfm -init %s/reg/highres2example_func.mat -interp nearestneighbour -o %sWM_200hpf.nii.gz" % (subjectfolder,featfile,featfile,featfile)
			call(command,shell = True)
			command = "fslmeants -i %s -m %sCSF_200hpf.nii.gz -o %s" % (inputfile,featfile,csffile)
			call(command,shell = True)
			command = "fslmeants -i %s -m %sWM_200hpf.nii.gz -o %s" % (inputfile,featfile,wmfile)
			call(command,shell = True)
		else: 
			print sectionColor2 + "Confound timeseries already created for %s. Moving on\n%s"  % (subject,mainColor)

##This cuts off the last 5 seconds from the file because it is long
def dotrim(task):
	print sectionColor + "----- Trimming data and confounds" + mainColor
	startvol = 0
	for run in runs:
		song = songs[run-1]
		timepoints = song_timepoints[run-1]
		run = str(run)
		timepoints = str(timepoints)

		#1. Trimming data file
		featfile = "%s%s_model1_pre_200hpf.feat/" % (funcfolder,song)

		if not args.noaroma:
			infile = "%sfiltered_func_data_aroma.nii.gz" % (featfile)
			outfile = "%sfiltered_func_data_aroma_trim.nii.gz" % (featfile)
			wmfile = "%swm_200hpf_aroma.txt" %(featfile)
			wmfiletrim = "%swm_200hpf_aroma_trim.txt" %(featfile)
			csffile = "%scsf_200hpf_aroma.txt" %(featfile)
			csffiletrim = "%scsf_200hpf_aroma_trim.txt" %(featfile)
			designfile = "%sresid_design_trim_200hpf_aroma.txt" % (featfile)

		if args.noaroma:
			infile = "%sfiltered_func_data.nii.gz" % (featfile)
			outfile = "%sfiltered_func_data_trim.nii.gz" % (featfile)
			wmfile = "%swm_200hpf.txt" %(featfile)
			wmfiletrim = "%swm_200hpf_trim.txt" %(featfile)
			csffile = "%scsf_200hpf.txt" %(featfile)
			csffiletrim = "%scsf_200hpf_trim.txt" %(featfile)
			designfile = "%sresid_design_trim_200hpf.txt" % (featfile)

		totalvols  = checkImageLength(infile)
		trimvols = totalvols - 5

		
		if not os.path.exists(outfile):
			#10 seconds have already been trimmed in preprocessing.
			#Now we want to remove an additional 5 from the end. 
			print yellow + "Trimming data from %d to %d for %s%s"  % (startvol,trimvols,song,mainColor)
			command = "fslroi %s %s %s %s" % (infile,outfile,startvol,trimvols)
			print command
			call(command,shell=True)
		else:
			print sectionColor2 + '%s already exists. Moving on' %(outfile)

		#2. Trim the movement files
		mcfile = "%smc/prefiltered_func_data_mcf.par" %(featfile)
		mcfiletrim = "%smc/prefiltered_func_data_mcf_trimmed.par" %(featfile)

		if not os.path.exists(mcfiletrim):
			print yellow + "Trimming motion parameters for %s%s"  % (song,mainColor)
			with open(mcfile,'r') as mc:
				lines = mc.readlines()
			trimmed_lines = lines[startvol:startvol+trimvols]
			with open(mcfiletrim,'w') as nf:
				nf.write(''.join(trimmed_lines))
		else: 
			print sectionColor2 + '%s already exists. Moving on' %(mcfiletrim)


		#3. Trim the white matter and CSF files
		
		if not os.path.exists(wmfiletrim):
			print yellow + "Trimming wm for %s%s"  % (song,mainColor)
			with open(wmfile,'r') as wm:
				lines = wm.readlines()
			trimmed_lines = lines[startvol:startvol+trimvols]
			with open(wmfiletrim,'w') as wnf:
				wnf.write(''.join(trimmed_lines))
		else: 
			print sectionColor2 + '%s already exists. Moving on' %(wmfiletrim)

		if not os.path.exists(csffiletrim):
			print yellow + "Trimming csf for %s%s"  % (song,mainColor)
			with open(csffile ,'r') as csf:
				lines = csf.readlines()
			trimmed_lines = lines[startvol:startvol+trimvols]
			with open(csffiletrim,'w') as cnf:
				cnf.write(''.join(trimmed_lines))
		else: 
			print sectionColor2 + '%s already exists. Moving on' %(csffiletrim)

		#4. If exists, trim the scrub file 
		scrubout = funcfolder + 'scrub_fd_%s' %(song)
		scrubouttrim = funcfolder + 'scrub_fd_%s_trim.txt' %(song)
		if os.path.exists(scrubout):
			if not os.path.exists(scrubouttrim): 
				print yellow + "Trimming scrub file for %s%s"  % (song,mainColor)
				with open(scrubout ,'r') as scrub:
					lines = scrub.readlines()
				trimmed_lines = lines[startvol:startvol+trimvols]
				with open(scrubouttrim,'w') as snf:
					snf.write(''.join(trimmed_lines))
			else: 
				print sectionColor2 + '%s already exists. Moving on' %(scrubouttrim)
		else: 
			print yellow + 'No scrubbing need for %s. Moving on' %(subject)


		############ 
		# CREATE DESIGN FILE: Take this out for now because more trimming coming IN LATER STEP 
		#################


		#Construct design.mat
		if not os.path.exists(designfile):
			mcparams = np.loadtxt(mcfiletrim)
			wmparams = np.loadtxt(wmfiletrim,ndmin=2)
			csfparams = np.loadtxt(csffiletrim,ndmin=2)
			if os.path.exists(scrubouttrim): 
				scrubparams = np.loadtxt(scrubouttrim,ndmin=2)
				allconfounds = np.hstack([mcparams,wmparams,csfparams,scrubparams])
				numconfounds = allconfounds.shape[1]
			else:
				allconfounds = np.hstack([mcparams,wmparams,csfparams])
				numconfounds = allconfounds.shape[1]

			numtimepoints = len(mcparams)
			heights = np.max(allconfounds,0)-np.min(allconfounds,0)
			print yellow + "Writing out trimmed deisgn file. Number of total regressors: %s%s" %(numconfounds,mainColor)

			df = open(designfile,'w')
			# df.write('/NumWaves\t%d\n' % numconfounds)
			# df.write('/NumPoints\t%s\n' % numtimepoints)
			# df.write('/PPheights\t')

			# for x in range(numconfounds):
			# 	df.write('%f\t' % heights[x])
			# df.write('\n\n/Matrix\n')
			# numtimepoints = int(numtimepoints)

			for row in range(numtimepoints):
				for column in range(numconfounds):
					df.write('%f\t' % allconfounds[row][column])
				df.write('\n')
			df.close()
		else: 
			print sectionColor2 + ' %s already exists. Moving on' %(designfile)


def checkImageLength(imagename):
	command = 'fslinfo %s' % imagename
	results = check_output(command,shell=True)
	TR = results.split()[9]
	return int(TR)

#Preprocessing steps
for subject in subjects:

	logfile = open(logfilename,'ab+')
	logfile.write("\n-----Analysis of subject %s started at %s\n" % (subject,starttime.strftime('%b %d %G %I:%M%p')))
	print sectionColor + '-----Analysis of subject %s started at %s\n%s' %(subject,starttime.strftime('%b %d %G %I:%M%p'),mainColor)
	if not args.nopre:

		#rename export folder to DICOM
		subjectfolder = "%s/%s/"  %(analysispath,subject)
		subjectnifti = "%s/%s/"  %(dicompath,subject)
		subjectdicom = dicompath + '/dicomdir/' + subject
		mpragebrain = subjectfolder + "mprage_brain.nii.gz"
		mprage = subjectfolder + "/mprage.nii.gz"
		t1 = subjectnifti + "anat/%s_T1w.nii.gz" %(subject)
		restfolder = subjectfolder + 'rest/'
		funcfolder = subjectfolder + 'music/'
		#slicetimefile = restfolder + "slicetime.txt"

		#create analysis folders: 
		if not os.path.exists(subjectfolder):
			print sectionColor + "Making analysis directory for %s...%s" %(subject,mainColor)
			os.mkdir(subjectfolder)

		#create restingstate and func folders
		if not os.path.exists(restfolder):
			print sectionColor + "Making resting state directory for %s...%s" %(subject,mainColor)
			os.mkdir(restfolder)

		if not os.path.exists(funcfolder):
			print sectionColor + "Making music task directory for %s...%s" %(subject,mainColor)
			os.mkdir(funcfolder)


		#convert from DICOM to NIFTI and put in neurobids in fmri_data
		if not args.nodcm:
			if not os.path.exists(subjectnifti):
				print sectionColor + "Converting %s's files to NIFTI%s"  % (subject,mainColor)
				#command = "dcm2niix -o %s -f '%%p' -z y %s" %(subjectfolder,dicomfolder) 
				command = "dcm2bids -d %s -p %s -c %s -o %s" %(subjectdicom,subject,configfile,dicompath)
				
				call(command,shell = True)
			else:
				print sectionColor2 + "Already converted %s files to NIFTI, moving on %s"  % (subject,mainColor)
		
		#Copy mprage
		if not os.path.exists(mprage):
			shutil.copy2(t1,mprage)
		
		#Skull stripping
		if not args.noskullstrip:
			if not os.path.exists(mpragebrain):
				command = '%s %s' % (skullstripscript,subjectfolder)
				call(command,shell=True)
			else:
				print sectionColor2 + "Already skullstripped, moving on %s"  % (mainColor)

		#Segmentation 
		if not args.noseg: 
			doseg(mpragebrain,subjectfolder)
		
		# Field maps preprocessing
		if not args.nofieldmap:	
			fmap = funcfolder + "/fieldmap_phase_rad.nii.gz"
			fmap_rest = restfolder + "fieldmap_phase_rad.nii.gz"

			
			#rename images
			if not os.path.exists(fmap_rest):
				print sectionColor + "Field map conversion for rest... %s"  %(mainColor)

				#Resting state field maps:
				phasesrc_rest = '%sfmap/%s_phase_run-01_fmap.nii.gz'  %(subjectnifti,subject)

				if os.path.exists(phasesrc_rest):
					phasedst_rest = '%sfieldmap_phase.nii.gz' %(restfolder)
					magsrc_rest = '%sfmap/%s_magnitude2_run-01_fmap.nii.gz'  %(subjectnifti,subject)
					magdst_rest = '%sfieldmap_mag.nii.gz' %(restfolder)

					shutil.copy2(phasesrc_rest,phasedst_rest)
					shutil.copy2(magsrc_rest,magdst_rest)

				#skull strip mag image
				command = "bet %sfieldmap_mag %sfieldmap_mag_brain" %(restfolder,restfolder)
				call(command,shell=True)

				#convert phase to rads
				command = "fsl_prepare_fieldmap SIEMENS %sfieldmap_phase.nii.gz %sfieldmap_mag_brain.nii.gz %sfieldmap_phase_rad 2.46" %(restfolder,restfolder,restfolder)
				call(command,shell=True)

			else:
				print sectionColor2 + "Field map conversion file already exist for rest moving on %s"  %(mainColor)

			if not os.path.exists(fmap):
				#Task-based field maps (always run-02)
				phasesrc_task = '%sfmap/%s_phase_run-02_fmap.nii.gz'  %(subjectnifti,subject)

				if os.path.exists(phasesrc_task):
					phasedst_task = '%sfieldmap_phase.nii.gz' %(funcfolder)
					magsrc_task = '%sfmap/%s_magnitude2_run-02_fmap.nii.gz'  %(subjectnifti,subject)
					magdst_task = '%sfieldmap_mag.nii.gz' %(funcfolder)

					shutil.copy2(phasesrc_task,phasedst_task)
					shutil.copy2(magsrc_task,magdst_task)

				#skull strip mag image
				command = "bet %sfieldmap_mag %sfieldmap_mag_brain" %(funcfolder,funcfolder)
				call(command,shell=True)

				#convert phase to rads
				command = "fsl_prepare_fieldmap SIEMENS %sfieldmap_phase.nii.gz %sfieldmap_mag_brain.nii.gz %sfieldmap_phase_rad 2.46" %(funcfolder,funcfolder,funcfolder)
				call(command,shell=True)

			else:
				print sectionColor2 + "Field map conversion files already exist for music, moving on %s"  %(mainColor)


		#Remove any non completed rest or task 
		funcpath = subjectnifti + 'func/'
		funcfiles = [elem for elem in os.listdir(funcpath) if ".nii.gz" in elem]
		funcfiles = [elem for elem in funcfiles if "run" in elem]
		for f in funcfiles: 
			funcfile = funcpath + f 
			length = (checkImageLength(funcfile))
			if length < 100:
				print red + 'File %s has %d TRs. Deleting...%s' %(f,length,mainColor) 
				os.remove(funcfile)
				os.remove(funcfile[:-6] + 'json')

		funcfiles = [elem for elem in os.listdir(funcpath) if "run" in elem]
		if len(funcfiles) == 2:
			for item in funcfiles:  
				funcoldname = funcpath + item
				funcnewname = funcoldname.replace(funcoldname[83:90],'') #This part doesn't work for rest, makes it a weird different name (different number of characters)
				print yellow + 'Renaming %s to %s%s' %(funcoldname,funcnewname,mainColor) 
				os.rename(funcoldname,funcnewname)

		elif len(funcfiles) > 2:
			print red + 'Still too many functional runs. Check file%s' %(mainColor)
			sys.exit()

		# #get slice timing files from afni using slicetimefiles.py script: 
		if not args.noslicetime:
			command = '%s %s' % (slicetimescript,subjectfolder)
			call(command,shell=True)

	else:
		print sectionColor2 + "Skipping preprocessing steps...%s" %(mainColor)


	#######################################################################
	## RUN THE FUNCTIONS
	#######################################################################

	if not args.nomusic:

		#Preprocessing
		if not args.nofeat:
			print sectionColor + "-----Preprocessing feat%s" %(mainColor)
			#modelname = "generic_pre_%s_200hpf" %song
			#dolowerlevels(modelname) #no registration, copy from resting state
			dolowerlevels("generic_pre_music_200hpf")

		if not args.noaroma:
			#song = songs[2-1]
			featfile = '%s%s_model1_pre_200hpf.feat/' % (funcfolder,song)
			icafolder = featfile + "/ICA_AROMA"
			icafiles = icafolder + '/melodic.ica'
			filteredfile = icafolder + "/denoised_func_data_nonaggr.nii.gz"

			if not os.path.exists(icafolder):
				#print red + "ICA-AROMA has not been completed for %s\n%s"  % (subject,mainColor)
				command = "%s -feat %s -out %s" % (icascript,featfile,icafolder)
				call(command,shell=True)
			else:
				print sectionColor2 + "ICA-AROMA has been completed for %s\n%s"  % (icafolder,mainColor)

			newoutput =  "%sfiltered_func_data_aroma.nii.gz" % (featfile)
			if not os.path.exists(newoutput):
				print sectionColor + "Copying denoised file to %s\n%s"  % (newoutput,mainColor)
				shutil.copy2(filteredfile,newoutput)

		#MExtract CSF and White matter
		if not args.noconfound:
			doconfound("music")

		#Cut off first 5 and last 10 
		if not args.notrim:
			dotrim("music")

		#Regress out confounds 
		#Full trim file (not cutting the first 20)
		trimdata = "%sfiltered_func_data_aroma_trim.nii.gz" % (featfile)
		residdata = "%sfiltered_func_data_200hpf_aroma_residuals.nii.gz" %(featfile)
		standata = "%sfiltered_func_data_200hpf_aroma_residuals_standard.nii.gz" % (featfile)
		designfile = "%sresid_design_trim_200hpf_aroma.txt" % (featfile)

		if not os.path.exists(trimdata):
			print red + 'Could not find %s. Moving on...%s' %(tffile, mainColor)
			continue

		if not os.path.exists(designfile):
			print red + 'Could not find %s. Moving on...%s' %(designfile, mainColor)
			continue

		if not os.path.exists(residdata):
			command = "fsl_glm -i %s --demean -m %smask.nii.gz -d %s --out_res=%s" % (trimdata,featfile,designfile,residdata)
			print sectionColor + 'Regressing out confounds for %s%s' %(trimdata, mainColor)
			call(command,shell=True)
		else:
			print sectionColor2 + 'Regressing out confounds already completed for %s%s' %(subject, mainColor)
		

		#Convert to standard space 
		restfile = "%sresting_pre.feat/" %(restfolder)
		standard = "%sreg/standard.nii.gz" % restfile
		premat = "%sreg/example_func2highres.mat" % (featfile)
		warpfile = "%sreg/highres2standard_warp.nii.gz" % (restfile)
		if not os.path.exists(standata):
			command = "applywarp --ref=%s --in=%s --out=%s --warp=%s --premat=%s" % (standard,residdata,standata,warpfile,premat)
			print sectionColor + 'Warping to standard space. Output is %s%s' %(standata, mainColor)
			call(command,shell=True)
		else:
			print sectionColor2 + 'Already converted to standard for %s%s' %(subject, mainColor)
		
		print 

	else:
		print sectionColor2 + "Skipping feats...%s" %(mainColor)
		
	endtime = datetime.now()
	delta = endtime - starttime
	logfile.write("-----Analysis of subject %s finished at %s, duration %s\n" % (subject,endtime.strftime('%b %d %G %I:%M%p'),str(delta)))
	logfile.close()

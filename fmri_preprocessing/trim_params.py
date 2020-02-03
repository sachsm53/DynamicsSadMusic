#!/usr/bin/env python

import sys,os
from subprocess import call, check_output
import argparse
from datetime import datetime
import re
#from scipy.interpolate import interp1d
#from scipy.interpolate import CubicSpline
#from scipy.signal import resample
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


#runs = [0,1,2,3]
runs = [2]
songs = ["rest","happy","sadln","sadsh"]
songrates = ['hnl','snl_s','snl_l']
song_timepoints =[360,183,515,271] #515 not 530 because already cut 10 from beginning and 5 from end
rate_time =[168,256,515]
conditions = ['emo','joy']

#command line options
parser = argparse.ArgumentParser()

parser.add_argument("--subjects",help="process listed subjects",nargs='+',action="store")
parser.add_argument("--all",help="process all subjects", action="store_true")
parser.add_argument("--half",help="process half subjects", action="store_true")
parser.add_argument("--nodata",help="skip trimming data", action="store_true")
parser.add_argument("--noresid",help="skip trimming residual params", action="store_true") #WM and CSF
parser.add_argument("--norating",help="skip trimming rating files", action="store_true") #WM and CSF
args = parser.parse_args()

def checkImageLength(imagename):
	command = 'fslinfo %s' % imagename
	results = check_output(command,shell=True)
	TR = results.split()[9]
	return int(TR)

#set paths
pathbase = "/Volumes/MusicProject/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-1/NaturalisticNetwork"
if not os.path.exists(pathbase):
	pathbase = "/Volumes/MusicProject-2/NaturalisticNetwork"

datapath = pathbase + "/fmri_analysis"	
dicompath = pathbase + "/fmri_data"

#Define number to cut from beginning
startvol = 20 ##CHANGE THIS HERE
cutvols = 20 ##CHANGE THIS HERE

#develop list of subjects
subjects = args.subjects

if args.all:
	#Get list of subjects
	subjects = os.listdir(datapath)
	subjects = [elem for elem in subjects if "sub" in elem]
	subjects = [elem for elem in subjects if "07" not in elem]
	subjects = [elem for elem in subjects if "pil" not in elem]
	subjects.sort()

if args.half:
	#Get list of subjects
	subjects = os.listdir(datapath)
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



if not args.nodata:

	if args. aroma:
		for subject in subjects:
			for run in runs:
				song = songs[run]
				subjectfolder = "%s/%s/music/%s_model1_pre_200hpf.feat" % (datapath,subject,song)
				#infile = "%s/filtered_func_data_200hpf_aroma_residuals_standard.nii.gz" % (subjectfolder)
				#outfile = "%s/filtered_func_data_200hpf_aroma_residuals_standard_cut20.nii.gz" % (subjectfolder)

				infile = "%s/filtered_func_data_200hpf_aroma_residuals.nii.gz" % (subjectfolder)
				outfile = "%s/filtered_func_data_200hpf_aroma_residuals_cut20.nii.gz" % (subjectfolder)
				addfile_new = "%s/filtered_func_data_200hpf_aroma_residuals_cut20_add.nii.gz" % (subjectfolder)
				mask = '%s/mask.nii.gz' %subjectfolder

				newdest = '%s/funpsy/data/%s_filtered_func_200hpf_aroma_residuals_standard_cut20.nii.gz' %(pathbase,subject)
				dstunzip = '%s/funpsy/data/%s_filtered_func_200hpf_aroma_residuals_standard_cut20.nii' %(pathbase,subject)

				if not os.path.exists(outfile): 
					totalvols  = checkImageLength(infile)
					trimvols = totalvols - cutvols
					command = "fslroi %s %s %s %s" % (infile,outfile,startvol,trimvols)
					print command
					call(command,shell=True)
				else:
					print red + '%s already exists. Moving on%s' %(outfile,mainColor)
					if not os.path.exists(addfile_new):
						additive = 10000
						command = "fslmaths %s -add %d -mul %s %s -odt float" % (outfile,additive,mask,addfile_new)
						print command
						call(command, shell = True)
					
					if not os.path.exists(newdest):
						shutil.copy2(outfile,newdest)

					if not os.path.exists(dstunzip):
						print 'Unzipping %s' %(newdest)
						command = 'gunzip %s' %(newdest)
						call(command,shell=True)
						continue
					else:
						print dstunzip
						continue

	else:
		for subject in subjects:
			for run in runs:
				song = songs[run]
				subjectfolder = "%s/%s/music/%s_model1_pre_200hpf.feat" % (datapath,subject,song)
				infile = "%s/filtered_func_data_trim.nii.gz" % (subjectfolder)
				outfile = "%s/filtered_func_data_cut20.nii.gz" % (subjectfolder)
				
				if not os.path.exists(outfile):
					#10 seconds have already been trimmed in preprocessing.
					#Now we want to remove an additional 5 from the end. 
					totalvols  = checkImageLength(infile)
					trimvols = totalvols - cutvols
					command = "fslroi %s %s %s %s" % (infile,outfile,startvol,trimvols)
					print command
					call(command,shell=True)
				else:
					print red + '%s already exists. Moving on%s' %(outfile,mainColor)

if not args.noresid:
	for subject in subjects:
		for run in runs:
			song = songs[run]

			subjectfolder = "%s/%s/music/%s_model1_pre_200hpf.feat/" % (datapath,subject,song)
			totalvols = song_timepoints[run] - cutvols

			#Trim the movement files
			mcfile = "%smc/prefiltered_func_data_mcf_trimmed.par" % (subjectfolder)
			mcfiletrim = "%smc/prefiltered_func_data_mcf_cut20.par" % (subjectfolder)

			if not os.path.exists(mcfiletrim):
				with open(mcfile,'r') as mc:
					lines = mc.readlines()
				trimmed_lines = lines[startvol:startvol+totalvols]
				with open(mcfiletrim,'w') as nf:
					nf.write(''.join(trimmed_lines))
			else: 
				print red + '%s already exists. Moving on' %(mcfiletrim)

			#Now trim the white matter and CSF files
			wmfile = "%swm_200hpf_trim.txt" %(subjectfolder)
			wmfiletrim = "%swm_200hpf_cut20.txt" % (subjectfolder)
			csffile = "%scsf_200hpf_trim.txt" % (subjectfolder)
			csffiletrim = "%scsf_200hpf_cut20.txt" % (subjectfolder)
			
			if not os.path.exists(wmfiletrim):
				with open(wmfile,'r') as wm:
					lines = wm.readlines()
				totalvols = song_timepoints[run] - cutvols
				trimmed_lines = lines[startvol:startvol+totalvols]
				with open(wmfiletrim,'w') as nf:
					nf.write(''.join(trimmed_lines))
			else: 
				print red + '%s already exists. Moving on' %(wmfiletrim)

			if not os.path.exists(csffiletrim):
				with open(csffile ,'r') as csf:
					lines = csf.readlines()
				totalvols = song_timepoints[run] - cutvols
				trimmed_lines = lines[startvol:startvol+totalvols]
				with open(csffiletrim,'w') as nf:
					nf.write(''.join(trimmed_lines))
			else: 
				print red + '%s already exists. Moving on' %(csffiletrim)


			#Now trim scrub file if it exists 
			scrubouttrim =  "%s/%s/music/scrub_fd_%s_trim.txt" %(datapath,subject,song)
			scrubcutout = "%sscrub_fd_%s_cut20.txt" %(subjectfolder,song)
			if os.path.exists(scrubouttrim):
				scrub = pd.read_table(scrubouttrim,header=None,delim_whitespace=True)
				print yellow + 'Found scrub file. There are %d slices to be regressed out' %(scrub.shape[1]) #scrub.count()[0]

				if not os.path.exists(scrubcutout):
					print yellow + "Trimming scrub file for %s%s"  % (song,mainColor)
					with open(scrubouttrim ,'r') as scrub:
						lines = scrub.readlines()
					trimmed_lines = lines[startvol:startvol+totalvols]
					with open(scrubcutout,'w') as snf:
						snf.write(''.join(trimmed_lines))
				else: 
					print red + '%s already exists. Moving on%s' %(scrubcutout,mainColor)
			else: 
				print yellow + 'No scrubbing need for %s. Moving on%s' %(subject,mainColor)

			####
			##Create design file	
			####
			designfile = "%s/resid_design_200hpf_cut20.txt" %(subjectfolder)

			#Construct design.mat
			if not os.path.exists(designfile):
				mcparams = np.loadtxt(mcfiletrim)
				wmparams = np.loadtxt(wmfiletrim,ndmin=2)
				csfparams = np.loadtxt(csffiletrim,ndmin=2)

				if os.path.exists(scrubcutout):
					scrubparams = np.loadtxt(scrubcutout,ndmin=2)
					scrub = pd.read_table(scrubouttrim,header=None,delim_whitespace=True)
					sums = scrub.sum(axis = 0)
					print yellow + 'Found scrub file. There are %d slices to be regressed out' %(sum(sums))

					allconfounds = np.hstack([mcparams,wmparams,csfparams,scrubparams])
					numconfounds = allconfounds.shape[1]
				else:
					allconfounds = np.hstack([mcparams,wmparams,csfparams])
					numconfounds = allconfounds.shape[1]

				numtimepoints = len(mcparams)
				heights = np.max(allconfounds,0)-np.min(allconfounds,0)
				print yellow + "Writing out trimmed design file. Number of total regressors: %s%s" %(numconfounds,mainColor)

				df = open(designfile,'w')

				for row in range(numtimepoints):
					for column in range(numconfounds):
						df.write('%f\t' % allconfounds[row][column])
					df.write('\n')
				df.close()
			else: 
				print sectionColor2 + '%s already exists. Moving on' %(designfile)

if not args.norating:
	

	for subject in subjects:
		for run in runs:
			#songrun = run
			#song = songrates[songrun]
			song = songs[run]

			subjectfolder = "%s/%s/music/%s_model1_pre_200hpf.feat/" % (datapath,subject,song)

			for c in conditions:
				ratefile = "%s/%s/music/rating_ev/%s_%s_ev_dmean.txt" % (datapath,subject,song,c)	
				#ratefile = "%s/%s/music/rating_ev/%s_min_ev.txt" % (datapath,subject,song)	
				outfile = '%s/%s_%s_ev_dmean_cut20.txt' %(subjectfolder,song,c)
				#outfile	= "%s/%s/music/rating_ev/%s_min_ev_cut20.txt" % (datapath,subject,song)	
				if not os.path.exists(ratefile):
					print red + "%s does not exist. Go check%s" %(ratefile,mainColor)
					continue


				if not os.path.exists(outfile):
					print yellow + 'Cutting ev file for song %s condition %s for %s' %(song,c,subject)
					with open(ratefile,'r') as r:
						lines = r.readlines()
					totalvols = song_timepoints[run] - cutvols
					trimmed_lines = lines[startvol:startvol+totalvols]
					with open(outfile,'w') as nf:
						len(outfile)
						nf.write(''.join(trimmed_lines))
				else: 
					print red + '%s already exists. Moving on%s' %(outfile,mainColor)



####################################

#CHECK TO SEE IF THE CUT REMOVED THE SCRUBBED TRs

####################################
# for subject in subjects:
# 	for run in runs:
# 		song = songs[run]

# 		subjectfolder = "%s/%s/music/%s_model1_pre.feat/" % (datapath,subject,song)
# 		scrubouttrim =  "%s/%s/music/scrub_fd_%s_trim.txt" %(datapath,subject,song)
# 		scrubcutout = "%sscrub_fd_%s_cut20.txt" %(subjectfolder,song)
		
# 		if os.path.exists(scrubcutout):
# 			scrub = pd.read_table(scrubcutout,header=None,delim_whitespace=True)
# 			sums = scrub.sum(axis = 0)
# 			if sum(sums) != scrub.shape[1]:
# 				print red + 'Found scrub file for %s. Originally there were %d slices to be regressed out, now there are %d' %(subject,scrub.shape[1],sum(sums)) #scrub.count()[0]
# 			else: 
# 				print yellow + 'Found scrub file for %s. There are %d slices to be regressed out' %(subject,sum(sums)) #scrub.count()[0]






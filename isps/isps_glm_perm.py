#!/usr/bin/env python

import sys,os
from subprocess import call, check_output
import argparse
from datetime import datetime
import re
import shutil
import numpy as np 
import pandas as pd
import brainiak.isfc
from brainiak import image, io
from brainiak.fcma.util import compute_correlation
import brainiak.utils.fmrisim
import matplotlib.pyplot as plt
import scipy.io
import statsmodels.api as sm
import statsmodels.stats.multitest
from scipy import stats
from nipype.algorithms import modelgen
from sklearn import preprocessing

#logging colors
sectionColor = "\033[94m"
sectionColor2 = "\033[96m"
groupColor = "\033[90m"
mainColor = "\033[92m"

#command line options
parser = argparse.ArgumentParser()
parser.add_argument("--song",help="select which song you want to run",nargs='+',action="store")
parser.add_argument("--reg",help="select which regressor to use",nargs='+',action="store")
args = parser.parse_args()

#Set path
basepath = '/Volumes/MusicProject/NaturalisticNetwork/funpsy'
if not os.path.exists(basepath):
	basepath = '/Volumes/MusicProject-1/NaturalisticNetwork/funpsy'

regpath = os.path.join(basepath,'regressors')

#Set up models
models = ['enjwrtemo','emowrtenj']
#models = ['emo','enj']
numperm = 5000

# DEFINITIONS
def prep_reg(song,reg):

	global v1
	global v2
	global v3
	global v4

	emopath = regpath + '/%s_%s_ratings_emo_cut.txt' %(song,reg)
	enjpath = regpath + '/%s_%s_ratings_enj_cut.txt' %(song,reg)

	emo_conv_path = regpath + '/%s_%s_ratings_emo_cut_hrfconv.txt' %(song,reg)
	enj_conv_path = regpath + '/%s_%s_ratings_enj_cut_hrfconv.txt' %(song,reg)

	emo_ortho_path = regpath + '/%s_%s_ratings_emo_cut_ortho_hrfconv.txt'%(song,reg)
	enj_ortho_path = regpath + '/%s_%s_ratings_enj_cut_ortho_hrfconv.txt'%(song,reg)

	#Load unpermuted original EV
	emo_ev = pd.read_csv(emopath,header = None)
	enj_ev = pd.read_csv(enjpath,header = None)


	#CONVOLVE WITH HRF
	if not os.path.exists(emo_conv_path):
		emo = emo_ev.values
		enj = enj_ev.values
		emo_conv = brainiak.utils.fmrisim.convolve_hrf(emo, 1, scale_function=True, temporal_resolution=1)
		enj_conv = brainiak.utils.fmrisim.convolve_hrf(enj, 1, scale_function=True, temporal_resolution=1)

		np.savetxt(emo_conv_path,emo_conv,fmt = '%0.4f')
		np.savetxt(enj_conv_path,enj_conv,fmt = '%0.4f')

	else:
		emo_conv = pd.read_csv(emo_conv_path,header = None)
		enj_conv = pd.read_csv(emo_conv_path,header = None)



	#ORTHOGONALIZE
	if not os.path.exists(emo_ortho_path):

		#Get in list form for ortho
		emo_list = []
		with open(emopath) as f:
			for line in f: 
				line = line.split()
				num = float(line[0])
				emo_list.append(num)

		enj_list = []
		with open(enjpath) as f:
			for line in f: 
				line = line.split()
				num = float(line[0])
				enj_list.append(num)

		#Orthogonalize both
		enjwrtemo = modelgen.orth(emo_list,enj_list) #orthogonalize enj wrt emo
		emowrtenj = modelgen.orth(enj_list,emo_list)
		enj_ortho = pd.DataFrame(enjwrtemo).values
		emo_ortho = pd.DataFrame(emowrtenj).values

		#Convolve with double gamma
		emo_ortho_conv = brainiak.utils.fmrisim.convolve_hrf(emo_ortho, 1, scale_function=True, temporal_resolution=1)
		enj_ortho_conv = brainiak.utils.fmrisim.convolve_hrf(enj_ortho, 1, scale_function=True, temporal_resolution=1)

		#Save
		np.savetxt(emo_ortho_path,emo_ortho_conv,fmt = '%0.4f\n')
		np.savetxt(enj_ortho_path,enj_ortho_conv,fmt = '%0.4f\n')
	
	else: 
		emo_ortho_conv = pd.read_csv(emo_ortho_path,header = None)
		enj_ortho_conv = pd.read_csv(enj_ortho_path,header = None)

	#Create dataframe
	# v1 = pd.DataFrame(stats.zscore(emo_conv))
	# v2 = pd.DataFrame(stats.zscore(enj_ortho_conv))
	# v3 = pd.DataFrame(stats.zscore(enj_conv))
	# v4 = pd.DataFrame(stats.zscore(emo_ortho_conv))

	v1 = pd.DataFrame(emo_conv)
	v2 = pd.DataFrame(enj_ortho_conv)
	v3 = pd.DataFrame(enj_conv)
	v4 = pd.DataFrame(emo_ortho_conv)

	#return v1,v2,v3,v4


def load_brain_data(song):
	global df
	global aal_labels
	print('Loading IPS data for', song,'...')
	# #LOAD IPSP data
	matdata = '%s/%s_ipsts_52rois.mat' %(basepath,song)
	ips = scipy.io.loadmat(matdata) 
	ipsts = ips['ipsts']
	df = pd.DataFrame(ipsts)

	#Load Label file
	aal_labels = pd.read_csv(os.path.join(basepath,'rois_aal_52.txt'),delimiter = '\t', header = None,names = ['count','region'])
	#return aal_labels



def run_glm(song,reg,m,nperm):
	global df_model1
	global df_model2

	if m == 'enjwrtemo':
		x = pd.concat([v1,v2],axis = 1)
		x.columns = ["emo_conv", "enj_ortho_conv"]
	elif m == 'emowrtenj':
		x = pd.concat([v3,v4],axis = 1)
		x.columns = ["enj_conv", "emo_ortho_conv"]

	all_coefs = []
	all_pvals = []
	all_coefs_ortho = []
	all_pvals_ortho = []
	#Run the GLM in python with both 
	for i in range(0,df.shape[1]):
		#print('Model %s. Running GLM for ROI Num %d R' %(m,i+1))
		y = df[[i]]
		x = sm.add_constant(x)
		model = sm.GLM(y,x).fit()
		#model.summary()
		#coef = model.params[[1]]
		coef = model.params
		all_coefs.append(coef.array[1])
		
		#If ortho, get r-value for that one
		all_coefs_ortho.append(coef.array[2])
		#print('%d: Rval is %.4f, Pval is %.4f' %(i+1,coef.array[0]))

		#Get a non-parametric beta values
		all_perm_coefs=[] 
		all_perm_coefs_ortho = []
		print('Running GLM for permutations...')
		for p in range(1,numperm+1):
			#Get in list form for ortho

			if reg == 'fused':
				emo_perm_fname = "%s/glm_perm/fft_%s_%s/emo_raw_ev%d.txt" %(basepath,song,reg,p)
				enj_perm_fname = "%s/glm_perm/fft_%s_%s/enj_raw_ev%d.txt" %(basepath,song,reg,p)
			elif reg == 'mean':
				emo_perm_fname = "%s/glm_perm/fft_%s_%s/emo_ev%d.txt" %(basepath,song,reg,p)
				enj_perm_fname = "%s/glm_perm/fft_%s_%s/enj_ev%d.txt" %(basepath,song,reg,p)
			
			emo_perm = []
			with open(emo_perm_fname) as f:
				for line in f: 
					line = line.split()
					num = float(line[0])
					emo_perm.append(num)

			enj_perm = []
			with open(enj_perm_fname) as f:
				for line in f: 
					line = line.split()
					num = float(line[0])
					enj_perm.append(num)

			#Orthogonalize both
			enjwrtemo_perm = modelgen.orth(emo_perm,enj_perm)
			emowrtenj_perm = modelgen.orth(enj_perm,emo_perm)
			enj_ortho = pd.DataFrame(enjwrtemo_perm).values
			emo_ortho = pd.DataFrame(emowrtenj_perm).values

			#Hemodynamic response function on all four 
			emo_perm_preconv = pd.DataFrame(emo_perm).values
			enj_perm_preconv = pd.DataFrame(enj_perm).values

			emo_perm_conv = brainiak.utils.fmrisim.convolve_hrf(emo_perm_preconv, 1, scale_function=True, temporal_resolution=1)
			enj_perm_conv = brainiak.utils.fmrisim.convolve_hrf(enj_perm_preconv, 1, scale_function=True, temporal_resolution=1)
			emo_perm_conv_ortho = brainiak.utils.fmrisim.convolve_hrf(emo_ortho, 1, scale_function=True, temporal_resolution=1)
			enj_perm_conv_ortho = brainiak.utils.fmrisim.convolve_hrf(enj_ortho, 1, scale_function=True, temporal_resolution=1)

			x1 = pd.DataFrame(emo_perm_conv)
			x2 = pd.DataFrame(enj_perm_conv_ortho)
			x3 = pd.DataFrame(enj_perm_conv)
			x4 = pd.DataFrame(emo_perm_conv_ortho)


			if m == 'enjwrtemo':
				xperm = pd.concat([x1,x2],axis = 1)
				xperm.columns = ["emo_conv", "enj_ortho_conv"]
			elif m == 'emowrtenj':
				xperm = pd.concat([x3,x4],axis = 1)
				xperm.columns = ["enj_conv", "emo_ortho_conv"]

			#Run GLM with perm
			xperm = sm.add_constant(xperm)
			pmodel = sm.GLM(y,xperm).fit()
			#est.summary()
			#model.summary()
			pcoef = pmodel.params
			all_perm_coefs.append(pcoef.array[1])
			
			#If ortho, get r-value for that one
			all_perm_coefs_ortho.append(pcoef.array[2])

		#Get Pvalue for each ROI 
		pval = (100 - stats.percentileofscore(all_perm_coefs,coef.array[1]))/100 
		all_pvals.append(pval) 
		print('Model: %s ROI %d %s: Rval is %.4f, Pval is %.4f' %(m,i+1,aal_labels.region[i],coef.array[1],pval))
		
		#If ortho, get p-value for that one
		pval_ortho = (100 - stats.percentileofscore(all_perm_coefs_ortho,coef.array[2]))/100 
		all_pvals_ortho.append(pval_ortho)
		#Print results	
		print('Model: %s ROI %d %s: Rval_ortho is %.4f, Pval_ortho is %.4f' %(m,i+1,aal_labels.region[i],coef.array[2],pval_ortho))

	## ALL ROIS
	# Put coefs in the right place
	if m == 'enjwrtemo':
		coefs_model1 = all_coefs
		coefs_model1_ortho = all_coefs_ortho
		pval_model1 = all_pvals
		pval_model1_ortho = all_pvals_ortho
		results1 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model1 = pd.DataFrame(results1)
		df_model1.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]
	elif m == 'emowrtenj':
		coefs_model2 = all_coefs
		pval_model2 = all_pvals
		coefs_model2_ortho = all_coefs_ortho
		pval_model2_ortho = all_pvals_ortho
		results2 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model2 = pd.DataFrame(results2)
		df_model2.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]




	#return coefs_model1,coefs_model2,pval_model1,pval_model2, coefs_model1_ortho,coefs_model2_ortho, pval_model1_ortho,pval_model2_ortho


# #ONCE BOTH MODELS HAVE BEEN RUN
def multi_correct(m):

	if m == 'enjwrtemo':
		coefs = df_model1.coefs
		coefs_ortho = df_model1.coefs_ortho
		pvals = df_model1.pvals
		pvals_ortho = df_model1.pvals_ortho
	elif m == 'emowrtenj':
		coefs = df_model2.coefs
		coefs_ortho = df_model2.coefs_ortho
		pvals = df_model2.pvals
		pvals_ortho = df_model2.pvals_ortho


	#MAIN IV
	pvals_fdr = statsmodels.stats.multitest.multipletests(pvals,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
	reject = pvals_fdr[0]
	adj_pvals = pvals_fdr[1]
	sigs = np.where(reject)[0] 

	if sigs.size > 0:
		print('Significant correlations for',m)
		for s in sigs:
			print('%d\t%s\t%.4f\t%.4f\t%.4f'%(s+1,aal_labels.region[s],coefs[s],adj_pvals[s],pvals[s]))
	else: 
		print('NOTHING SIGNIFICANT FOR',m)

	print()

	#ORTHO IV
	pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_ortho,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
	reject = pvals_fdr[0]
	adj_pvals = pvals_fdr[1]
	sigs = np.where(reject)[0] 

	if sigs.size > 0:
		print('Significant correlations for',m,'ORTHO')
		for s in sigs:
			print('%d\t%s\t%.4f\t%.4f\t%.4f'%(s+1,aal_labels.region[s],coefs_ortho[s],adj_pvals[s],pvals_ortho[s]))
	else: 
		print('NOTHING SIGNIFICANT FOR',m,'ORTHO')
	print()


	##All pvals together 
	pvals_all = df_model1.pvals + df_model2.pvals #enjoy
	coefs_all = df_model1.coefs + df_model2.coefs #enjoy

	#MAIN IV
	pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_all,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
	reject = pvals_fdr[0]
	adj_pvals = pvals_fdr[1]
	sigs = np.where(reject)[0] 

	if sigs.size > 0:
		print('Significant correlations when combining across 2 models')
		for s in sigs:
			print('%d\t%s\t%.4f\t%.4f\t%.4f'%(s+1,aal_labels.region[s],coefs_all[s],adj_pvals[s],pvals_all[s]))
	else: 
		print('NOTHING SIGNIFICANT when combining across 2 models')

	print()

	# #MAIN IV
	# pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_emo,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
	# reject = pvals_fdr[0]
	# adj_pvals = pvals_fdr[1]
	# sigs = np.where(reject)[0] 

	# if sigs.size > 0:
	# 	print('Significant correlations for EMOTION when combining across 2 models')
	# 	for s in sigs:
	# 		print('%d\t%s\t%.4f\t%.4f\t%.4f'%(s+1,aal_labels.region[s],coefs_emo[s],adj_pvals[s],pvals_emo[s]))
	# else: 
	# 	print('NOTHING SIGNIFICANT for EMOTION when combining across 2 models')

	# print()


def network_glm(song,reg,m,nperm):
	global df_model1
	global df_model2
	global df_model3
	global df_model4
	global df_model5
	global df_model6
	all_coefs = []
	all_pvals = []
	all_coefs_ortho = []
	all_pvals_ortho = []
	networks = ['auditory','striatal','limbic','orbito','dmn']
	#network_label = np.concatenate([np.repeat('aud',8),np.repeat('striat',8),np.repeat('limbic',12),np.repeat('orbito',12),np.repeat('dmn',12)])

	network_label = pd.read_csv(os.path.join(basepath,'networks_aal_5.txt'),delimiter = '\t', header = None,names = ['count','region','network'])

	if m == 'enjwrtemo':
		x = pd.concat([v1,v2],axis = 1)
		x.columns = ["emo_conv", "enj_ortho_conv"]
	elif m == 'emowrtenj':
		x = pd.concat([v3,v4],axis = 1)
		x.columns = ["enj_conv", "emo_ortho_conv"]
	elif m == 'rms_emo':
		x = pd.concat([v1,v7],axis = 1)
		x.columns = ["emo_conv", "rms_ortho_conv"]
	elif m == 'br_emo':
		x = pd.concat([v1,v8],axis = 1)
		x.columns = ["emo_conv", "br_ortho_conv"]
	elif m == 'rms_enj':
		x = pd.concat([v3,v9],axis = 1)
		x.columns = ["emo_conv", "rms_ortho_conv"]
	elif m == 'br_enj':
		x = pd.concat([v3,v10],axis = 1)
		x.columns = ["emo_conv", "br_ortho_conv"]
	elif m == 'rms':
		x = v5
	elif m == 'br':
		x = v6


	for i,n in enumerate(networks):
		netnodes = np.where(network_label.network == n)
		network = df[netnodes[0]]
		y = network.mean(axis = 1)
		
		x = sm.add_constant(x)
		model = sm.GLM(y,x).fit()
		coef = model.params
		all_coefs.append(coef.array[1])
		
		#If ortho, get r-value for that one
		all_coefs_ortho.append(coef.array[2])
		#print('%d: Rval is %.4f, Pval is %.4f' %(i+1,coef.array[0]))

		#Get a non-parametric beta values
		all_perm_coefs=[] 
		all_perm_coefs_ortho = []
		print('Running GLM for permutations...')
		for p in range(1,nperm+1):
			#Get in list form for ortho

			if reg == 'fused':
				emo_perm_fname = "%s/glm_perm/fft_%s_%s/emo_raw_ev%d.txt" %(basepath,song,reg,p)
				enj_perm_fname = "%s/glm_perm/fft_%s_%s/enj_raw_ev%d.txt" %(basepath,song,reg,p)
			elif reg == 'mean':
				emo_perm_fname = "%s/glm_perm/fft_%s_%s/emo_ev%d.txt" %(basepath,song,reg,p)
				enj_perm_fname = "%s/glm_perm/fft_%s_%s/enj_ev%d.txt" %(basepath,song,reg,p)
			
			emo_perm = []
			with open(emo_perm_fname) as f:
				for line in f: 
					line = line.split()
					num = float(line[0])
					emo_perm.append(num)

			enj_perm = []
			with open(enj_perm_fname) as f:
				for line in f: 
					line = line.split()
					num = float(line[0])
					enj_perm.append(num)

			#Orthogonalize both
			enjwrtemo_perm = modelgen.orth(emo_perm,enj_perm)
			emowrtenj_perm = modelgen.orth(enj_perm,emo_perm)
			enj_ortho = pd.DataFrame(enjwrtemo_perm).values
			emo_ortho = pd.DataFrame(emowrtenj_perm).values

			#Hemodynamic response function on all four 
			emo_perm_preconv = pd.DataFrame(emo_perm).values
			enj_perm_preconv = pd.DataFrame(enj_perm).values

			emo_perm_conv = brainiak.utils.fmrisim.convolve_hrf(emo_perm_preconv, 1, scale_function=True, temporal_resolution=1)
			enj_perm_conv = brainiak.utils.fmrisim.convolve_hrf(enj_perm_preconv, 1, scale_function=True, temporal_resolution=1)
			emo_perm_conv_ortho = brainiak.utils.fmrisim.convolve_hrf(emo_ortho, 1, scale_function=True, temporal_resolution=1)
			enj_perm_conv_ortho = brainiak.utils.fmrisim.convolve_hrf(enj_ortho, 1, scale_function=True, temporal_resolution=1)

			#GET RMS AND BR
			rms_ev_fname = "%s/glm_perm/fft_%s_acoustic/rms_ev%d.txt" %(basepath,song,p)
			rms_perm = []
			with open(rms_ev_fname) as f:
				for line in f: 
					line = line.split()
					num = float(line[0])
					rms_perm.append(num)

			#Orthogonalize both
			rmspermwrtenj = modelgen.orth(enj_perm,rms_perm) #orthogonalize enj wrt emo
			rms_perm_enj_ortho = pd.DataFrame(rmspermwrtenj).values
			rms_perm_enj_ortho_conv = brainiak.utils.fmrisim.convolve_hrf(rms_perm_enj_ortho, 1, scale_function=True, temporal_resolution=1)

			rmspermwrtemo = modelgen.orth(emo_perm,rms_perm) #orthogonalize enj wrt emo
			rms_perm_emo_ortho = pd.DataFrame(rmspermwrtenj).values
			rms_perm_emo_ortho_conv = brainiak.utils.fmrisim.convolve_hrf(rms_perm_emo_ortho, 1, scale_function=True, temporal_resolution=1)

			#Get brightness 
			br_ev_fname = "%s/glm_perm/fft_%s_acoustic/br_ev%d.txt" %(basepath,song,p)
			br_perm = []
			with open(br_ev_fname) as f:
				for line in f: 
					line = line.split()
					num = float(line[0])
					br_perm.append(num)

			#Orthogonalize both
			rmspermwrtenj = modelgen.orth(enj_list,rms_perm) #orthogonalize enj wrt emo
			rms_perm_ortho = pd.DataFrame(rmspermwrtenj).values
			rms_perm_ortho_conv = brainiak.utils.fmrisim.convolve_hrf(rms_perm_ortho, 1, scale_function=True, temporal_resolution=1)

			#Orthogonalize both
			brpermwrtenj = modelgen.orth(enj_perm,br_perm) #orthogonalize enj wrt emo
			br_perm_enj_ortho = pd.DataFrame(brpermwrtenj).values
			br_perm_enj_ortho_conv = brainiak.utils.fmrisim.convolve_hrf(br_perm_enj_ortho, 1, scale_function=True, temporal_resolution=1)

			brpermwrtemo = modelgen.orth(emo_perm,br_perm) #orthogonalize enj wrt emo
			br_perm_emo_ortho = pd.DataFrame(brpermwrtemo).values
			br_perm_emo_ortho_conv = brainiak.utils.fmrisim.convolve_hrf(br_perm_emo_ortho, 1, scale_function=True, temporal_resolution=1)

			#Hrmsdynamic response function on all four 
			rms_perm_preconv = pd.DataFrame(rms_perm).values
			rms_perm_conv = brainiak.utils.fmrisim.convolve_hrf(rms_perm_preconv, 1, scale_function=True, temporal_resolution=1)

			#Hrmsdynamic response function on all four 
			br_perm_preconv = pd.DataFrame(br_perm).values
			br_perm_conv = brainiak.utils.fmrisim.convolve_hrf(br_perm_preconv, 1, scale_function=True, temporal_resolution=1)


			x1 = pd.DataFrame(emo_perm_conv)
			x2 = pd.DataFrame(enj_perm_conv_ortho)
			x3 = pd.DataFrame(enj_perm_conv)
			x4 = pd.DataFrame(emo_perm_conv_ortho)
			x5 = pd.DataFrame(rms_perm_conv)
			x6 = pd.DataFrame(br_perm_conv)
			x7 = pd.DataFrame(rms_perm_emo_ortho_conv)
			x8 = pd.DataFrame(br_perm_emo_ortho_conv)
			x9 = pd.DataFrame(rms_perm_enj_ortho_conv)
			x10 = pd.DataFrame(br_perm_enj_ortho_conv)


			if m == 'enjwrtemo':
				xperm = pd.concat([x1,x2],axis = 1)
				xperm.columns = ["emo_conv", "enj_ortho_conv"]
			elif m == 'emowrtenj':
				xperm = pd.concat([x3,x4],axis = 1)
				xperm.columns = ["enj_conv", "emo_ortho_conv"]
			if m == 'rms_emo':
				xperm = pd.concat([x1,x7],axis = 1)
				xperm.columns = ["emo_conv", "enj_ortho_conv"]
			elif m == 'br_emo':
				xperm = pd.concat([x1,x8],axis = 1)
				xperm.columns = ["enj_conv", "emo_ortho_conv"]
			if m == 'rms_enj':
				xperm = pd.concat([x3,x9],axis = 1)
				xperm.columns = ["emo_conv", "enj_ortho_conv"]
			elif m == 'br_enj':
				xperm = pd.concat([x3,x10],axis = 1)
				xperm.columns = ["enj_conv", "emo_ortho_conv"]
			if m == 'rms':
				xperm = x5
			elif m == 'br':
				xperm = x6

			#Run GLM with perm
			xperm = sm.add_constant(xperm)
			pmodel = sm.GLM(y,xperm).fit()
			#est.summary()
			#model.summary()
			pcoef = pmodel.params
			all_perm_coefs.append(pcoef.array[1])
			
			#If ortho, get r-value for that one
			all_perm_coefs_ortho.append(pcoef.array[2])

		#Get Pvalue for each ROI 
		pval = (100 - stats.percentileofscore(all_perm_coefs,coef.array[1]))/100 
		all_pvals.append(pval) 
		print('Model: %s ROI %d %s: Rval is %.4f, Pval is %.4f' %(m,i+1,n,coef.array[1],pval))
		
		#If ortho, get p-value for that one
		pval_ortho = (100 - stats.percentileofscore(all_perm_coefs_ortho,coef.array[2]))/100 
		all_pvals_ortho.append(pval_ortho)
		
		#Print results	
		print('Model: %s ROI %d %s: Rval_ortho is %.4f, Pval_ortho is %.4f' %(m,i+1,n,coef.array[2],pval_ortho))

		## ALL ROIS
	# Put coefs in the right place
	if m == 'enjwrtemo':
		coefs_model1 = all_coefs
		coefs_model1_ortho = all_coefs_ortho
		pval_model1 = all_pvals
		pval_model1_ortho = all_pvals_ortho
		results1 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model1 = pd.DataFrame(results1)
		df_model1.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]
	elif m == 'emowrtenj':
		coefs_model2 = all_coefs
		pval_model2 = all_pvals
		coefs_model2_ortho = all_coefs_ortho
		pval_model2_ortho = all_pvals_ortho
		results2 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model2 = pd.DataFrame(results2)
		df_model2.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]
	elif m == 'rms_emo':
		pval_model3_ortho = all_pvals_ortho
		results3 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model3 = pd.DataFrame(results3)
		df_model3.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]
	elif m == 'br_emo':
		results4 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model4 = pd.DataFrame(results4)
		df_model4.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]
	elif m == 'rms_enj':
		results5 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model5 = pd.DataFrame(results5)
		df_model5.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]
	elif m == 'br_enj':
		results6 = np.column_stack([all_coefs,all_pvals,all_coefs_ortho,all_pvals_ortho])
		df_model6 = pd.DataFrame(results6)
		df_model6.columns = ["coefs", "pvals","coefs_ortho","pvals_ortho"]


#### RUN 
song = 'sadln'
reg = 'mean'
#prep_reg(song,reg)
load_brain_data(song)
models = ['rms_emo','br_emo','rms_enj','br_enj']
for n,m in enumerate(models):
	#run_glm(song,reg,m,5000)
	network_glm(song,reg,m,5000)

for n,m in enumerate(models):
	multi_correct(m)

##All pvals together 
pvals_all = np.concatenate([df_model1.pvals,df_model2.pvals]) #emo followed by enjoy
coefs_all = np.concatenate([df_model1.coefs,df_model2.coefs]) #emo followed by enjoy

pvals_all = np.concatenate([df_model3.pvals,df_model5.pvals]) #RMS
coefs_all = np.concatenate([df_model3.coefs,df_model5.coefs]) #RMS

pvals_all = np.concatenate([df_model4.pvals,df_model6.pvals]) #BR
coefs_all = np.concatenate([df_model4.coefs,df_model6.coefs]) #BR

pvals_all = np.concatenate([df_model3.pvals,df_model4.pvals,df_model5.pvals,df_model6.pvals]) #ALL
coefs_all = np.concatenate([df_model3.coefs,df_model4.coefs,df_model6.coefs,df_model6.coefs]) #BR

pvals_all = np.concatenate([df_model_emo_br.pvals,df_model_enj_br.pvals]) #emo followed by enjoy
coefs_all = np.concatenate([df_model_emo_br.coefs,df_model_enj_br.coefs]) #emo followed by enjoy

#MAIN IV
pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_all,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
reject = pvals_fdr[0]
adj_pvals = pvals_fdr[1]
sigs = np.where(reject)[0] 

networks = ['auditory','striatal','limbic','orbito','dmn']
networks2 = networks*2
if sigs.size > 0:
	print('Significant correlations when combining across 2 models')
	for s in sigs:
		if s < len(networks2):
			model = 'emotion'
		else:
			model = 'enjoyment'

		print('%s\t%s\t%.4f\t%.4f\t%.4f'%(model,networks2[s],coefs_all[s],adj_pvals[s],pvals_all[s]))
else: 
	print('NOTHING SIGNIFICANT when combining across 2 models')

print()




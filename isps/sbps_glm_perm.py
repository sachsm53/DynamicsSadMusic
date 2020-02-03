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


def load_brain_data(song,kind): #type is either isps or sbps
	global df
	global aal_labels
	global labels
	global sbps
	print('Loading %s data for %s' %(kind,song))

	if kind == 'isps':
		# #LOAD IPSP data
		matdata = '%s/%s_ipsts_52rois.mat' %(basepath,song)
		ips = scipy.io.loadmat(matdata) 
		ipsts = ips['ipsts']
		df = pd.DataFrame(ipsts)

		#Load Label file
		aal_labels = pd.read_csv(os.path.join(basepath,'rois_aal_52.txt'),delimiter = '\t', header = None,names = ['count','region'])
		#return aal_labels

	elif kind == 'sbps':
		#LOAD IPSP data
		matdata = '%s/sadln_sbps_52rois.mat' %(basepath)
		sbpsmat = scipy.io.loadmat(matdata) 
		sbps = pd.DataFrame(sbpsmat['data'])
		labels = pd.DataFrame(sbpsmat['labels'],columns = ["roi1", "roi2"])
		labels['roi1'] = labels['roi1'].map(lambda x: str(x)[1:-1])
		labels['roi2'] = labels['roi2'].map(lambda x: str(x)[1:-1])

		#count = 0
		# for i,r in enumerate(labels.roi1):
		# 	#if roi1 == "'Frontal_Sup_Orb_L'" or roi2 == "'Frontal_Sup_Orb_L'" or roi1 == "'Frontal_Sup_Orb_R'" or roi2 == "'Frontal_Sup_Orb_R'":
		# 	if labels.roi1[i] in lofc and labels.roi2[i] not in lofc or labels.roi2[i] in lofc and labels.roi1[i] not in lofc:
		# 		count = count + 1 
		# 		print(count,i,r,labels.roi2[i]) #174, should be 192



def cross_network(song,reg,m,nperm,kind):
	global df_results
	global df_results_emo
	global df_results_enj


	networks = ['auditory','striatal','limbic','orbito','dmn']
	network_label = pd.read_csv(os.path.join(basepath,'networks_aal_5.txt'),delimiter = '\t', header = None,names = ['count','region','network'])
	network_names = pd.read_csv(os.path.join(basepath,'rois_aal_52.txt'),delimiter = '\t', header = None,names = ['count','region'])
	#network_names['network'] = network_label

	all_roi_names = ['Heschl_L', 'Heschl_R', 'Temporal_Sup_L', 'Temporal_Sup_R',
	   'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R',
	   'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Caudate_L',
	   'Caudate_R', 'Putamen_L', 'Putamen_R', 'Pallidum_L', 'Pallidum_R',
	   'Thalamus_L', 'Thalamus_R', 'Insula_L', 'Insula_R',
	   'Cingulum_Ant_L', 'Cingulum_Ant_R', 'Cingulum_Mid_L',
	   'Cingulum_Mid_R', 'Hippocampus_L', 'Hippocampus_R',
	   'ParaHippocampal_L', 'ParaHippocampal_R', 'Amygdala_L',
	   'Amygdala_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R',
	   'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Med_Orb_L',
	   'Frontal_Med_Orb_R', 'Rectus_L', 'Rectus_R', 'Frontal_Sup_L',
	   'Frontal_Sup_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R',
	   'Cingulum_Post_L', 'Cingulum_Post_R', 'Parietal_Inf_L',
	   'Parietal_Inf_R', 'Angular_L', 'Angular_R', 'Precuneus_L',
	   'Precuneus_R', 'Temporal_Mid_L', 'Temporal_Mid_R',
	   'Temporal_Inf_L', 'Temporal_Inf_R']

	all_roi_names = ["'" + item + "'" for item in all_roi_names]

	if m == 'enjwrtemo':
		x = pd.concat([v1,v2],axis = 1)
		x.columns = ["emo_conv", "enj_ortho_conv"]
	elif m == 'emowrtenj':
		x = pd.concat([v3,v4],axis = 1)
		x.columns = ["enj_conv", "emo_ortho_conv"]

    
	#Create empty arrays
	all_coefs_main = []
	all_pvals_main = []
	all_coefs_ortho = []
	all_pvals_ortho = []

	counter = 0
	for i,n in enumerate(networks):
        
		netnames = [all_roi_names[i] for i in np.where(network_label == n)[0]]
		inx_roi1_l = [i for i, j in enumerate(labels.roi1) if j in netnames]
		inx_roi1_r =  [i for i, j in enumerate(labels.roi2) if j in netnames]
		inx_roi1 = np.concatenate([inx_roi1_l,inx_roi1_r])
        
        
		if kind == 'within':
			netnodes = np.intersect1d(inx_roi1_l,inx_roi1_r)
			print()
			print('Running SBPS for %d %s network within: %d nodes' %(i,n,len(netnodes)))

			network = sbps[netnodes]
			y = network.mean(axis = 1)

			xmod = sm.add_constant(x)
			model = sm.GLM(y,xmod).fit()
			coef = model.params
			#Get r-value for that ROI for both enjoyment and emotion
			all_coefs_main.append(coef.array[1])
			all_coefs_ortho.append(coef.array[2])

			#Get a non-parametric beta values with permutation
			all_perm_coefs_main = [] 
			all_perm_coefs_ortho = []

			print('Running GLM for',numperm, 'permutations...')

            #Loop through permutations
			for p in range(1,numperm+1):
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


				xperm = sm.add_constant(xperm)

				#Run GLM with perm
				pmodel = sm.GLM(y,xperm).fit()
				pcoef = pmodel.params

				all_perm_coefs_main.append(pcoef.array[1])
				all_perm_coefs_ortho.append(pcoef.array[2])

			#Get Pvalue for each network for first regressor 
			pval_main = (100 - stats.percentileofscore(all_perm_coefs_main,coef.array[1]))/100 
			all_pvals_main.append(pval_main) 
			print('Main: Network %s: Rval is %.4f, Pval is %.4f' %(n,coef.array[1],pval_main))

			#Get Pvalue for each network for first second regressor 
			pval_ortho = (100 - stats.percentileofscore(all_perm_coefs_ortho,coef.array[2]))/100 
			all_pvals_ortho.append(pval_ortho)
			print('Ortho: Network %s: Rval is %.4f, Pval is %.4f' %(n,coef.array[2],pval_ortho))

			#Put into dataframe
			results = np.column_stack((all_coefs_main,all_pvals_main,all_coefs_ortho,all_pvals_ortho))
			df_results = pd.DataFrame(results)
			df_results.columns = ["main_coefs", "main_pvals","ortho_coefs","ortho_pvals"]


		elif kind == 'between':
			for i2,n2 in enumerate(networks):
				if i2 > i:
					counter = counter + 1
					netnames2 = [all_roi_names[i] for i in np.where(network_label == n2)[0]]
					inx_roi2_l = [i for i, j in enumerate(labels.roi1) if j in netnames2]
					inx_roi2_r =  [i for i, j in enumerate(labels.roi2) if j in netnames2]
					inx_roi2 = np.concatenate([inx_roi2_l,inx_roi2_r])

					netnodes = np.intersect1d(inx_roi1,inx_roi2)
					print()
					print('%d SBPS for %d %s network to %d %s network: %d nodes' %(counter,i,n,i2,n2,len(netnodes)))

					network = sbps[netnodes]
					y = network.mean(axis = 1)

					xmod = sm.add_constant(x)
					model = sm.GLM(y,xmod).fit()
					coef = model.params
					#Get r-value for that ROI for both enjoyment and emotion
					all_coefs_main.append(coef.array[1])
					all_coefs_ortho.append(coef.array[2])

					#Get a non-parametric beta values with permutation
					all_perm_coefs_main = [] 
					all_perm_coefs_ortho = []

					print('Running GLM for',numperm, 'permutations...')

                    #Loop through permutations
					for p in range(1,numperm+1):
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


						xperm = sm.add_constant(xperm)

						#Run GLM with perm
						pmodel = sm.GLM(y,xperm).fit()
						pcoef = pmodel.params

						all_perm_coefs_main.append(pcoef.array[1])
						all_perm_coefs_ortho.append(pcoef.array[2])

					#Get Pvalue for each network for first regressor 
					pval_main = (100 - stats.percentileofscore(all_perm_coefs_main,coef.array[1]))/100 
					all_pvals_main.append(pval_main) 
					print('Main: Network %s: Rval is %.4f, Pval is %.4f' %(n,coef.array[1],pval_main))

					#Get Pvalue for each network for first second regressor 
					pval_ortho = (100 - stats.percentileofscore(all_perm_coefs_ortho,coef.array[2]))/100 
					all_pvals_ortho.append(pval_ortho)
					print('Ortho: Network %s: Rval is %.4f, Pval is %.4f' %(n,coef.array[2],pval_ortho))

					#Put into dataframe
					results = np.column_stack((all_coefs_main,all_pvals_main,all_coefs_ortho,all_pvals_ortho))
					df_results = pd.DataFrame(results)
					df_results.columns = ["main_coefs", "main_pvals","ortho_coefs","ortho_pvals"]

			if m == 'enjwrtemo': 
				df_results_emo = df_results
			elif m == 'emowrtenj':
				df_results_enj = df_results


	if m == 'enjwrtemo': 
		df_results_emo = df_results
	elif m == 'emowrtenj':
		df_results_enj = df_results





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
	pvals_enj = df_model2.pvals + df_model1.pvals_ortho #enjoy
	coefs_enj = df_model2.coefs + df_model1.coefs_ortho #enjoy
	pvals_emo = df_model1.pvals + df_model2.pvals_ortho #emo
	coefs_enj = df_model1.coefs + df_model2.coefs_ortho #enjoy

	#MAIN IV
	pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_enj,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
	reject = pvals_fdr[0]
	adj_pvals = pvals_fdr[1]
	sigs = np.where(reject)[0] 

	if sigs.size > 0:
		print('Significant correlations for ENJOY when combining across 2 models')
		for s in sigs:
			print('%d\t%s\t%.4f\t%.4f\t%.4f'%(s+1,aal_labels.region[s],coefs_enj[s],adj_pvals[s],pvals_enj[s]))
	else: 
		print('NOTHING SIGNIFICANT for ENJOY when combining across 2 models')

	print()

	#MAIN IV
	pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_emo,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
	reject = pvals_fdr[0]
	adj_pvals = pvals_fdr[1]
	sigs = np.where(reject)[0] 

	if sigs.size > 0:
		print('Significant correlations for EMOTION when combining across 2 models')
		for s in sigs:
			print('%d\t%s\t%.4f\t%.4f\t%.4f'%(s+1,aal_labels.region[s],coefs_emo[s],adj_pvals[s],pvals_emo[s]))
	else: 
		print('NOTHING SIGNIFICANT for EMOTION when combining across 2 models')

	print()




#### RUN 
song = 'sadln'
reg = 'mean'
prep_reg(song,reg)
load_brain_data(song,'sbps')
models = ['emowrtenj','enjwrtemo']
for n,m in enumerate(models):
	#run_glm(song,reg,m,5000)
	#network_glm(song,reg,m,5000)
	cross_network(song,reg,m,5000,'between')

# for n,m in enumerate(models):
# 	multi_correct(m)

#For network
networks_all = np.append(networks,networks)
pvals_all = np.concatenate([df_results_emo.main_pvals, df_results_enj.main_pvals]) 
coefs_all = np.concatenate([df_results_emo.main_coefs, df_results_enj.main_coefs])
#MAIN IV
pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_all,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
reject = pvals_fdr[0]
adj_pvals = pvals_fdr[1]
sigs = np.where(reject)[0] 

if sigs.size > 0:
	print('Significant correlations for ENJOY when combining across 2 models')
	for s in sigs:
		print('%d\t%s\t%.4f\t%.4f\t%.4f'%(s+1,networks_all[s],coefs_all[s],adj_pvals[s],pvals_all[s]))
else: 
	print('NOTHING SIGNIFICANT for ENJOY when combining across 2 models')


#FDR Correction, combining across both models
pvals_combine = np.concatenate([df_results_emo.main_pvals,df_results_enj.main_pvals])
coefs_combine = np.concatenate([df_results_emo.main_coefs,df_results_enj.main_coefs])
networks = ['auditory','striatal','limbic','orbito','dmn']
networks_combine = networks*2
pvals_fdr = statsmodels.stats.multitest.multipletests(pvals_combine,alpha = 0.05,method = 'fdr_bh',is_sorted=False)
reject = pvals_fdr[0]
adj_pvals = pvals_fdr[1]
sigs = np.where(reject)[0]
if sigs.size > 0:
	if s < len(networks):
		model = 'emotion'
	else:
		model = 'enjoyment'
	print('Significant correlations across two models')
	for s in sigs:
		print('%s\t%s\t%.4f\t%.4f\t%.4f'%(model,networks_combine[s],pvals_combine[s],adj_pvals[s],coefs_combine[s]))
else: 
	print('NOTHING SIGNIFICANT FOR',df_results.columns[i])





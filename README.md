# Scripts and analysis steps for *Dynamic intersubject synchronization in response to sad music*

All scripts used to analyze and create images for data presented in Sachs, M. E., Habibi, A., Damasio, A., & Kaplan, J. T. (2019). Dynamic intersubject neural synchronization reflects affective responses to sad music. NeuroImage, 116512.


### Prerequisites

Programs you will need downloaded: 

* FSL
* dcm2bids
* ICA-AROMA
* ISCToolbox
* MIRToolbox
* AFNI

The following python packages: 

* brainiak
* nibabel
* nipype

Other files/scripts that are called: 

`generic_resting_pre.fsf`
`generic_stats_design_enj_wrt_emo.fsf`
`generic_stats_design_emo_wrt_enj.fsf`
`slicetimefiles.py`
`ICA_AROMA.py`
`skullstrip.py`
`ICA_AROMA.py`
`copyreg.sh`

### Stimuli

The three normalized, .mp3 files used during fMRI. For this publication, only the piece entitled "snl_l.norm" which is entitled "Discovery of the Camps" from the *Band of Brothers* soundtrack

### Preprocessing Scripts and steps (in folder fmri_preprocessing)

  1) `scrubbing.py` gets the scrubbed regressors using FD greater than 0.5

  2) `preprocessing_snl_200hpf_aroma.py` performs all preprocessing steps on raw BIDS data. The steps, in the following order, include: 
  		- converts from DICOM to NIFTI using dcm2bids
  		- skullstripping using `skullstrip.py` that calls FSL (if reg folder already exists, it copies using copyreg.sh it rather than redoing this)
  		- doseg: performs segmentation on MPRAGE using FSL FAST (if reg folder already exists, it copies using copyreg.sh it rather than redoing this)
  		- Field map conversion using FSL 
  		- checkImageLength: function to makes sure image is correct length
  		- calls `slicetimefiles.py` to get slice timing files using afni code (dicom_hdr -slice_times)
	  	- dolowerlevel using `generic_resting_pre.fsf`: 
			- smoothing (5mm)
			- MCFLIRT
			- no slice time correction
			- no temporal filtering filtering (neither high nor low)
			- deletes first 5 volumes (silence)
			- field map unwarping
			- functional skull stripping
			- motion correction (MCFLIRT)
			- registration, BBR to the highres, then linear + nonlinear to the MNI space
				- OR copies reg folder from the rest preprocessed folder, if already completed
		- docheckfeat: checks the length of the output of FEAT (filtered_func_data.nii.gz), preprocessed data, to make sure it ran correctly
			- output file is filtered_func_data.nii.gz
		- doaroma: ICA AROMA (for denoising) using ICA_AROMA.py
			- output file is filtered_func_data_aroma.nii.gz
		- doconfound: creates a MAT file with all nuissance regressors (WM, CSF, motion)
			- output file is WM_200hpf.nii.gz, CSF_200hpf.nii.gz
		- dotrim: cuts off the last 5 seconds from the file and the confound regressors, also creates a design file with all regressors
			- output file is filtered_func_data_aroma_trim.nii.gz and resid_design_trim_200hpf_aroma.txt which includes: 
			- CSF regressor, WM regressors, 6 motion parameters, and scrubbing regressor (if needed)
		- regress out the residuals using fsl_glm and resid_design_trim_200hpf_aroma.txt
			- output is filtered_func_data_200hpf_aroma_residuals.nii.gz
		- converts BOLD data to standard space using FSL applwarp 

Input arguments include specific subjects (--subjects ['sub-01','sub-02']) or you can run all subjects using --all. Additionally, you can turn off any steps (see beginning of script) 

```
./preprocessing.py —all
```

3) `trim_param.py` is used to additionally cut more time points from the BOLD data and from the ratings data (set to 20 for this paper)

4) `trim_pre_regres.py` is used to trim an additional 20s from the beginning of the song after AROMA from all data, then redo the regressing out of confounds, standardizing, and then zipping for the ISC toolbox


### Analysis of continuous affective ratings scripts (in folder affect_ratings)

1) `cont_rating_preproc.m`: MATLAB function for loading and preprocessing (interpolating, downsampling to 10hz, and demeaning) raw ratings collected from fader in psychtoolbox; output is a text file with the first column being the time stamps and the second column being the rating value (between 0 10)

2) `cont_rating_avg.m`: MATLAB function for getting average ratings value (mean, median, t-statistic) across participants at each moment in time; output is a MATLAB variable called allfiles_emo and allfiles_enj with all preprocessed subjectlevel ratings together as well as average across all participants (ev_sadln_emo_mean_ratings_38ss_cut20.txt and ev_sadln_enj_mean_ratings_38ss_cut20.txt), which are used as regressors


### Behavioral/survey measures sript (in folder isc)

  `survey_score.py`: script for scoring the IRI from raw Qualtrics survey


### Voxelwise regression with affective ratings scripts (in voxelwise_regress)

  `fsl_glm_statsmodel.py`: script for setting up and running first-level, subjectwise stats models for voxelwise regression with ratings (uses generic .fsf files)  


### Whole-brain ISC differences between high and low empathy groups scripts (in folder isc)

1) `pairwise_cor_afni.py`: calculates voxelwise pairwise correlations of the entire BOLD signal using AFNI function 3dTcorrelate; outputs are a text file with the order of pairwise correlations; output is AFNI BRIC file

	```
	./pairwise_cor_afni.py --aroma --stand
	```

2) 3dLME_2grp.R: R script for using AFNI code to conduct linear mixed-effects (LME) model with a crossed random-effects formulation to identify voxels that had higher ISC values within the high Fantasy group as compared to the low Fantasy group as well as higher ISC values within rather than across groups.


### Intersubject phase synchronization (ipsp) scripts (in folder isps)

1) `isps_calc.m`: use the FUNSPY scripts (now in the ISCToolbox) to calculate phase synchronization (both ISPS and SBPS)

2) `isps_roi.m`: take the results from `isps_calc.m` and get the ISPS measure at each time point for each AAL-atlas defined ROI and between ROIS and then select the pre-defined ROIS we want from the 5 networks; output is a matlab file of ISPS at each TR for 52 ROIs (sadln_ipsts_52rois.mat)


### Predicting dynamic synchronization from continuous affective rating scripts (in folder isps)

1) `permute_regressor.m` : MATLAB function for getting a number of phase-scrambled regressors for calculating a null-distribution of regression coefficients; saves each one as a text file 

2) `isps_glm_perm.py` : Functions for running the GLM between continuous regressors (emotion ratings and acoustic features) and ISPS values for each ROI and collapsed into five networks of interest. The following functions are used
		- prep_reg: loads the enjoyment and sadness continuous ratings, convolves with HRF, and orthogonalizes one to the other 
		- load_brain_data: loads the Matlab variable with the ISPS data for each of 52 ROIs
		- run_glm: for 5000 permutations, calculates the regression coefficients between ratings (or phase randomized ratings)
		- network_glm: same as above but collapsed across 5 networks of interest (rather than 52 regions)
		- multi_correct: Benjamin-Hochberg FDR correction at p = 0.05 across regions/networks tested

3) 'sbps_glm_perm.py': Functions for running the GLM between continuous behavioral ratings and SBPS values for each ROI and collapsed into five networks of interest. In addition to all functions above, this script also includes:
		- cross-network: calculates the GLM between ratings and across-network based SBPS

### Continuous measure of musical features scripts (in mus_features folder)

  `mus_features.m`: matlab script for extracting brightness and RMS (using MIRtoolbox code) and normalizing it and cutting it for use as continuous regressor of neuroimaging data

## Authors

* **Matthew Sachs**


## Acknowledgments

* Jonas Kaplan for FSL-based fMRI preprocessing scripts and general analysis support
* Enrico Glerean for phase-scrambling scripts

%
% This is a demo file to load the parameters of your analysis. Please refer to the manual.
%


funpsy();

% EDIT BELOW HERE

%% INPUT DATA // An array of strings of valid files with location
% For data requirements look at the manual
song = 'happy'; 
basepath = '/Volumes/MusicProject/NaturalisticNetwork/funpsy/'; 
if ~exist(basepath)
    basepath = '/Volumes/MusicProject-1/NaturalisticNetwork/funpsy/'; 
end

cfg=[];  % this will contain the parameters of our expriment
scount = 0; 
for s=1:40
    
    if s == 7 || s == 18 || s == 28 || s == 30 
        continue
    end
    scount = scount + 1; 
    if length(num2str(s)) > 1
        subnum = num2str(s); 
    else
        subnum = sprintf('0%d',s); 
    end
    
    cfg.indata{scount}=[basepath 'data_' song '/sub-' subnum '_' song '_filtered_func_200hpf_standard_aroma.nii'];
end

% OUTPUT FOLDER // A folder where the software will write all data. 
% Disk space: you need at least twice the space occupied by the original datasets
outfolder = sprintf('%sfunpsy_out_%s_36ss_aroma_cut/',basepath,song)
if ~exist(outfolder)
    fprintf('Making %s outfolder\n', outfolder)
    mkdir(outfolder)
end
cfg.outpath=outfolder;   % Please, create the folder before using it

    

% SAMPLING RATE // In Hertz (1/TR)
cfg.Fs=1/1; 

% FILTER BAND // In Hertz, it has to be a four element vector
cfg.F=[0.025 0.04 0.07 0.09];

% FILTER DEVIATION // see 'help firpmord'
cfg.DEV=[0.05 0.01 0.05];

% BRAIN MASKS // masks that the software will use
cfg.coregistered_mask='./atlases/masks/MNI152_T1_2mm_brain_mask.nii';
cfg.compute_group_mask = 1;     % if = 1, it computes a group mask based on the power of each voxel
cfg.compute_spectrum = 0;       % if = 1, it computes a group frequency spectrum for each voxel

% NAME OF YOUR ANALYSIS SESSION
cfg.session_name = sprintf('funpsy_%s_36ss_aroma_cut',song);

% RUN PRE-ANALYSIS
sessionfile=funpsy_makepsess(cfg);  % validates input parameters
                                    % sessionfile will be the only variable needed in next function calls
                                    % it's a string with the path to the matlab file with all informations about
                                    % this phase analysis session


%% MAKE THE DATA // creates analytic signal
cfg=[];
cfg.sessionfile=sessionfile;
% cfg.compute_group_mask=1;     % overrides session settings
% cfg.compute_spectrum=1;       % overrides session settings
out = funpsy_makedata(cfg);     % filters, compute masks and creates the analytic signal

%% Whole brain analysis of intersubject synch.
cfg=[];
cfg.sessionfile=sessionfile;
out = funpsy_ips(cfg);          % compute whole brain intersubject phase synchrony
                                % results stored in out.results.ips

%% Pairwise ROIs analysis
cfg=[];
load atlases/aal_2mm_rois.mat
cfg.sessionfile=sessionfile;
cfg.rois=rois;
cfg.usemean = 1; 				% set usemean =1 if you want to just do a mean of the voxels in the region. Default is usemean =0, which uses first principal component
out = funpsy_makeroidata(cfg);  % extract roi time series based on the 1st principal component or the mean

% SBPS (Seed Based Phase Synchrony)
cfg=[];
cfg.sessionfile=sessionfile;
%cfg.useppc=1;					% Pairwise phase consitency is not implemented for SBPS since it needs testing
out = funpsy_sbps(cfg);         % takes a list of seeds/rois and computes full differential functional phase synchrony
                                % between each pair of seeds/rois.
                                % results stored in out.results.sbps


% ISBPS (Intersubject Seed Based Phase Synchrony)
cfg=[];
cfg.sessionfile=sessionfile;
%cfg.useppc=1;                  % Pairwise phase consitency is not implemented for ISBPS since it needs testing
out = funpsy_isbps(cfg);       	% takes a list of seeds/rois and computes full differential functional phase synchrony
                                % between each pair of seeds/rois.
                                % results stored in out.results.isbps


error('stop')
 
%%% Code for statistics is currently consuming too many resources. It needs a bit of work. See the readme file.
                        
%% Statistics
% SBPS (pairwise ROI analysis)
cfg=[];
cfg.sessionfile=sessionfile;
cfg.nonparam=1;     % recommended. If 0 uses parametric tests. 
cfg.parallel = 1;   % Experimental feature - uses parallel computing
cfg.perm=1000;     % for each ROI pair, does a non parametric test
                    % NOTE: the rois are already specified in data creation
                    % To modify them you need to rerun the analysis
cfg.statstype='sbps';
out = funpsy_stats(cfg);    % results in out.sbps_stats


% IPS (Whole)
cfg.statstype='ips';
out = funpsy_stats(cfg);    % results in out.ips_stats

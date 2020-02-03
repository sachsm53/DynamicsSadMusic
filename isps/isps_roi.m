clear all
close all

%% Analysis with FUNSPY Data
% Three measures IPS, SBPS, ISBPS.
song = 'happy';
basepath = '/Volumes/MusicProject/NaturalisticNetwork/funpsy'; 
if ~exist(basepath)
    basepath = '/Volumes/MusicProject-1/NaturalisticNetwork/funpsy'; 
end
analysis = sprintf('funpsy_out_%s_36ss_aroma_cut',song); 
analysispath = sprintf('%s/%s',basepath,analysis);
nTR = 495; 
TR = 1; 
numperms = 10000;

% load the psess file
resultmat = sprintf('%s/funpsy_%s_36ss_aroma_cut.mat',analysispath,song); 
load(resultmat)

% Load the variable variable ips
load([analysispath '/results/ips.mat']); % variable ips
%load([basepath psess.results.ips(2:end)]); 

% Define the networks of interest
aud = [79:84 87 88];
orbito = [5 6 9 10 25 26 27 28];
striatal = [71:78];
limbic = [29:34 37:42]; %insular, hippo, amygdala
dmn = [3 4 23 24 35 36 61 62 65:68 85 86 89 90]; %posterior cingulate, inferior parietal gyrus,precunues, angular (from FOX)

rr= horzcat(aud,striatal,limbic,orbito,dmn); %Get regions from all networks of interest - 50 ROIS
%rr = sort(rr);
%rr = [1:116]; %if you want to do all 116 ROIS
load(sprintf('%s/atlases/AAL/aal_labels.mat',basepath)); aal_labels; %116 ROIS
aal_labels{25} = 'Frontal_Med_Orb_L'; 
aal_labels{26} = 'Frontal_Med_Orb_R';

%Get labels and network name
for r=1:length(rr)
    fprintf('%d\t%s\n',r,aal_labels{rr(r)})
    %disp(r,aal_labels{rr(r)})
end


% load IPS result 4D matrix and compute a ROI-averaged IPS time series 
ipsts=[];
labels = cell(length(rr),1);
for r=1:length(rr)
    fprintf('%d %s\n',r,aal_labels{rr(r)}) 
    rID=rr(r);
    map=psess.rois(rID).map;
    temp=0;
    for i=1:size(map,1)
        temp=temp+squeeze(ips(map(i,1),map(i,2),map(i,3),:));
    end
    temp=temp/size(map,1);
    ipsts=[ipsts temp];
    %disp(aal_labels{rr(r)})
    labels{r,1} = aal_labels{rr(r)};
end

%cut first 20 seconds (DONT DO THIS FOR SADSH AND HAPPY)
ipsts_cut = ipsts(21:end,:);
numrois = size(rr,2);


%Time series from ROI 
% for r=1:length(rr)
%     rID=rr(r);
%     for s=1:psess.Nsubj
%        load([psess.roidata{s} '/' num2str(rID)]); %variable roits
%        allrois(:,s,r)=angle(roits);
%     end  
% end

% load IPS result 4D matrix and compute a ROI-averaged IPS time series 
ipsts=[];
labels = cell(length(rr),1);
for r=1:length(rr)
    fprintf('%d %s\n',r,aal_labels{rr(r)}) 
    rID=rr(r);
    map=psess.rois(rID).map;
    temp=0;
    for i=1:size(map,1)
        temp=temp+squeeze(ips(map(i,1),map(i,2),map(i,3),:));
    end
    temp=temp/size(map,1);
    ipsts=[ipsts temp];
    %disp(aal_labels{rr(r)})
    labels{r,1} = aal_labels{rr(r)};
end

%%save sbps output
ipsts_52.labels = labels; 
ipsts_52.data = ipsts;
ipstsout = sprintf('%s/%s_ipsts_52rois.mat',basepath,song);
if ~exist(ipstsout)
    fprintf('Saviing %s\n',ipstsout)
    save(ipstsout,'-struct','ipsts_52')
end


%Script for extracting features from MIRtoolbox

directory = pwd; 
stimdir = [pwd '/stimuli']; 
files = dir(stimdir);
filefolder = files(arrayfun(@(x) x.name(1), files) ~= '.');
filenames = {filefolder(:).name};
sort(filenames)
resolution = .1;
songlengths = [168,515,256];
         
%% Extracting and preprocessing features (RMS and brightness) for SADLN: 
for f = 1:length(filenames)
    clipname = filenames{f};
    clip = clipname(1:end-4); 
    songlength = songlengths(f); 
    CurrentAudioFile = sprintf('%s/%s',stimdir,clipname); 
    fprintf('%d\t%s\n',f,CurrentAudioFile);
    
    % Load the audio file
    a = miraudio(CurrentAudioFile);
    l = mirlength(a);
    L = mirgetdata(l);  

    % RMS
    rms = mirrms(a,'Frame',1,'s',1,'/1');
    RMS = mirgetdata(rms2);

    fprintf('%s\tRMS mean: %.4f\tRMS length: %d\n', clipname,mean(RMS),length(RMS))

    rms = mirrms(a,'Frame',1,'s',1,'/1');
    RMS = mirgetdata(rms);

    %Cut first 20s and demean
    rms_cut = RMS(21:end);
    rms_cut_dmean = detrend(rms_cut,'constant');

    % Brightness:
    spectrum = mirspectrum(a,'Frame',0.1,'s',1,'/1');
    br = mirbrightness(spectrum,'MinRMS',0.00001,'Cutoff',1500); 
    BR = mirgetdata(br);

    fprintf('%s\tBR mean: %.4f\tBR length: %d\n', clipname,mean(BR),length(BR))

    %Cut first 20s and demean
    br_cut = BR(21:end);
    br_cut_dmean = detrend(br_cut,'constant');

    %Save
    txtname = sprintf('%s/rms_ev_%s_cut.txt',directory,clipname);
    fileID = fopen(txtname,'w');
    fprintf(fileID,'%.4f\n',rms_win);
    fclose(fileID);

    txtname = sprintf('%s/br_ev_%s_cut.txt',directory,clipname);
    fileID = fopen(txtname,'w');
    fprintf(fileID,'%.4f\n',rms_win);
    fclose(fileID);
end



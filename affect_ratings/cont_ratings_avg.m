function outputdata = ratings_fmri_all(song,windowsize,windowspacing,meas,half,cut,lagname,detrend)

%Plot? Yes = 1, No = 0
%Save individual ratings? Yes = 1, No = 0

%Half can be 'half' - for first half; 'all' for all; 'test' for just second
%half

%Measure can be median or mean (the measure of the sliding window),
%EX: ratings_fmri_all('sadln',30,1,'mean','all','nocut','nolag','nodetend')
% ratings_fmri_all('sadln',30,1,'mean','all','cut','nolag')
% ratings_fmri_all('happy',30,1,'mean','all','nocut','nolag')

basefolder = '/Volumes/MusicProject/NaturalisticNetwork'; 
outfolder = sprintf('%s/fmri_ratings_summary/',basefolder);
ratingfolder = sprintf('%s/fmri_ratings_summary/ratings_fmri_preprocess',basefolder);

if strcmp(song, 'sadln')
    songlength = 515;
    if strcmp(detrend, 'nodetrend')
        files = dir([ratingfolder filesep '*sadln*raw.txt']); %not detrended
    else
        files = dir([ratingfolder filesep 'sub*none.txt']); %detrended
    end
    emoremovelist = {'sub-18','sub-30'};
    enjremovelist = {'sub-30','sub-31'};
elseif strcmp(song,'happy')
    songlength = 168;
    files = dir([ratingfolder filesep '*happy*raw.txt']); %detrended
    emoremovelist = {'sub-18'};
    enjremovelist = {'sub-28','sub-34'};
elseif strcmp(song,'sadsh')
    songlength = 256;
    files = dir([ratingfolder filesep '*sadsh*dmean.txt']); %detrended
    emoremovelist = {'sub-04','sub-30'};
    enjremovelist = {'sub-30','sub-15'};
end

%New removelist (with second half)
removelist = {'sub-07','sub-pi','sub-28','sub-30'}; %People to remove for bad brain data 

values = {'mean','sqmean','median','tstat'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Average files across people
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(cut,'cut')
    songlength = songlength-20; 
end


if strcmp(half,'all') 
    numfiles = length(files);  
elseif strcmp(half,'half')
    numfiles = 38;
elseif strcmp(half,'test')
    numfiles = 40;
end

%check first to see if exist 

jcount = 0; 
ecount = 0;
allfiles_enj = zeros(songlength,(numfiles/2) - length(removelist));
allfiles_emo = zeros(songlength,(numfiles/2) - length(removelist));
subjectlist_enj = cell(1,(numfiles/2) - length(removelist));
subjectlist_emo = cell(1,(numfiles/2) - length(removelist));
for i = 1:numfiles
    
    %Do second half (test) or first half (train) or all 
    if strcmp(half,'test') 
        ratefilename = files(i+38).name;
    else
        ratefilename = files(i).name;
    end
    
    condition = ratefilename(14:16); 
    subject = ratefilename(1:6);
    
    %Skip if subject in remove list
    if ~isempty(find(~cellfun(@isempty,strfind(removelist,subject))))
        fprintf('Skipping subject %s\n',subject);
        continue
    end

    subdata = load([ratingfolder filesep ratefilename]);
    
    if strcmp(cut, 'cut')
        if strcmp(lagname,'lag')
            subdata_lag = subdata(16:end-5); 
            subdata = subdata_lag;
        else
            subdata = subdata(21:end); %CUT IS DONE HERE
        end
    end
    
    %Put all ratings in one matrix
    if ~isempty(strfind(condition, 'joy'))
        if ~isempty(find(~cellfun(@isempty,strfind(enjremovelist,subject))))
            fprintf('Skipping subject %s for joy\n',subject);
            continue
        end
        jcount = jcount + 1;
        fprintf('Subject %s, joy %d\n',subject,jcount)
        allfiles_enj(:,jcount) = squeeze(subdata); 
        subjectlist_enj{jcount} = subject;
    elseif ~isempty(strfind(condition, 'emo'))
        if ~isempty(find(~cellfun(@isempty,strfind(emoremovelist,subject))))
            fprintf('Skipping subject %s for emo\n',subject);
            continue
        end
        ecount = ecount + 1; 
        fprintf('Subject %s, emo %d\n',subject,ecount)
        allfiles_emo(:,ecount) = squeeze(subdata);
        subjectlist_emo{ecount} = subject;
    end
end

%save variables
allfiles_emo = allfiles_emo(:,1:ecount);
allfiles_enj = allfiles_enj(:,1:jcount);

if strcmp(detrend, 'nodetrend')
    emo_var = sprintf('%s%s_emo_all_ratings_%s_%dss_%s_%s_%s',outfolder,song,detrend,ecount,cut,lagname,half);
    enj_var = sprintf('%s%s_enj_all_ratings_%s_%dss_%s_%s_%s',outfolder,song,detrend,jcount,cut,lagname,half);
else
    emo_var = sprintf('%s%s_emo_all_ratings_detrend_%dss_%s_%s_%s',outfolder,song,ecount,cut,lagname,half);
    enj_var = sprintf('%s%s_enj_all_ratings_detrend_%dss_%s_%s_%s',outfolder,song,jcount,cut,lagname,half); 
end

if ~exist(emo_var)
    save(emo_var,'allfiles_emo')
end
if ~exist(enj_var)
    save(enj_var, 'allfiles_enj')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Get average across people
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tstats_enj = [];
tstats_emo = [];
med_joy = [];
med_emo = [];
avg_enj = []; 
avg_emo = []; 

fprintf('\tGetting average, median, and tstat across people\n');
for v = 1:songlength
    [h,p,ci,stats_emo] = ttest(allfiles_emo(v,:));
    [h,p,ci,stats_enj] = ttest(allfiles_enj(v,:));
    tstats_enj = [tstats_enj stats_enj.tstat];
    tstats_emo = [tstats_emo stats_emo.tstat];
    
    m.e = median(allfiles_emo(v,:));
    m.j = median(allfiles_enj(v,:));
    med_joy = [med_joy m.j];
    med_emo = [med_emo m.e];
    
    avg_enj = [avg_enj mean(allfiles_enj(v,:))]; 
    avg_emo = [avg_emo mean(allfiles_emo(v,:))]; 
    
    %fprintf('Time point %d: Med Emo: %.4f, Med Enj: %.4f\n',v,m.e,m.j);
end

%% Get 2nd order as well 
%Demean all for plotting 
tstats_enj_norm = detrend(tstats_enj,'constant'); %REMOVES THE MEAN 
tstats_emo_norm = detrend(tstats_emo,'constant'); 
med_enj_norm = detrend(med_joy,'constant'); 
med_emo_norm = detrend(med_emo,'constant');
avg_enj_norm = detrend(avg_enj,'constant'); 
avg_emo_norm = detrend(avg_emo,'constant');
avg_emo_sq = avg_emo_norm.^2; 
avg_enj_sq = avg_enj_norm.^2;
avg_rating_emo = horzcat(avg_emo_norm',avg_emo_sq',med_emo_norm',tstats_emo_norm');
avg_rating_enj = horzcat(avg_enj_norm',avg_enj_sq',med_enj_norm',tstats_enj_norm');

%Save just the mean ratings
f.emo = sprintf('%s%s_emo_avg_measure_ratings_%dss_%s',outfolder,song,ecount,half);
f.joy = sprintf('%s%s_enj_avg_measure_ratings_%dss_%s',outfolder,song,jcount,half);
if ~exist(f.emo)
    save(f.emo,'avg_rating_emo')
    save(f.joy,'avg_rating_enj')
end

%Save just the mean ratings as text
emo_fid = sprintf('%s/ev_%s_emo_mean_ratings_%dss_%s.txt',basefolder,song,ecount,cut); %this was saved to funpsy and used as regressor
enj_fid = sprintf('%s/ev_%s_enj_mean_ratings_%dss_%s.txt',basefolder,song,jcount,cut);
if ~exist(emo_fid)
    fID = fopen(emo_fid,'w');
    fprintf(fID,'%.4f\n',avg_emo_norm);
    fclose(fID);
end
if ~exist(enj_fid)
    fID = fopen(enj_fid,'w');
    fprintf(fID,'%.4f\n',avg_enj_norm);
    fclose(fID);
end


%% Plot all measures (averaged across participants)
figure
subplot(2,1,1)
hold on
plot(avg_rating_emo(:,1)','g') %mean
plot(avg_rating_emo(:,2)','b') %squared
plot(avg_rating_emo(:,3)','r') %median
plot(avg_rating_emo(:,4)','black') %t test 
axis([0 songlength -6 6]);
xlabel 'TR', ylabel 'Normalized Values'
title('Raw Emotion Ratings and Measures')

subplot(2,1,2)
hold on
plot(avg_rating_enj(:,1)','g')
plot(avg_rating_enj(:,2)','b')
plot(avg_rating_enj(:,3)','r')
plot(avg_rating_enj(:,4)','black')
axis([0 songlength -6 6]);
xlabel 'TR', ylabel 'Normalized Values'
title('Raw Enjoyment Ratings and Measures')
hold off

for i = 1:4
    tenj = find(avg_rating_enj(:,i) == max(avg_rating_enj(:,i)));
    temo = find(avg_rating_emo(:,i) == max(avg_rating_emo(:,i))); 
    fprintf('Max emo %s located at %d, max enj %s located at %d\n',values{i},temo,values{i},tenj) 
en

function ratings_ev(savefile,plotit)

if ~isnumeric(plotit)
    error('Bad type!');
    return
end

%songlengths = [515,256,168];
%songname = {'snl_l','snl_s','hnl'};
songlengths = [515];
songname = {'snl_l'};
fmrisongname = {'sadln','sadsh','happy'};
lineColor = [0 0 .5];
textColor = [0 0 .5];
textSize = 12;

basefolder = '/Volumes/MusicProject/NaturalisticNetwork'; 
ratingfilepath = '/Volumes/MusicProject/NaturalisticNetwork/ratings_fmri'; 

removelist = {'sub-07','sub-pi'}; 

d= dir([ratingfilepath filesep 'sub*']);
d=d(~ismember({d.name},{'.','..'}));

%All data array 1 = SNL_L enjoy
%All data array 2 = SNL_L emo
%All data array 3 = SNL_S enjoy
%All data array 4 = SNL_S emo
%All data array 5 = HNL enjoy
%All data array 6 = HNL emo

resolution = .1;
fs = 10; %sample rate 100hz
fc = .0125; %cutoff frequency

for t = 1:length(songname)
    song = songname{t}; 
    fmrisong = fmrisongname{t};
    csnl = 0;
    songlength = songlengths(t); 
    qp = [0:resolution:songlength-1];
    for i = 1:(length(d)) %Run all 20 participants this way
        ratefilename = d(i).name;
        condition = ratefilename(end-10:end-8); 
        subject = ratefilename(1:6);
        %fprintf('%s\n',subject)

        if ~isempty(find(~cellfun(@isempty,strfind(removelist,subject))))
            fprintf('Skipping subject %s\n',subject);
            continue
        end
        

        if ~isempty(strfind(ratefilename, song))
            
            featfolder = sprintf('%s/fmri_analysis/%s/music/%s_model1_pre_200hpf.feat',basefolder,subject,fmrisong); 
            if exist(featfolder,'dir') ~= 7
                fprintf('%s does not exist. Moving on...\n',featfolder)
                continue
            end
            
            outfolder = sprintf('%s/fmri_analysis/%s/music/rating_ev',basefolder,subject);
            if exist(outfolder,'dir') ~= 7
                fprintf('%s\n \tMaking ev directory...\n', subject)
                mkdir(outfolder)
            end
            
            outname_norm = sprintf('%s/%s_%s_ev_dmean.txt',outfolder,fmrisong,condition); 
            if exist(outname_norm,'file') == 2
                fprintf('Found %s ev file. Moving on...\n',outname_norm)
                %delete(outname_norm)
                %continue
            end
         
            outname_preprocess_folder = sprintf('%s/fmri_ratings_summary/ratings_fmri_preprocess/%s_%s_%s_dmean.txt',basefolder,subject,fmrisong,condition); 
            if exist(outname_preprocess_folder,'file') == 2
                fprintf('Found %s ev file. Moving on...\n',outname_preprocess_folder)
                %delete(outname_norm)
                %continue
            end
            
            csnl = csnl + 1; 
            logfile = sprintf('%s/%s',ratingfilepath,ratefilename); 
            log = load(logfile);
            

            if isempty(log)
                fprintf('%d: subject %s song %s condition %s BLANK LOG FOUND\n',csnl,subject,song,condition);
                continue
            end

            %Interpolate data
            time = log(:,1);
            value = log(:,2);
            ratings = value/12.7;
            interpdata = interp1(time,ratings,qp,'spline')'; 
            fprintf('%s\t%s\t%s\t%d\t%d\n',subject,condition,song,length(ratings),length(interpdata))
            
            %Downsample data to 10hz
            ratings_ds = downsample(interpdata, 10);
            
            %make negatives zeros: 
            if min(ratings_ds) < 0
                fprintf('\tAt first min is %.6f,nax is %.2f\n',min(ratings_ds),max(ratings_ds))
                ratings_ds(ratings_ds < 0) = 0; 
                fprintf('\tNow, min is %.6f, max is %.2f\n',min(ratings_ds),max(ratings_ds))
            end   
            
            %Scale by detrending with mean
            ratings_norm = detrend(ratings_ds,'constant'); %2) scale by detrending with mean
            
            
            %if end of loop, plot
            if i == 40
                figure
                plot(ratings_ds)
                hold
                plot(ratings_norm,'g')
            end
            
            %Save file
            if savefile == 1
                if exist(outname_norm, 'file') ~= 2
                    fprintf('\tSaving %s %s ev file for %s. Number of ratings: %d\n',subject,condition,fmrisong,length(ratings_norm))
                    fileID = fopen(outname_norm,'w');
                    if fileID==-1
                        error('Cannot open file for writing: %s', outname_norm);
                    end
                    fprintf(fileID,'%.4f\n',ratings_norm);
                    fclose(fileID);
                else 
                    fprintf('\t%s ev file already exists.\n',outname_norm)
                end
                
                
               if exist(outname_preprocess_folder, 'file') ~= 2
                    fprintf('\tSaving %s %s ev file for %s. Number of ratings: %d\n',subject,condition,fmrisong,length(ratings_preprocess_sq))
                    fileID = fopen(outname_preprocess_folder,'w');
                    if fileID==-1
                        error('Cannot open file for writing: %s', outname_preprocess_folder);
                    end
                    fprintf(fileID,'%.4f\n',ratings_norm);
                    fclose(fileID);
                else 
                    fprintf('\t%s ev file already exists.\n',outname_preprocess_folder)
               end
            end
        end
    end
end

                
                
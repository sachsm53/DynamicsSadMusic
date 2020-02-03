function permute_regressor(song,reg,nPerm)
    
    clear all
    close all

    basepath = '/Volumes/MusicProject/NaturalisticNetwork/funpsy'; 
    if ~exist(basepath)
        basepath = '/Volumes/MusicProject-1/NaturalisticNetwork/funpsy'; 
    end


    %Output
    outfolder = sprintf('%s/glm_perm/fft_%s_%s',basepath,song,reg);
    if ~exist(outfolder)
        mkdir(outfolder)
    end


    %cut the happy and sadsh
    emo_cut = dlmread(sprintf('%s/regressors/%s_%s_ratings_emo_cut.txt',basepath,song,reg),'\n');
    enj_cut = dlmread(sprintf('%s/regressors/%s_%s_ratings_enj_cut.txt',basepath,song,reg),'\n');
    
    if strcmp(song,'happy') || strcmp(song,'sadsh')
        %Load raw ratings 
        emo = dlmread(sprintf('%s/ev_%s_emo_mean_ratings_36ss_nocut.txt',basepath,song),'\n');
        emo = emo';
        enj = dlmread(sprintf('%s/ev_%s_enj_mean_ratings_36ss_nocut.txt',basepath,song),'\n');
        enj = enj';
        fprintf('Cutting %s ratings\n',song)
        emo_cut = emo(21:end); 
        enj_cut = enj(21:end);
    else
        %Load raw ratings 
        emo = dlmread(sprintf('%s/ev_%s_emo_mean_ratings_36ss_cut.txt',basepath,song),'\n');
        emo = emo';
        enj = dlmread(sprintf('%s/ev_%s_enj_mean_ratings_36ss_cut.txt',basepath,song),'\n');
        enj = enj';
    end


    %% fake ratings via FFT
    clear fake_ratings_fft_emo
    clear fake_ratings_fft_enj
    rng(0);
    for i = 1:nPerm
        X=fft(zscore(emo_cut));
        AX=abs(X);
        PH=angle(X);
        randPH=PH(randperm(length(PH)));
        fake_ratings_fft_emo(:,i)=real(ifft(AX.*exp(j*randPH)));
        fname = sprintf('%s/emo_raw_ev%d.txt',outfolder,i); 
        if ~exist(fname)
            fprintf('Saving %s\n',fname)
            fid = fopen(fname,'w');
            fprintf(fid,'%.4f\n',fake_ratings_fft_emo(:,i));
            fclose(fid);
        end
    end

    for i = 1:nPerm
        X=fft(zscore(enj_cut));
        AX=abs(X);
        PH=angle(X);
        randPH=PH(randperm(length(PH)));
        fake_ratings_fft_enj(:,i)=real(ifft(AX.*exp(j*randPH)));
        fname = sprintf('%s/enj_raw_ev%d.txt',outfolder,i);
        if ~exist(fname)
            fprintf('Saving %s\n',fname)
            fid = fopen(fname,'w');
            fprintf(fid,'%.4f\n',fake_ratings_fft_enj(:,i));
            fclose(fid);
        end
    end
end
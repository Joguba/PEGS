close all
clear all

%for these datasets
for step=[1]

    %directory = sprintf( '/eno/jekollme/DATA/arneReanalyzedAgain/data%02d/', step); %where is the data right now
    %directory = '/Users/jekollmer/20160711/all2/';
    %directory =   sprintf( '/eno/jekollme/DATA/UniaxialCompression/20160701/Steps/step%02d/',step);
    %datapath =  sprintf( 'uniaxial/20160711/step%02d/', step); %where will the data be on the HPC
    %datapath =  sprintf( '/gpfs_common/share03/kdaniel/jekollme/uniaxial/20160701/step%02d/',step);
    %dataname =  sprintf( 'Step%02d*preprocessing.mat', step); %where will the data be on the HPC
    %dataname = '*Test.mat'; % for newtonizing
    %dataname =  '*preprocessing.mat'; %where will the data be on the HPC
    dataname = '*.jpg';
    %directory = sprintf( '/eno/cbkirber/DATA/20170208/1/');
    %dataname = '*preprocessing*.mat';
    directory = 'C:\Users\Mille\Desktop\Work\BD_Track\DATA\';
    datapath = sprintf( 'gpfs_common/share03/kdaniel/jsmille9/DATA/');
    
    
        
    files = dir([directory, dataname]);  %which files 
    nFrames = length(files); %how many files are we processing ?
    
    timeRequest = 1400; % how much compute time do we request ? 1440 minutes is a day and if you request more than a day you end up in a less favourable queue %ideally this should be proportional to nFrames

    nBatches = 1;
    perBatch = round(nFrames/nBatches);
    
    for batch = 1:nBatches
        jobname =  sprintf( 'jobfileBatchTest%02d-batch%02d.txt', step,batch); %where will the data be on the HPC

        fileID = fopen([directory,jobname],'w'); %open a jobfile (will erase old jobfile)

        fprintf(fileID, ['#! /bin/csh \n#BSUB -W ',num2str(timeRequest),' \n#BSUB -n 8 \n#BSUB -R span[ptile=8] \nsource /usr/local/apps/MATLAB/matlab2016a.csh \n']); %print jobfile header

        if batch < nBatches
            for frame = (batch-1)*perBatch+1:batch*perBatch %loop over these frames 
                fprintf(fileID, ['matlab -nojvm -nodisplay -singleCompThread -r "HPC_Conv(''',files(frame).name,''')"  > myout_',files(frame).name(1:end-18),'.txt \n']); %add a frame to process to jobfile
            end
        else
            for frame = (batch-1)*perBatch+1:nFrames %loop over these frames 
                fprintf(fileID, ['matlab -nojvm -nodisplay -singleCompThread -r "HPC_Conv(''',files(frame).name,''')"  > myout_',files(frame).name(1:end-18),'.txt \n']); %add a frame to process to jobfile
            end
        end

        fprintf(fileID, ['#BSUB -J "',files(frame).name(1:end-18),'" \n']);
        fprintf(fileID, ['#BSUB -o ',datapath,'out',num2str(batch),'.%%J \n']);
        fprintf(fileID, ['#BUSB -e ',datapath,'err',num2str(batch),'.%%J \n']);
        fclose(fileID);
    end
end
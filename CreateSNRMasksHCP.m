%parpool('local', 20) %% Name of cluster profile for batch job

clear all
disp(sprintf('Job Submitted: %s', datestr(now)));

%% NEED TO CHECK WHICH ON OF THESE SHOULD BE 1
SNR = 0;  %% Toggles whether to calculate SNR with signal measure only
tSNR = 1;  %% Toggles whether to calculate tSNR incorporating noise
tmasks = 0;  %% Toggles whether to apply tmask to each session for SNR maps
ConcatenateSessions = 1;  %% Toggles whether to concatenate subs for a group mask
MSCtemplate = 1;  %% Toggles whether to use MSC template or generic template
SessionTemplate = 0;  %% Toggles whether to use a separate template for each session

<<<<<<< HEAD
outdir = '/projects/p31161/';
datadir = '/projects/b1081/Lifespan/derivatives/';
subs = ['LS03'];
sessions = [1:5];
runs = [9,9,11,8,9];
addpath(genpath('/projects/b1081/Scripts'));

=======
outdir = '/projects/p31161/lateralizationVariants/';
%outdir = '/projects/b1081/member_directories/dperez/lateralizationVariants/';
subs = {'113922'};
session = {'REST1' 'REST2'};
acqType = {'LR' 'RL'};
addpath(genpath('/home/dcr8536/lateralizationVariants/'));
addpath(genpath('/projects/b1081/Scripts/'));
>>>>>>> fa1fa9d4ee8fb884fb36ba95e607d8aa4885a18f

disp(sprintf('Job Started: %s', datestr(now)));

if SessionTemplate == 1
    outdir= strcat(outdir, 'SNR_Maps/');
end

catData = [];

for i=1:numel(subs)
<<<<<<< HEAD
    for j=1:numel(sessions)
        for k=1:runs(j)    
=======
    for j=1:numel(session)
        for k=1:numel(acqType)    
            % Create template path for resting data
            HCPciftidir = ['/projects/p31161/lateralizationVariants/'];
            if tmasks == 1    
                % Load rest tmasks
                load(['/projects/p31161/lateralizationVariants/' subs{i} '_rfMRI_' session{j} '_' acqType{k} '_QC.mat'])
            end
    
            data = '/projects/p31161/lateralizationVariants/113922_rfMRI_REST1_LR/rfMRI_REST1_LR.4dfp.img';
            HCPciftidir = ['/home/dcr8536/lateralizationVariants/'];
>>>>>>> fa1fa9d4ee8fb884fb36ba95e607d8aa4885a18f
            if tmasks == 1    
                % Load rest tmasks
                load(sprintf('%s/preproc_FCProc/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_QC.mat',datadir,subs(i),sessions(j),subs(i),sessions(j)));
            end
<<<<<<< HEAD
            % sub-LS03_ses-1_task-rest_run-01_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii
            %data = sprintf('%s/postFCproc_CIFTI/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%02d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii',datadir,subs(i),sessions(j),runs(k));
            data = sprintf('%s/preproc_fmriprep-20.0.6/fmriprep/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%02d_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz',datadir,subs(i),sessions(j),subs(i),sessions(j),runs(k));
            [inputdata, f, v, etype] = read_4dfpimg_HCP(data);
=======
    
            %data = '/home/dcr8536/lateralizationVariants/Data/101006_rfMRI_REST2_LR/rfMRI_REST2_LR.4dfp.img';
            %data = [HCPciftidir 'Data/' subs{i} '_rfMRI_' session{j} '_' acqType{k} '/rfMRI_' session{j} '_' acqType{k} '.4dfp.img'];
            [inputdata] = read_4dfpimg_HCP(data);
>>>>>>> fa1fa9d4ee8fb884fb36ba95e607d8aa4885a18f
            [v, f, etype] = fcimage_attributes(data); %v= voxel size, f=number of frames or fourth dimension, etype = endian type (way data is stored)

            if tmasks == 1    
                resttmask = QC(k).tmask;
                disp(sprintf('tmask for rest file has %i good sample points, %s', sum(resttmask), datestr(now)));
                inputdata = inputdata(:,logical(resttmask));
            end
        
            if MSCtemplate == 1 && SessionTemplate == 1            
                disp('Calculating SNR')            
                if SNR == 1            
                    MeanSNR = mean(inputdata,2);    
                elseif tSNR == 1
                    signal = mean(inputdata,2);        
                    noise = std(inputdata,0,2);    
                    MeanSNR = signal./noise;                   
                end                
                if SNR == 1 && tmasks == 1        
                    outname = ([subs{i} '_' session{j} '_' acqType{k} '_SNRMap_tmasks.4dfp.img']);        
                elseif SNR == 1        
                    outname = ([subs{i} '_' session{j} '_' acqType{k} '_SNRMap.4dfp.img']);        
                elseif tSNR == 1 && tmasks == 1        
                    outname = ([subs{i} '_' session{j} '_' acqType{k} '_tSNRMap_tmasks.4dfp.img']);        
                elseif tSNR == 1        
                    outname = ([subs{i} '_' session{j} '_' acqType{k} '_tSNRMap.4dfp.img']);        
                end
            
                out_data = MeanSNR;
                %fout = [outdir outname];    
                fout= ['/Users/dianaperez/Desktop/' outname];
                disp('Writing .4dfp file')    
                write_4dfpimg(out_data,fout,etype);
                write_4dfpifh(fout,1,etype); %note that the 1 denotes this is only 1 volume large; etype should be the same as when the data was loaded    
                disp('Mapping volume to surface')                

                %map_vol_to_surface(fout, 'both', 'ribbon-constrained', 'MNI', subs(i),session(j), acqType(k));
                %map_vol_to_surface_HCP(fout,subs(i),session(j), acqType(k))            
                %5map_vol_to_surface(fout,'both','ribbon-constrained','MNI')            
                clear out_data    
                niftiout = strrep(fout, '.4dfp.img', '.nii');    
                disp('Creating NIFTI from .4dfp')    
                system(['nifti_4dfp -n ' fout ' ' niftiout]);                                
            else            
                disp(sprintf('Adding %i sample points to catData, %s', size(inputdata,2), datestr(now)));                
                catData = [catData inputdata];        
            end        
        end
%end        %% End subject loop for concatenated data

%      if SessionTemplate == 0
% %         
% %         disp('Calculating SNR')
% % 
% %         if SNR == 1
% %             
% %             MeanSNR = mean(catData,2);
% %     
% %         elseif tSNR == 1
% % 
% %             signal = mean(catData,2);
% %         
% %             noise = std(catData,0,2);
% %     
% %             MeanSNR = signal./noise;
% %                    
% %         end
% %     
% %         if SNR == 1 && tmasks == 1 && ConcatenateSubs == 1
% %         
% %             outname = ('AllSubs_SNRMap_tmasks_REST_AllSessions.4dfp.img');
% %     
% %         elseif SNR == 1 && ConcatenateSubs == 1
% %         
% %             outname = ('AllSubs_SNRMap_REST_AllSessions.4dfp.img');
% %         
% %         elseif tSNR == 1 && tmasks == 1 && ConcatenateSubs == 1
% %         
% %             outname = ('AllSubs_tSNRMap_tmasks_REST_AllSessions.4dfp.img');
% %         
% %         elseif tSNR == 1 && ConcatenateSubs == 1
% %         
% %             outname = ('AllSubs_tSNRMap_REST_AllSessions.4dfp.img');
% %         
% %         elseif SNR == 1 && tmasks == 1 && MSCtemplate == 1
% %             
% %             outname = ([subs{i} '_' 'SNRMap_tmasks_REST_MSCTemplate_AllSessions.4dfp.img']);
% %             
% %         elseif SNR == 1 && tmasks == 1
% %         
% %             outname = ([subs{i} '_' 'SNRMap_tmasks_REST_AllSessions.4dfp.img']);
% %         
% %         elseif SNR == 1 && MSCtemplate == 1
% %         
% %             outname = ([subs{i} '_' 'SNRMap_REST_MSCTemplate_AllSessions.4dfp.img']);
% %             
% %         elseif SNR == 1
% %         
% %             outname = ([subs{i} '_' 'SNRMap_REST_AllSessions.4dfp.img']);
% %             
% %         elseif tSNR == 1 && tmasks == 1 && MSCtemplate == 1
% %         
% %             outname = ([subs{i} '_' 'tSNRMap_tmasks_REST_MSCTemplate_AllSessions.4dfp.img']);
% %         
% %         elseif tSNR == 1 && tmasks == 1
% %         
% %             outname = ([subs{i} '_' 'tSNRMap_tmasks_REST_AllSessions.4dfp.img']);
% %             
% %         elseif tSNR == 1 && MSCtemplate == 1
% %         
% %             outname = ([subs{i} '_' 'tSNRMap_REST_MSCTemplate_AllSessions.4dfp.img']);    
% %         
% %         elseif tSNR == 1
% %         
% %             outname = ([subs{i} '_' 'tSNRMap_REST_AllSessions.4dfp.img']);
% %         
% %         end
% %     
%     
%         out_data = MeanSNR;
%         fout = [outdir outname];
%     
%         disp('Writing .4dfp file')
%     
%         write_4dfpimg(out_data,fout,etype);
%         write_4dfpifh(fout,1,etype); %note that the 1 denotes this is only 1 volume large; etype should be the same as when the data was loaded
%     
%         disp('Mapping volume to surface')
%     
%     
%         if MSCtemplate == 1     %% Use MSC-specific template
%         
%             sessionvox = strsplit(vcidlist(1).name, '_');
%         
%             map_vol_to_surface_MSCspecific(fout,subs{i},char(sessionvox(1)))
%         
%         
%         elseif MSCtemplate == 0     %% Use generic templated
%     
%             map_vol_to_surface(fout,'both','ribbon-constrained','711-2B'); %%% this function maps volume data to the the surface using a group estimate (it's not as precise as the individualized method we usually use for the MSC, but is good for a quick look at the data) - it can be found in scripts/WorkbenchScripts/map_vol_to_surface.m;
%             %both = both hemispheres or just one
%             %ribbon-constrained = how the interpolation is done (within the gray matter ribbon, or other options, including no interpolation/averaging
%             %711-2B = group template space that the data is in; WashU data is usually in 711-2B, but other datasets are often in MN
%         
%         end
%     
%         clear out_data
%     
%         niftiout = strrep(fout, '.4dfp.img', '.nii');
%     
%         disp('Creating NIFTI from .4dfp')
%     
%         system(['nifti_4dfp -n ' fout ' ' niftiout]);
%     
%     end

    catData = [];
    end
end
    
       %% End subject loop for individual data
    
    
    
    
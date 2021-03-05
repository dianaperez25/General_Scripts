clear all
disp(sprintf('Job Submitted: %s', datestr(now)));

%% NEED TO CHECK WHICH ON OF THESE SHOULD BE 1
SNR = 1;  %% Toggles whether to calculate SNR with signal measure only
tSNR = 0;  %% Toggles whether to calculate tSNR incorporating noise
tmasks = 1;  %% Toggles whether to apply tmask to each session for SNR maps
ConcatenateSessions = 1;  %% Toggles whether to concatenate subs for a group mask
MSCtemplate = 1;  %% Toggles whether to use MSC template or generic template
SessionTemplate = 1;  %% Toggles whether to use a separate template for each session

outdir = '/projects/p31161/';
datadir = '/projects/b1081/Lifespan/derivatives/';
subs = {'LS03'};
sessions = [1:5];
runs = [9,9,11,8,9];
addpath(genpath('/projects/b1081/Scripts'));
outname = 'LS03_tSNRmask.nii.gz';
outdir = '/projects/p31161/';

disp(sprintf('Job Started: %s', datestr(now)));

if SessionTemplate == 1
    outdir= strcat(outdir, 'SNR_Maps/');
end

inputdata = [];

for i=1:numel(subs)
    for j=1:numel(sessions)
        for k=1:runs(j)
            sub = subs{i};
            %% step 1: load the data I need -- tmask, functional data (unsure if I should be loading cifti or bold run) 
            
            % sub-LS03_ses-1_task-rest_run-01_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii
            %data = sprintf('%s/postFCproc_CIFTI/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%02d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii',datadir,subs(i),sessions(j),runs(k));
            data = sprintf('%s/preproc_fmriprep-20.0.6/fmriprep/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_run-%02d_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz',datadir,sub,sessions(j),sub,sessions(j),k);
            
            %% step 2: extract inputdata variable -- since I'm not using 4dfp.img files, I think I'm going to use the load_nii functions,
            % unclear if I will need endian type info at this point...
            % need to find voxel size and number of frames, can probably
            % easily find that info...
            %[inputdata, f, v, etype] = read_4dfpimg_HCP(data);
            %[v, f, etype] = fcimage_attributes(data); %v= voxel size, f=number of frames or fourth dimension, etype = endian type (way data is stored)
            [run_data] = load_untouch_nii_wrapper(data);
            
            inputdata = [inputdata run_data];
        end
            %% step 3: apply tmask to data -- will likely need to edit this to match 
            %dimensions bc I think inputdata is reshaped while in load_nii
            %func    
            if tmasks == 1 
                load(sprintf('%s/preproc_FCProc/sub-%s/ses-%d/func/sub-%s_ses-%d_task-rest_QC.mat',datadir,sub,sessions(j),sub,sessions(j)));
                resttmask = QCsub.tmask;
                disp(sprintf('tmask for rest file has %i good sample points, %s', sum(resttmask), datestr(now)));
                inputdata = inputdata(:,logical(resttmask));
            end
            
            %% step 4: calculate Mean SNR 
            if MSCtemplate == 1 && SessionTemplate == 1            
                disp('Calculating SNR')            
                if SNR == 1            
                    MeanSNR = mean(inputdata,2);    
                elseif tSNR == 1
                    signal = mean(inputdata,2);        
                    noise = std(inputdata,0,2);    
                    MeanSNR = signal./noise;                   
                end                
                
            %% step 5: make nifti and save it
            % use load_untouch_nii, not wrapper...
                outdat = load_untouch_nii(data);
                tempimg_out = int32(MeanSNR);
                d = [91,109,91,270];
                dims_bold = [d(1) d(2) d(3) size(tempimg_out,2)];
                outdat.img= reshape(tempimg_out,dims_bold);
                fout = [outdir outname];    
                disp('Writing nifti file') 
                [nifti] = make_nii(outdat);
                save_nii(nifti, fout);
                %write_4dfpimg(out_data,fout,etype);
                %write_4dfpifh(fout,1,etype); %note that the 1 denotes this is only 1 volume large; etype should be the same as when the data was loaded    
            %% step 6: map to surface -- This is going to be a lot of work probably....
            % function works with nifti yay!
            
                disp('Mapping volume to surface')                
                map_vol_to_surface(fout,'both','ribbon-constrained','MNI')            
                clear out_data    
                niftiout = strrep(fout, '.4dfp.img', '.nii');    
                disp('Creating NIFTI from .4dfp')    
                system(['nifti_4dfp -n ' fout ' ' niftiout]);                                
            else            
                disp(sprintf('Adding %i sample points to catData, %s', size(inputdata,2), datestr(now)));                
                catData = [catData inputdata];        
            end        

    end
end        %% End subject loop for concatenated data



    
    
    
    
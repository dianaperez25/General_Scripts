atlas = 'Seitzman300';
outDir_top = '/Users/dianaperez/Desktop/';
atlas_dir = '/Users/dianaperez/Box/Quest_Backup/Atlases/';
outDir_top = [FCdir '/corrmats_' atlas '/'];
num_smppts = 1000;
subject = 'LS03';
sessions = [1:5];
runs = [9,9,11,8,9];
data_folder = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/postFCproc_CIFTI/';
tmask_folder = '/Volumes/GRATTONLAB/Lifespan/BIDS/Nifti/derivatives/preproc_fmriprep-20.2.0/fmriprep/';

atlas_params = atlas_parameters_GrattonLab(atlas,atlas_dir);
roi_data = load_nii_wrapper(atlas_params.MNI_nii_file); %vox by 1

%% concatenate data and tmask

for i = 1:numel(sessions)
    for j = 1:runs(i)
        disp(sprintf('Loading data for session %d run %02d...', sessions(i), j))
        data_file = sprintf('%s/sub-%s/ses-%d/cifti_timeseries_normalwall/sub-%s_ses-%d_task-rest_run-%d_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii',data_folder, subject, sessions(i), subject,sessions(i),j);
       data = ft_read_cifti_mod(data_file);
        input_data = data.data;
        clear data;

        disp(sprintf('Loading tmask for session %d run %02d...', sessions(i), j))
        tmask = sprintf('%s/sub-%s/ses-%d/func/FD_outputs/sub-%s_ses-%d_task-rest_run-%d_desc-tmask_fFD.txt',tmask_folder,subject,sessions(i),subject,sessions(i),j);
        tmask_data = table2array(readtable(tmask)); 
        
        % apply tmask
        %masked_data = input_data(:,logical(tmask_data));
        %clear tmask_data;
        %clear input_data;
        
        % concatenate
        disp('concatenating data')
        catData = [catData input_data];        
        catTmask = [catTmask tmask_data'];
        disp(sprintf('data size is now %d by %d', size(catData,1), size(catData,2)))
        %clear masked_data;
        clear input_data;
        clear tmask_data;
    end
end

%% apply tmask
data = catData(:,logical(catTmask));
true_half = data(:,1:5454);
other_half = data(:,5455:7631);

corrmat = paircorr_mod();
fout_str = sprintf('%s/sub-%s_sess-%d_task-%s_corrmat_%s',outDir,subject,subInfo(i).session,subInfo(i).condition,atlas);
    
    figure_corrmat_GrattonLab(corrmat,atlas_params,-1,1);
    saveas(gcf,[fout_str '.tiff'],'tiff');
    close(gcf);
    
    % save out files
    save([fout_str '.mat'],'sess_roi_timeseries','sess_roi_timeseries_concat','tmask','tmask_concat','corrmat');
    





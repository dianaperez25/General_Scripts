%% script to apply brain mask to functional run

%% load whole brain mask, functional run, and QC file
% brain mask is logical that indicates which voxels are brain and which are
% non-brain; functional run is all signal values, need to load all runs; QC
% file has global signal calculated by fmriprep I think

wholebrain = load_untouch_nii('/projects/b1081/Lifespan/derivatives/preproc_fmriprep-20.0.6/fmriprep/sub-LS05/ses-1/func/sub-LS05_ses-1_task-rest_run-01_space-MNI152NLin6Asymres-2_desc-WholeBrain_Mask_FCProcess.nii.gz');
load('/projects/b1081/Lifespan/derivatives/preproc_FCProc/sub-LS05/ses-1/func/sub-LS05_ses-1_task-rest_QC.mat');

wb_mask = wholebrain.img;
globalsignalfp = [];

for j=1:numel(QCsub.runs)
    
    func_run = load_untouch_nii(sprintf('/projects/b1081/Lifespan/derivatives/preproc_FCProc/sub-LS05/ses-1/func/sub-LS05_ses-1_task-rest_run-%02d_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz', QCsub.runs(j)));
    func = func_run.img;
    func_mask = zeros(91,109,91,270);

    %% Apply whole brain mask to functional run
    for x = 1:size(wb_mask, 1)
        for y = 1:size(wb_mask, 2)
            for z = 1:size(wb_mask, 3)
                if wb_mask(x,y,z) == 1
                    func_mask(x,y,z,:) = func(x,y,z,:);
                end
            end
        end
    end

    %% Now average them

    globalSignalrun = mean(func_mask,1);
    globalSignalrun = mean(globalSignalrun,2);
    globalSignalrun = mean(globalSignalrun,3);
    globalSignalrun = reshape(globalSignalrun, 270, 1);
    globalsignalfp = [globalsignalfp; globalSignalrun];
    
end

R = corr(QCsub.global_signal, globalsignalfp);

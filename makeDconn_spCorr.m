function makeDconn_spCorr(subject, sessions, runs, writeDconn, make_spCorr)
%% Make Dconn and spCorr Map
% under construction
% need to finalize output name, how to deal with group average and making
% it part of output name for spcorr map. Looking into newStr =
% strrep(str,old,new)

%% PATHS
addpath(genpath('/projects/b1081/Scripts'))
data_folder = '/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/';
tmask_folder = '/projects/b1081/Lifespan/Nifti/derivatives/preproc_fmriprep-20.2.0/fmriprep/';
template_fname = '/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI/cifti_timeseries_normalwall/sub-INET003_ses-1_task-rest_run-01_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii';
groupAvgLoc = '/projects/b1081/Atlases/';
groupAvgName = '120_allsubs_corr';

%% VARIABLES
%subject = 'LS03';
%sessions = [2,5]; %[1,2,3,4], can be 1+
%runs = [9,9]; % assumes these are numbered 1:runnum, in order of sessions above
%output_file = sprintf('%s/dconn_cifti_normalwall/sub-%s_allsess_tmasked.dconn.nii', data_folder, subject); %CHANGE IF NEEDED
output_file = '/projects/p31161/sub-LS03_concat_timeseries.dtseries.nii';
writeDconn = 0;
make_spCorr = 1;

%% LOAD AND CONCATENATE DATA
catData = [];
catTmask = [];
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

%% MAKE DCONN
[dconn] = make_dconn(catData,catTmask,template_fname,writeDconn,output_file);

%% MAKE SPATIAL CORRELATION MAP
if make_spCorr
    make_spCorr_Map(1, dconn, output_file, groupAvgLoc, groupAvgName)
end

end 


function make_spCorr_Map(cortexOnly, dconn_dat, outputname, groupAvgLoc, groupAvgName)

if cortexOnly
    voxnum = 59412;
end

if ischar(dconn_dat)
    dconn_dat = ft_read_cifti_mod('/projects/b1081/Lifespan/Nifti/derivatives/postFCproc_CIFTI/dconn_cifti_normalwall/sub-LS03_allsess_tmasked.dconn.nii');
    cifti_corrmap = dconn_dat.data(1:voxnum,1:voxnum);
else
    cifti_corrmap = dconn_dat(1:voxnum, 1:voxnum);
end

clear dconn_dat

tic;
disp('making spatial correlations...')

cifti_corrmap(isnan(cifti_corrmap)) = 0;
cifti_corrmap1 = cifti_corrmap(1:voxnum,1:(voxnum/4));
disp(['First quarter of Cifti corrmap is ' num2str(size(cifti_corrmap1,1)) ' by ' num2str(size(cifti_corrmap1,2)) ': ' datestr(now)])
% Remove NaNs (produced if vertices have no data)


for dconnrow = 1:size(cifti_corrmap1,1)            
    cifti_corrmap1(dconnrow,:) = single(FisherTransform(cifti_corrmap1(dconnrow,:)));
end

disp(sprintf('First quarter of Fisher Transform Finished: %s', datestr(now)));
                    
disp(sprintf('Loading Template: %s', datestr(now)));

% Load group-average corrmat (assumes it's a dconn)
group = ft_read_cifti_mod([groupAvgLoc '/' groupAvgName '.dconn.nii']);
group = group.data(1:voxnum,1:voxnum);
group1 = group(1:voxnum,1:(voxnum/4));
sizedata = size(group1,2);

if cortexOnly==1
    % Load template variants file
    template = ft_read_cifti_mod('/projects/b1081/member_directories/bkraus/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
    template.data = [];
else
    template = ft_read_cifti_mod('/projects/b1081/MSC/TaskFC/FCProc_MSC05_motor_pass2/cifti_timeseries_normalwall_native_freesurf/vc39006_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii');
    template.data = [];
    template.hdr(6).dim = 1;
end

disp(['First quarter of Cifti group data is ' num2str(size(group1,1)) ' by ' num2str(size(group1,2)) ': ' datestr(now)])

disp(sprintf('Template Loaded: %s', datestr(now)));
            
            % Compare to group average
for i=1:sizedata
    template.data(i,1) = paircorr_mod(group1(:,i),cifti_corrmap1(:,i));
end
            
clear cifti_corrmap1
clear group1

disp(sprintf('First quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap2 = cifti_corrmap(1:voxnum,(voxnum/4)+1:(voxnum/2));

disp(['Second quarter of Cifti corrmap is ' num2str(size(cifti_corrmap2,1)) ' by ' num2str(size(cifti_corrmap2,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap2,1)
    cifti_corrmap2(dconnrow,:) = single(FisherTransform(cifti_corrmap2(dconnrow,:)));
end

disp(sprintf('Second quarter of Fisher Transform Finished: %s', datestr(now)));

group2 = group(1:voxnum,(voxnum/4)+1:(voxnum/2));
sizedata = size(group2,2);

disp(['Second quarter of Cifti group data is ' num2str(size(group2,1)) ' by ' num2str(size(group2,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data((i+(voxnum/4)),1) = paircorr_mod(group2(1:voxnum,i),cifti_corrmap2(1:voxnum,i));
end

clear cifti_corrmap2
clear group2

disp(sprintf('Second quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap3 = cifti_corrmap(:,(voxnum/2)+1:(voxnum*.75));

disp(['Third quarter of Cifti corrmap is ' num2str(size(cifti_corrmap3,1)) ' by ' num2str(size(cifti_corrmap3,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap3,1)
    cifti_corrmap3(dconnrow,:) = single(FisherTransform(cifti_corrmap3(dconnrow,:)));
end

disp(sprintf('Third quarter of Fisher Transform Finished: %s', datestr(now)));

group3 = group(:,(voxnum/2)+1:(voxnum*.75));
sizedata = size(group3,2);

disp(['Third quarter of Cifti group data is ' num2str(size(group3,1)) ' by ' num2str(size(group3,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data(i+(voxnum/2),1) = paircorr_mod(group3(:,i),cifti_corrmap3(:,i));
end

clear cifti_corrmap3
clear group3

disp(sprintf('Third quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap4 = cifti_corrmap(:,(voxnum*.75)+1:end);

disp(['Fourth quarter of Cifti corrmap is ' num2str(size(cifti_corrmap4,1)) ' by ' num2str(size(cifti_corrmap4,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap4,1)
    cifti_corrmap4(dconnrow,:) = single(FisherTransform(cifti_corrmap4(dconnrow,:)));
end

disp(sprintf('Fourth quarter of Fisher Transform Finished: %s', datestr(now)));

group4 = group(:,(voxnum*.75)+1:end);
sizedata = size(group4,2);

disp(['Fourth quarter of Cifti group data is ' num2str(size(group4,1)) ' by ' num2str(size(group4,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data(i+(voxnum*.75),1) = paircorr_mod(group4(:,i),cifti_corrmap4(:,i));
end

clear cifti_corrmap4
clear group4

disp(sprintf('Fourth quarter of Correlation Finished: %s', datestr(now)));


disp(['Final size of spatial correlation map is ' num2str(size(template.data,1)) ' by ' num2str(size(template.data,2)) ': ' datestr(now)])


% Write out the variants
if cortexOnly == 1
    out_fname = strrep(outputname, '_allsess_tmasked.dconn.nii', ['_vs_' groupAvgName '_cortex_corr']);
    ft_write_cifti_mod(out_fname,template)
    template.data = [];

elseif cortexOnly == 0
    out_fname = strrep(outputname, '_allsess_tmasked.dconn.nii', ['_vs_' groupAvgName '_subcort_corr']);
    ft_write_cifti_mod(out_fname,template)
    template.data = [];

end
 
disp('done')
toc
end
 

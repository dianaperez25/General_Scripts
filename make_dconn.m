function dconn_dat = make_dconn(input_data,tmask_data,save_dconn,varargin)
% 
%
% Script to take CIFTI input and output a dconn
% Toggle whether it is saved or just output
%
% EXAMPLE w filenames:
% make_dconn('/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI/cifti_timeseries_normalwall/sub-INET003_ses-1_task-rest_run-01_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii','/projects/b1081/iNetworks/Nifti/derivatives/preproc_fmriprep-20.0.6/fmriprep/sub-INET003/ses-1/func/FD_outputs/sub-INET003_ses-1_task-rest_run-01_desc-tmask_fFD.txt',1,'/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI/dconn_cifti_normalwall/sub-INET003_ses-1_task-rest_run-01_desc-tmasked.dconn.nii');
%
% INPUTS
% input_data: either fname string to input dtseries or already loaded
% matrix
% tmask_data: either fname string to input tmask or already loaded logical
% vector to use as tmask
% save_dconn: 1 or 0, save out dconn for viewing/use (if 1, varargin =
% string with output file name (full path))
%
% OUTPUTS
% dconn_dat: output dconn matrix for use
%
% Notes: requires scripts/files in our shared directory - make sure these
% are in your path
%
% CG, Based partially on BTK script from task/rest data - 06/2020
% 


% INPUT parameters
% dconn template: start with dtseries and modify accordingly (this assumes
% cortex + subcortex)
%template_fname = '/projects/b1081/MSC/TaskFC/FCProc_MSC01_mem_pass2/cifti_timeseries_normalwall_native_freesurf/vc38671_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii';
template_fname = '/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI/cifti_timeseries_normalwall/sub-INET003_ses-1_task-rest_run-01_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii';
voxnum = 65625; % Cortex + subcortex/cerebellum %Cortex only: voxnum = 59412;

%% FIND DATA
data_folder = '/projects/b1081/Lifespan/derivatives/postFCproc_CIFTI/cifti_timeseries_normalwall/';
df = readtable(input_data);
numdatas=size(df.sub,1);
    
%% LOAD DATA
disp('Loading data...')

for x = 1:numdatas
    runs = str2double(regexp(df.runs{i},',','split'));
    for y = 1:

    for x = 1:numel(data_files.name)
        data = ft_read_cifti_mod(data_files(x).name);
        input_data = data.data;
%% LOAD TMASK
% check if tmask data is a file name or an input matrix string
% (i) if filename, load tmask
% (ii) if input matrix string, proceed
disp('Loading tmask...')
if ischar(tmask_data)
    tmask_data = table2array(readtable(tmask_data));
end
tmask_data = logical(tmask_data); 

%% APPLY tmask to data
disp('Masking data with tmask ...')
data_masked = input_data(:,tmask_data); %assume input data is grayoordinate X time
clear input_data;

%% Take correlations among vertices
% a single run takes about 1 min
disp('Starting dconn correlations...')
tic;
dconn_dat = paircorr_mod(data_masked');
disp('Finished correlations')
toc

%% save out correlations (if desired)
if save_dconn
    
    % load template:
    template = ft_read_cifti_mod(template_fname);
    
    % modify as needed
    template.data = [];
    template.dimord = 'pos_pos';
    template.hdr.dim(7) = voxnum;
    template.hdr.dim(6) = template.hdr.dim(7);
    template.hdr.intent_code = 3001;
    template.hdr.intent_name = 'ConnDense';
    
    % substitute data
    template.data = dconn_dat;
    
    % write dconn
    % a single run takes about 30 seconds
    tic;
    disp('writing dconn...');
    ft_write_cifti_mod(varargin{1},template);
    disp('done');
    toc
end
    
 

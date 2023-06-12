% script to get the number of runs per session

data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/iNetworks/BIDS/Nifti/derivatives/preproc_fmriprep-20.2.0/fmriprep/';
subject = 'INET032';

subj_dir = [data_dir 'sub-' subject];
cd(subj_dir)
folder_contents = dir(pwd);
folders = {folder_contents.name}.';
folders = folders(contains(folders,'ses'));
sessions = numel(folders);

for ses = 1:sessions
    ses_dir = sprintf('%s/ses-%d/func/', subj_dir, ses);
    cd(ses_dir)
    bold_runs = dir(['*', 'rest', '*', 'preproc_bold.nii.gz']);
end

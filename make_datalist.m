% script to get the number of runs per session

%% EDIT THESE ACCORDING TO YOUR NEEDS
subject = 'INET060';
output_dir = '/Users/dianaperez/Desktop/datalists'; 
output_fname = sprintf('%s/sub-%s_datalist.xlsx', output_dir, subject);
%edit these if needed
data_dir = '/Volumes/fsmresfiles/PBS/Gratton_Lab/iNetworks/BIDS/Nifti/derivatives/preproc_fmriprep-20.2.0/fmriprep/';
task = 'rest';
TR = 1.1;
dropFr = 5;
topDir = '/scratch/dcr8536/iNetworks/Nifti/';
dataFolder = 'preproc_fmriprep-20.2.0';
confoundsFolder = 'preproc_fmriprep-20.2.0';
FDtype = 'fFD';

% creates the first row for the datalist, the title for each column
datalist = {'sub', 'sess', 'task', 'TR', 'dropFr', 'topDir', 'dataFolder', 'confoundsFolder', 'FDtype', 'runs'};

if ~exist(output_dir) % makes output directory if it doesn't exist
    mkdir(output_dir)
end

% below we will find how many sessions exist for this subject...
subj_dir = [data_dir 'sub-' subject];
cd(subj_dir)
folder_contents = dir(pwd); % ...by examining the contents of the subject's folders
folders = {folder_contents.name}.';
folders = folders(contains(folders,'ses'));
sessions = numel(folders); % here we get the number of folders that start with "ses-"


for ses = 1:sessions % for each session we will create a row...
    datalist{end+1,1} = subject; %... with subject ID ...
    datalist{end,2} = ses;% ... session number ...
    datalist{end,3} = task; % ... the task ...
    datalist{end,4} = TR; % ... TR ...
    datalist{end,5} = dropFr; % ... number of frames that will be dropped from the start ...
    datalist{end,6} = topDir; % ... top directory where necessary data is found ...
    datalist{end,7} = dataFolder; % ... directory where functional data is found ...
    datalist{end,8} = confoundsFolder; % ... directory where confound data is found ...
    datalist{end,9} = FDtype; % ... type of frame displacement measure we want ...
    % ... and last is number of runs, but we have to find how many runs there are
    % in each session first)
    % we'll do this by examining the contents of each session folder
    ses_dir = sprintf('%s/ses-%d/func/', subj_dir, ses);
    cd(ses_dir)
    bold_runs = dir(['*', 'rest', '*', 'preproc_bold.nii.gz']); % we find the files that have this string (the rest bold runs)
    runs = cellstr(extractBefore(extractAfter({bold_runs.name}.', 'run-'), '_space')); % and extract the numbers after the string "run-"
    %these lines below will create the string that will go in the last
    %cell, with the run numbers separated by commas
    new_cell = [];
    for run = 1:numel(runs); new_cell = sprintf('%s,%s', new_cell, runs{run}); end
    new_cell = extractAfter(new_cell, ',');
    datalist{end,10} = new_cell;
end

% now write the datalist to a spreadsheet
writecell(datalist, output_fname, "Filetype", "spreadsheet")
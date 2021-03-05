%% This script is testing using globus cli to transfer files from RDSS to 
% quest before running fmriprep preproc and transferring the outputs back
% to RDSS after the fmriprep is done

% subject ID for paths
sub = 'LS10';

% endpoint IDs
sourceID = 'c4e4d798-640e-11ea-960d-0afc9e7dd773';

% login to globus
%system('globus login')

% transfer command
[status, cmdout] = system(['globus transfer ' sourceID ':/Lifespan/BIDS/Nifti/sub-' sub ' 01f757e6-5ff0-11ea-960c-0afc9e7dd773:/Lifespan/sub-' sub ' --recursive --label "Testing CLI with fmriprep script"'])

% extract task ID

taskID = cmdout(103:end);

system(['globus task wait ' taskID]) 
% cd to location of batch script
cd /projects/b1081/member_directories/dperez

% submit job 
system('sbatch fmriprep_20.0.6')




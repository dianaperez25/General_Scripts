%% Script to get a variable with XYZ coordinates - subject specific
function [ordinates] = get_brainordinates(subject, type, dataLoc, output_dir)
cd /projects/p31161/Scripts/gifti-master/@gifti/private

%% VARIABLES/PATHS
workbenchdir = '/projects/b1081/Scripts/workbench2/bin_linux64/';
surface_R = [dataLoc subject '.R.' type '.32k_fs_LR.surf.gii'];
surface_L = [dataLoc subject '.L.' type '.32k_fs_LR.surf.gii'];
outfile_L = [output_dir subject '_' type '_metric_file_left.func.gii'];
outfile_R = [output_dir subject '_' type '_metric_file_right.func.gii'];

%create structure for gifti 
this = [];
this.data = [];
this.metadata = [];
this.label = [];

%% Step 1 - Make metric files from fsaverage surf.gii files
system([workbenchdir 'wb_command -surface-coordinates-to-metric ' surface_L ' ' outfile_L])
system([workbenchdir 'wb_command -surface-coordinates-to-metric ' surface_R ' ' outfile_R])

%% Step 2 - Load metric files 
left_metric = gifti_read(outfile_L, this);
right_metric = gifti_read(outfile_R, this);

%% Step 3 - Put coordinates together in one matrix
ordinates = [[left_metric.data{1,1}.data; right_metric.data{1,1}.data] [left_metric.data{1,2}.data; right_metric.data{1,2}.data] [left_metric.data{1,3}.data; right_metric.data{1,3}.data]]; 
end
 
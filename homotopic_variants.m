%% Script to quantify  how many vertices have variants in homotopic regions

%-------------------------------------------------------------
%% PATHS && VARIABLES
addpath(genpath('/projects/p31161/Scripts/'))
addpath(genpath('/projects/b1081/Scripts'))
workbenchdir = '/projects/b1081/Scripts/workbench2/bin_linux64/';
dataLoc = '/projects/p31161/lateralizationVariants/';
output_dir = '/projects/p31161/lateralizationVariants/';
subject = {'113619'}; %140117, 100206, 113619, 103010
%filename = [dataLoc subject '_ThresholdedVariantMap_SNRExclude_7.5.dtseries.nii'];
% ------------------------------------------------------------
for s = 1:numel(subject)
%% get brainordinates
[ordinates] = get_brainordinates(subject{s}, 'midthickness', dataLoc, output_dir);
filename = [dataLoc subject{s} '_ThresholdedVariantMap_SNRExclude_7.5.dtseries.nii'];

%% load cifti, if vertex is variant, find homotopic pair
cifti = ft_read_cifti_mod(filename);
vertices = cifti.data;
ht_pairs = [];
count = 1;
LH_vertices = cifti.data(cifti.brainstructure==1);

for vertex = 14:length(LH_vertices)
    %if vertices(vertex) > 0
        [left_cifti, left_metric, left_coordinates, right_cifti, right_metric, right_coordinates] = homotopic_pair(vertex, ordinates);
        ht_pairs(count,1) = left_cifti;
        ht_pairs(count,3) = left_metric;
        %ht_pairs(count,3) = [num2str(left_coordinates(1)) ', ' num2str(left_coordinates(2)) ', ' num2str(left_coordinates(3))];
        ht_pairs(count,2) = right_cifti;
        ht_pairs(count,4) = right_metric;
        %ht_pairs(count,6) = right_coordinates;
        count = count + 1;
    %end
end

filename2 = ['/projects/p31161/lateralizationVariants/ht_pairs_' subject{s} '.mat'];
save(filename2, 'ht_pairs')
end
%% HERE NEED TO WRITE SOMETHING TO QUANTIFY NUMBER OF VARIANTS 

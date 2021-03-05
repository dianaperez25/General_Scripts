%% Overlap Map Making Script for subgroups of left handers
% input: variant maps
% output: a map of proportion of subjects that have overlapping variants at given vertex

clear all
%--------------------------------------------------------------------------
%% PATHS
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'))
addpath(genpath('/Users/dianaperez/Box/Research/lateralizationVariants'))
dataLoc = '/Users/dianaperez/Box/Research/lateralizationVariants/variantMaps/'; %location of variant maps

%% VARIABLES
threshold = [5, 7.5, 10]; 
wholeGroup = 1; % 1 for 752 subs (includes sets of relatives) 0 for 384 subs (excludes relatives of subjects)
template = ft_read_cifti_mod([dataLoc '100206_ThresholdedVariantMap_SNRExclude_5.dtseries.nii']);

if wholeGroup == 1 % 
    load('goodSubs752.mat')
    subs = goodSubs752;
elseif wholeGroup == 0
    load('goodSubs384.mat')
    subs = goodSubs384;
end

LH = [];
MH = [];
RH = [];

for s = 1:length(subs)
    if subs(s,2) < -28
        LH = [LH; subs(s,:)];
    elseif subs(s,2) < 48
        MH = [MH; subs(s,:)];
    elseif subs(s,2) <= 100
        RH = [RH; subs(s,:)];
    end
end

[B, I] = sort(LH(:,2));
sortedLH = LH(I,:);

LHgroup1 = sortedLH(1:length(sortedLH)/2,:);
LHgroup2 = sortedLH((length(sortedLH)/2)+1:end,:);

for t = 1:length(threshold)
    disp(['Threshold ' num2str(threshold(t)) '%'])
    fileName = ['_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
       
        [LH1_files] = getfiles(dataLoc, fileName, LHgroup1);
        [LH2_files] = getfiles(dataLoc, fileName, LHgroup2);

        outfile1 = ['/Users/dianaperez/Desktop/HCP_LHgroup1_OverlapMap_752GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'];
        outfile2 = ['/Users/dianaperez/Desktop/HCP_LHgroup2_OverlapMap_752GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'];
        
        [overlap_map1] = makemap(LH1_files, template);
        [overlap_map2] = makemap(LH2_files, template);
        
        ft_write_cifti_mod(outfile1, overlap_map1);
        ft_write_cifti_mod(outfile2, overlap_map2);
        
end
    

function [files] = getfiles(dataLoc, fileName, subs)
%% Script to make list of paths to files
files = [];    
    for x = 1:length(subs)
        file = [dataLoc num2str(subs(x)) fileName];
        files{x,1} = file;
    end
end 

function [overlap_map] = makemap(files, template)
allsubs = [];
groupmap = [];
nfiles = length(files);
disp(['Calculating overlap for ' num2str(nfiles) ' subjects in this group...'])

for x = 1:nfiles
    % load subject data
    cifti = ft_read_cifti_mod(files{x});

    for q = 1:length(cifti.data)
        if cifti.data(q) > 0
            allsubs(q,x) = 1;
        elseif cifti.data(q) == 0
            allsubs(q,x) = 0;
        end
    end
    clear cifti
end

groupsum = sum(allsubs,2);
disp(['The maximum number of subjects that overlap in this group is ' num2str(max(groupsum))]) 
groupmap = sum(allsubs,2)/nfiles;
disp(['The maximum proportion of subjects that overlap in this group is ' num2str(max(groupmap))]) 
overlap_map = template;
overlap_map.data = groupmap;
end

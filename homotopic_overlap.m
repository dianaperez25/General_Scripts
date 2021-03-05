%% HOMOTOPIC Overlap Map Making Script
% input: variant maps
% output: a map of proportion of subjects that have overlapping variants at given vertex

clear all
%--------------------------------------------------------------------------
%% PATHS
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'))
addpath(genpath('/Users/dianaperez/Box/Research/lateralizationVariants'))
dataLoc = '/Users/dianaperez/Box/Research/lateralizationVariants/variantMaps/'; %location of variant maps

%% VARIABLES
threshold = [7.5]; 
wholeGroup = 0; % 1 for 752 subs (includes sets of relatives) 0 for 384 subs (excludes relatives of subjects)
matchGroups = 0; % match num of RH subs to num of LH subs
groups = {'LH', 'RH', 'MH'};

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

if matchGroups == 1
    RH = datasample(RH, length(LH), 'Replace', false);
end


for t = 1:length(threshold)
    disp(['Threshold ' num2str(threshold(t)) '%'])
    fileName = ['_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
    [LH_files, RH_files, MH_files] = makingdatafiles(dataLoc, fileName, LH, MH, RH, 2);
    allSubs_files = [LH_files; MH_files; RH_files];
    datafiles = {LH_files RH_files MH_files};
    template = ft_read_cifti_mod(LH_files{1});
    [overlap_map] = makemap('allSubs', allSubs_files, template);
    ft_write_cifti_mod(['/Users/dianaperez/Desktop/HCP_allSubs_OverlapMap_384GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'], overlap_map);
%     for g = 1:length(groups)
%         if matchGroups == 1 && wholeGroup == 1
%             outfile = ['/Users/dianaperez/Desktop/HCP_' groups{g} '_OverlapMap_752GoodSubs_matchedGroups_threshold_' num2str(threshold(t)) '.dtseries.nii'];
%         elseif matchGroups == 1
%             outfile = ['/Users/dianaperez/Desktop/HCP_' groups{g} '_OverlapMap_matchedGroups_threshold_' num2str(threshold(t)) '.dtseries.nii'];
%         elseif wholeGroup == 1
%             outfile = ['/Users/dianaperez/Desktop/HCP_' groups{g} '_OverlapMap_752GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'];
%         else
%             outfile = ['/Users/dianaperez/Desktop/HCP_' groups{g} '_OverlapMap_384GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'];
%         end
%         [overlap_map] = makemap(groups{g}, datafiles{g}, template);
%         ft_write_cifti_mod(outfile, overlap_map);
%     end
end

function [overlap_map] = makemap(group, files, template)
allsubs = [];
groupmap = [];
nfiles = length(files);
%disp(['Making overlap map for ' num2str(nfiles) ' subjects in ' group ' group...'])
hems = template.brainstructure(template.brainstructure>0);

for x = 1:nfiles
    % load subject data
    cifti = ft_read_cifti_mod(files{x});
    variants = [];
    ht_variants = [];
    count1 = 1;
    count2 = 1;
    left_hem = cifti.data(hems==1);
    right_hem = cifti.data(hems==2);
    for v = 1:64984
        if cifti.brainstructure(v) == -1
            variants(v,1) = 0;
        elseif cifti.brainstructure(v) == 1
            variants(v,1) = left_hem(count1);
            count1 = count1 + 1;
        elseif cifti.brainstructure(v) == 2
            variants(v,1) = right_hem(count2);
            count2 = count2 + 1;
        end
    end
    left_hem = variants(1:32492);
    right_hem = variants(32493:end);
    for i = 1:length(left_hem)
        if left_hem(i) > 0 && right_hem(i) > 0
            ht_variants(i,1) = 1;
        else
            ht_variants(i,1) = 0;
        end
    end
        
    for q = 1:length(ht_variants)
        if ht_variants(q) > 0
            allsubs(q,x) = 1;
        elseif ht_variants(q) == 0
            allsubs(q,x) = 0;
        end
    end
    clear cifti
end

groupsum = sum(allsubs,2);
disp(['The maximum number of subjects that overlap in the ' group ' group is ' num2str(max(groupsum))]) 
groupmap = sum(allsubs,2)/nfiles;
disp(['The maximum proportion of subjects that overlap in the ' group ' group is ' num2str(max(groupmap))]) 
overlap_map = template;
groupmap = [groupmap; groupmap];
groupmap = groupmap(template.brainstructure>0);
overlap_map.data = groupmap;
end

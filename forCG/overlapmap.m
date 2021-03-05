%% Overlap Map Making Script
% input: variant maps
% output: a map of proportion of subjects that have overlapping variants at given vertex


clear all
% addpath(genpath('/Users/dianaperez/Desktop'))
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'))
addpath(genpath('/Users/dianaperez/Box/Research/lateralizationVariants'))
dataLoc = '/Users/dianaperez/Box/Research/lateralizationVariants/variantMaps/'; %location of variant maps
threshold = [5,10]; 
wholeGroup = 1; % 1 for 752 subs (includes sets of relatives) 0 for 384 subs (excludes relatives of subjects)

if wholeGroup == 1
    nLH = 55;
    nMid = 43;
    nRH = 651;
    allSubs = xlsread('goodSubs_allInfo', 7); % sheet 7 for 752goodsubs
elseif wholeGroup == 0
    nLH = 27;
    nMid = 26;
    nRH = 330;
    allSubs = xlsread('goodSubs_allInfo', 3); % sheet 3 for 384goodsubs
end

LH = allSubs(1:nLH,1);
middle = allSubs(1:nMid,2);
RH = allSubs(1:nRH,3);

for t = 1:length(threshold)
    fileName = ['_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
    [LH_files, RH_files, middle_files] = makingdatafiles(dataLoc, fileName,LH,middle,RH,2);


    %create groupmap variable
    allsubs_RH = [];
    groupmap_RH = [];
    nfiles = size(RH_files,1);
    outfile_RH = ['/Users/dianaperez/Desktop/HCP_RightHand_OverlapMap_752GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'];
    outfile_LH = ['/Users/dianaperez/Desktop/HCP_LeftHand_OverlapMap_752GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'];
    outfile_middle = ['/Users/dianaperez/Desktop/HCP_middle_OverlapMap_752GoodSubs_threshold_' num2str(threshold(t)) '.dtseries.nii'];

    template = ft_read_cifti_mod(RH_files{1});

    for x = 1:nfiles
        % load subject data
        cifti_RH = ft_read_cifti_mod(RH_files{x});
        for q = 1:length(cifti_RH.data)
            if cifti_RH.data(q) > 0
                allsubs_RH(q,x) = 1;
            elseif cifti_RH.data(q) == 0
                allsubs_RH(q,x) = 0;
            end
        end
        clear cifti_RH
    end

    groupsumRH = sum(allsubs_RH, 2);
    groupmap_RH = sum(allsubs_RH, 2)/nRH;
    overlapmap_RH = template;
    overlapmap_RH.data = groupmap_RH;
    ft_write_cifti_mod(outfile_RH, overlapmap_RH);

    clear x
    clear q
    clear nfiles

    %% make group map for left handers

    %creat groupmap variable
    allsubs_LH = [];
    groupmap_LH = [];
    nfiles = size(LH_files,1);

    for x = 1:nfiles
        % load subject data
        cifti_LH = ft_read_cifti_mod(LH_files{x});

        for q = 1:length(cifti_LH.data)
            if cifti_LH.data(q) > 0
                allsubs_LH(q,x) = 1;
            elseif cifti_LH.data(q) == 0
                allsubs_LH(q,x) = 0;
            end
        end
        clear cifti_LH
    end

    groupsumLH = sum(allsubs_LH, 2);
    groupmap_LH = sum(allsubs_LH, 2)/nLH;
    overlapmap_LH = template;
    overlapmap_LH.data = groupmap_LH;
    ft_write_cifti_mod(outfile_LH, overlapmap_LH);


    %% make group map for middle subjects
    %creat groupmap variable
    allsubs_middle = [];
    groupmap_middle = [];
    nfiles = size(middle_files,1);

    for x = 1:nfiles
        % load subject data
        cifti_middle = ft_read_cifti_mod(middle_files{x});

        for q = 1:length(cifti_middle.data)
            if cifti_middle.data(q) > 0
                allsubs_middle(q,x) = 1;
            elseif cifti_middle.data(q) == 0
                allsubs_middle(q,x) = 0;
            end
        end
        clear cifti_middle
    end

    groupsum_mid = sum(allsubs_middle, 2);
    groupmap_middle = sum(allsubs_middle, 2)/nMid;
    overlapmap_middle = template;
    overlapmap_middle.data = groupmap_middle;
    ft_write_cifti_mod(outfile_middle, overlapmap_middle);
end
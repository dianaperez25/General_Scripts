%% PERMUTATION TESTS

clear all


%parpool('local', 1)
%--------------------------------------------------------------------------
%% PATHS
addpath(genpath('/projects/b1081/Scripts'));
%addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'));
%rootDir = '/Users/dianaperez/Box/Research/Lateralization_Variants/';
rootDir = '/projects/p31161/lateralizationVariants/';
addpath(genpath(rootDir))
varMapLoc = [rootDir 'variantMaps/']; %location of variant maps 
outputDir = [rootDir 'Permutation_Tests/'];
template = ft_read_cifti_mod([varMapLoc '109325_ThresholdedVariantMap_SNRExclude_7.5.dtseries.nii']);
%% VARIABLES
threshold = [7.5];1
wholeGroup = 0; % 1 for 752 subs (includes sets of relatives) 0 for 384 subs (excludes relatives of subjects)
matchGroups = 0; % match num of RH subs to num of LH subs
groups = {'allSubs'}; %'LH', 'RH', 'MH'
numperms = 200;% number of permutations
spCorrs=zeros(numperms,1);% initialize mat of values i'm permuting
LH = []; MH = []; RH = [];% initialize variables for each group - will contain subs files

if wholeGroup == 1 % 
    load([rootDir '/subjectData/goodSubs752.mat'])
    subs = goodSubs752;
elseif wholeGroup == 0
    load([rootDir '/subjectData/goodSubs384.mat'])
    subs = goodSubs384;
end



for s = 1:length(subs)
    if subs(s,2) < -28
        LH = [LH; subs(s,:)];
    elseif subs(s,2) < 48
        MH = [MH; subs(s,:)];
    elseif subs(s,2) <= 100
        RH = [RH; subs(s,:)];
    end
end

%% HERE LOAD MASK - need to make a mask of areas where we would expect to get variants
% and excludes areas where we never get variants


spCorrs=zeros(numperms,1);% initialize mat of values i'm permuting

if matchGroups == 1
    RH = datasample(RH, length(LH), 'Replace', false);
end
flip_switch = zeros(length(subs),1);
flip_switch(1:(length(flip_switch)/2)) = 1;
    
for t = 1:length(threshold)
    disp(['Threshold ' num2str(threshold(t)) '%'])
    fileName = ['_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
    
    %% HERE - make matrix same length as allSubs_files with 1's and 0's
    
%     true_left = [];
%     right2left = [];
    
    for p=1:numperms
        tic;
        rng('shuffle');
        ind = randperm(length(flip_switch))';
        
        for s = 1:length(flip_switch)
            rand_flip_switch(s,1) = flip_switch(ind(s));
        end

        left = [];
        flip_left = [];

        for x = 1:length(subs)
            % load subject data
            file = [varMapLoc num2str(subs(x)) fileName];
            cifti = ft_read_cifti_mod(file);
            new_cifti = [];
            count = 1;

            for i = 1:length(cifti.brainstructure)
                if cifti.brainstructure(i) > 0 
                    if cifti.data(count) >0
                        new_cifti(i,1) = 1;
                    else new_cifti(i,1) = 0;
                    end
                    count = count + 1;
                else
                    new_cifti(i,1) = cifti.brainstructure(i);
                end
            end

             left_hem = new_cifti(1:(length(new_cifti)/2));
             right_hem = new_cifti((length(new_cifti)/2)+1:end);
             indices = cifti.brainstructure(1:32492);
%             true_left(x,1) = left(cifti.brainstructure==1);
%             right2left(x,1) = right(cifti.brainstructure==1);
            
            if rand_flip_switch(x) == 0
                left(:,x) = left_hem(indices==1);
                flip_left(:,x) = right_hem(indices==1);
            elseif rand_flip_switch(x) == 1
                left(:,x) = right_hem(cifti.brainstructure==1);
                flip_left(:,x) = left_hem(cifti.brainstructure==1);
            end

        end

        left_sum = sum(left,2);
        %disp(['Permutation #' num2str(p) ': The maximum number of subjects that overlap in the left hemisphere is ' num2str(max(left_sum))]) 
        flip_left_sum = sum(flip_left,2);
        %disp(['Permutation #' num2str(p) ': The maximum number of subjects that overlap in the flipped left hemisphere is ' num2str(max(flip_left_sum))]) 
        %groupmap = sum(allsubs,2)/nfiles;
        overlap_left = sum(left,2)/length(subs);
        %disp(['Permutation #' num2str(p) ': The maximum proportion of subjects that overlap in the left hemisphere is ' num2str(max(overlap_left))]) 
        overlap_flip_left = sum(flip_left,2)/length(subs);
        %disp(['Permutation #' num2str(p) ': The maximum proportion of subjects that overlap in the flipped left hemisphere is ' num2str(max(overlap_flip_left))]) 
        
        spCorrs(p,t) = corr(overlap_left, overlap_flip_left);
        disp(['Permutation #' num2str(p) ': correlation ' num2str(spCorrs(p,t))]);
        toc
    end
end

save([outputDir '/VariantsvsFlippedVariants_spCorr_' num2str(numperms) 'permutations_3.mat'],'spCorrs')

% function [overlap_map] = makemap_flipped(group, files, template)
% allsubs = [];
% groupmap = [];
% nfiles = length(files);
% disp(['Making overlap map for ' num2str(nfiles) ' subjects in ' group ' group...'])
% 
% for x = 1:nfiles
%     % load subject data
%     cifti = ft_read_cifti_mod(files{x});
%     new_cifti = [];
%     count = 1;
%     
%     for i = 1:length(cifti.brainstructure)
%         if cifti.brainstructure(i) > 0 
%             new_cifti(i,1) = cifti.data(count);
%             count = count + 1;
%         else
%             new_cifti(i,1) = cifti.brainstructure(i);
%         end
%     end
%     
%     left_hem = new_cifti(1:(length(new_cifti)/2));
%     right_hem = new_cifti((length(new_cifti)/2)+1:end);
% 
%     new_cifti_all_vertices = [right_hem; left_hem]; %flipped cifti
%     new_cifti = new_cifti_all_vertices(cifti.brainstructure>0); %should I flip .brainstructure too? To account for more right hem vertices
%     
%     for q = 1:length(new_cifti)
%         if new_cifti(q) > 0
%             allsubs(q,x) = 1;
%         elseif new_cifti(q) == 0
%             allsubs(q,x) = 0;
%         end
%     end
%     
%     
% end
% 
% groupsum = sum(allsubs,2);
% disp(['The maximum number of subjects that overlap in the flipped ' group ' group is ' num2str(max(groupsum))]) 
% %groupmap = sum(allsubs,2)/nfiles;
% overlap_map = sum(allsubs,2)/nfiles;
% disp(['The maximum proportion of subjects that overlap in the flipped ' group ' group is ' num2str(max(groupmap))]) 
% % overlap_map = template;
% % overlap_map.data = groupmap;
% end
% 
% function [overlap_map] = makemap(group, files, template)
% allsubs = [];
% groupmap = [];
% nfiles = length(files);
% disp(['Making overlap map for ' num2str(nfiles) ' subjects in ' group ' group...'])
% 
% for x = 1:nfiles
%     % load subject data
%     cifti = ft_read_cifti_mod(files{x});
% 
%     for q = 1:length(cifti.data)
%         if cifti.data(q) > 0
%             allsubs(q,x) = 1;
%         elseif cifti.data(q) == 0
%             allsubs(q,x) = 0;
%         end
%     end
%     clear cifti
% end
% 
% groupsum = sum(allsubs,2);
% disp(['The maximum number of subjects that overlap in the ' group ' group is ' num2str(max(groupsum))]) 
% overlap_map = sum(allsubs,2)/nfiles;
% disp(['The maximum proportion of subjects that overlap in the ' group ' group is ' num2str(max(groupmap))]) 
% %overlap_map = template;
% %overlap_map.data = groupmap;
% end
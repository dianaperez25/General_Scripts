%% PERMUTATION TESTING - True Left hem vs Pseudo Left hem
% Pseudo left hem is the right hem flipped onto left hem
% input is variant maps 
% output is a matrix with r - values for n permutations, with the option of
% creating a plot, too.
% Written by Diana P -- January 2021

clear all

%--------------------------------------------------------------------------
%% PATHS
%--------------------------------------------------------------------------
%% for quest

%addpath(genpath('/projects/b1081/Scripts'));
%rootDir = '/projects/p31161/lateralizationVariants/';

%% for local computer

addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'));
rootDir = '/Users/dianaperez/Box/Research/Lateralization_Variants/';

%% for both

addpath(genpath(rootDir))
varMapLoc = [rootDir 'variantMaps/']; %location of variant maps 
outputDir = [rootDir 'Permutation_Tests/'];

%% HERE WILL NEED TO LOAD MASK - need to make a mask of areas where we would expect to get variants
% and excludes areas where we never get variants

%--------------------------------------------------------------------------
%% VARIABLES
%--------------------------------------------------------------------------

threshold = [7.5]; % specify thresholds for variants that want tested
wholeGroup = 0; % 1 for 752 subs (includes sets of relatives) 0 for 384 subs (excludes relatives of subjects)
matchGroups = 0; % 1 will match num of RH subs to num of LH subs
groups = {'MH'}; % specify groups that will undergo permutation testing
numperms = 1000; % number of permutations
spCorrs=zeros(numperms,1);% initialize mat of values i'm permuting
LH = []; MH = []; RH = []; % initialize mat of groups of subs; LH = left handers, RH = right handers, MH = 'middle' handers
allSubs = 0; % 1 will do permutation testing for all subjects across groups together, in addition to each group separately
plot_results = 0; % 1 will output a plot of the r values, in addition to the default matrix

%--------------------------------------------------------------------------
%% DIVIDING SUBJECTS INTO GROUPS
%--------------------------------------------------------------------------

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

if matchGroups == 1
    RH = datasample(RH, length(LH), 'Replace', false);
end

%--------------------------------------------------------------------------
%% BEGIN ANALYSIS
%--------------------------------------------------------------------------

for t = 1:length(threshold) % Will run this loop for every threshold specified under variables
        
    disp(['Threshold ' num2str(threshold(t)) '%'])
    fileName = ['_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii']; %how variant maps are named after the subject ID
    [LH_files, RH_files, MH_files] = makingdatafiles(varMapLoc, fileName, LH, MH, RH); %makes matrices inc file paths for all subjects in each group
    files = {MH_files}; % structure will all file paths
    
    
    for g = 1:(numel(groups)) % will go through this loop for each group specified under variables
        
        disp(['Group:' groups{g}])
        
        group_files = files{g};
        
        % matrix with 1's and 0's, will determine if sub's maps will be flipped
        % -- will be random
        flip_switch = zeros(length(group_files),1);
        flip_switch(1:(length(flip_switch)/2)) = 1;
        
    %--------------------------------------------------------------------------
    %% RUN PERMUTATIONS - MAKE RANDOMIZED MAPS
    %--------------------------------------------------------------------------
    
        for p = 1:numperms % will go through this loop for each permutation
            
            disp(['Permutation #' num2str(p)])
            rng('shuffle'); % need seed for randomizer
            
%             ind = randperm(length(files))';
%             rand_files = [];
%             for s = 1:length(files)
%                 rand_files{s,1} = files{ind(s)};
%             end
        
            % randomize 1's and 0's
            ind = randperm(length(flip_switch))';
            for s = 1:length(flip_switch)
                rand_flip_switch(s,1) = flip_switch(ind(s));
            end
            
                     
            [overlap_left, overlap_flip_left] = makemaps(groups{g}, group_files, rand_flip_switch);

            % get difference between left hem and flipped left hem
            %pseudo_diff = overlap_left - overlap_flip_left;
            
            % get correlation
            spCorrs(p,t) = corr(overlap_left, overlap_flip_left);
            disp(['Correlation between left hem and flipped left hem for group ' groups{g} ' is ' num2str(spCorrs(p,t))]) 

      
%         % Make new LH_files, MH files, and RH files variable with randomized subs
%         %datafiles = {rand_files(1:length(LH_files)) rand_files(length(LH_files)+1:length(RH_files)) rand_files(end-length(MH_files):end)};
%         template = ft_read_cifti_mod(LH_files{1});
%         pseudo_LH_files = rand_files(1:length(LH_files));
%         pseudo_RH_files = rand_files(end-length(RH_files)+1:end);
%         
%         [pseudo_RH_map] = makemap('pseudoRH', pseudo_RH_files, template);
%         [pseudo_LH_map] = makemap('pseudoLH', pseudo_LH_files, template);
%         
%         %perm_diff = RH_overlap_map - LH_overlap_map;
%         spCorrs(p,t) = corr(pseudo_RH_map, pseudo_LH_map);
%         disp(['Permutation #' num2str(p) ': correlation between left handers and right handers overlap map is ' num2str(spCorrs(p,t))]) 
        end 
        
        %save([outputDir '/RightvsLeftHanders_spCorr_' num2str(numperms) '_thresh_' num2str(threshold(t)) '_permutations.mat'],'spCorrs')
        save([outputDir '/VariantsvsFlippedVariants_' groups{g} '_spCorr_' num2str(numperms) '_thresh_' num2str(threshold(t)) '_permutations.mat'],'spCorrs')
        
        if plot_results == 1
            %% HERE, write code to output a plot of the results
            %% TO PLOT
            %plot(1:1000, pseudo_maps, 'bo', 'MarkerFaceColor', 'b')
            %hold on
            %plot(500, true_map, 'ro', 'MarkerFaceColor', 'r')
        end

    end
end


function [overlap_left overlap_flip_left] = makemaps(group, files, rand_flip_switch)
 % initialize variables to include left hem and flipped left hem
left = [];
flip_left = [];

allsubs = [];
groupmap = [];
nfiles = length(files);
disp(['Making overlap map for ' num2str(nfiles) ' subjects in ' group ' group...'])

for x = 1:nfiles
    % load subject data
    cifti = ft_read_cifti_mod(files{x});
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
    
    if rand_flip_switch(x) == 0
        left(:,x) = left_hem(cifti.brainstructure==1);
        flip_left(:,x) = right_hem(cifti.brainstructure==1);
    elseif rand_flip_switch(x) == 1
        left(:,x) = right_hem(cifti.brainstructure==1);
        flip_left(:,x) = left_hem(cifti.brainstructure==1);
    end
end


    %     new_cifti_all_vertices = [right_hem; left_hem]; %flipped cifti
%     new_cifti = new_cifti_all_vertices(cifti.brainstructure>0); %should I flip .brainstructure too? To account for more right hem vertices
%     
%     
%     for q = 1:length(new_cifti)
%         if new_cifti(q) > 0
%             allsubs(q,x) = 1;
%         elseif new_cifti(q) == 0
%             allsubs(q,x) = 0;
%         end
%     end
    
    
% will save these for when I add the mask
%left_masked = left(mask==1,:); 
%flip_left_masked = flip_left(mask==1,:);
left_masked = left;
flip_left_masked = flip_left;

% sum across subs to get overlap for each vertex
left_sum = sum(left_masked,2);
%disp(['Permutation #' num2str(p) ': The maximum number of subjects that overlap in the left hemisphere is ' num2str(max(left_sum))]) 
flip_left_sum = sum(flip_left_masked,2);
%disp(['Permutation #' num2str(p) ': The maximum number of subjects that overlap in the flipped left hemisphere is ' num2str(max(flip_left_sum))]) 

% divide by number of subs to obtain proportion
overlap_left = sum(left,2)/length(files);        
%disp(['Permutation #' num2str(p) ': The maximum proportion of subjects that overlap in the left hemisphere is ' num2str(max(overlap_left))])         
overlap_flip_left = sum(flip_left,2)/length(files);
%disp(['Permutation #' num2str(p) ': The maximum proportion of subjects that overlap in the flipped left hemisphere is ' num2str(max(overlap_flip_left))]) 
end

function [LH_files, RH_files, middle_files] = makingdatafiles(dataLoc, fileName, LH, middle, RH)
%% function to make list of paths to files

LH_files = [];
middle_files = [];
RH_files = [];
    
for x = 1:length(LH)
    file = [dataLoc num2str(LH(x)) fileName];
    LH_files{x,1} = file;
end
for x = 1:length(middle)
    file = [dataLoc num2str(middle(x)) fileName];
    middle_files{x,1} = file;
end
for x = 1:length(RH)
    file = [dataLoc num2str(RH(x)) fileName];
    RH_files{x,1} = file;
end
    
    
end


%% PERMUTATION TESTS

clear all
%--------------------------------------------------------------------------
%% PATHS
%addpath(genpath('/projects/b1081/Scripts'));
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'));
rootDir = '/Users/dianaperez/Box/Research/Lateralization_Variants/';
%rootDir = '/projects/p31161/lateralizationVariants/';
addpath(genpath(rootDir))
varMapLoc = [rootDir 'variantMaps/']; %location of variant maps 
outputDir = [rootDir 'Permutation_Tests/'];
template = ft_read_cifti_mod([varMapLoc '109325_ThresholdedVariantMap_SNRExclude_7.5.dtseries.nii']);
%% VARIABLES
threshold = [7.5]; 
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
%full_mask = ft_read_cifti_mod('/Users/dianaperez/Desktop/HCP_brainmask.dtseries.nii');
%mask = full_mask.data(1:29696);

spCorrs=zeros(numperms,1);% initialize mat of values i'm permuting

if matchGroups == 1
    RH = datasample(RH, length(LH), 'Replace', false);
end

for t = 1:length(threshold)
    disp(['Threshold ' num2str(threshold(t)) '%'])
    fileName = ['_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
    [LH_files, RH_files, MH_files] = makingdatafiles(varMapLoc, fileName, LH, MH, RH, 2);
    files = [LH_files; MH_files; RH_files];
    %template = ft_read_cifti_mod(LH_files{1});
    %% HERE - make matrix same length as allSubs_files with 1's and 0's
    flip_switch = zeros(length(files),1);
    flip_switch(1:(length(flip_switch)/2)) = 1;
%     true_left = [];
%     right2left = [];
    tic;
    for p=1:numperms
        rng('shuffle');
        ind = randperm(length(flip_switch))';
        
        for s = 1:length(flip_switch)
            rand_flip_switch(s,1) = flip_switch(ind(s));
        end

        left = [];
        flip_left = [];

        for x = 1:length(files)
            % load subject data
            cifti = ft_read_cifti_mod(files{x});
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
%             true_left(x,1) = left(cifti.brainstructure==1);
%             right2left(x,1) = right(cifti.brainstructure==1);
            
            if rand_flip_switch(x) == 0
                left(:,x) = left_hem(cifti.brainstructure==1);
                flip_left(:,x) = right_hem(cifti.brainstructure==1);
            elseif rand_flip_switch(x) == 1
                left(:,x) = right_hem(cifti.brainstructure==1);
                flip_left(:,x) = left_hem(cifti.brainstructure==1);
            end

        end
        
        %left_masked = left(mask==1,:);
        %flip_left_masked = flip_left(mask==1,:);
        left_masked = left;
        flip_left_masked = flip_left;
        left_sum = sum(left_masked,2);
        disp(['Permutation #' num2str(p) ': The maximum number of subjects that overlap in the left hemisphere is ' num2str(max(left_sum))]) 
        flip_left_sum = sum(flip_left_masked,2);
        disp(['Permutation #' num2str(p) ': The maximum number of subjects that overlap in the flipped left hemisphere is ' num2str(max(flip_left_sum))]) 
        %groupmap = sum(allsubs,2)/nfiles;
        overlap_left = sum(left,2)/length(files);
        disp(['Permutation #' num2str(p) ': The maximum proportion of subjects that overlap in the left hemisphere is ' num2str(max(overlap_left))]) 
        overlap_flip_left = sum(flip_left,2)/length(files);
        disp(['Permutation #' num2str(p) ': The maximum proportion of subjects that overlap in the flipped left hemisphere is ' num2str(max(overlap_flip_left))]) 
        
        spCorrs(p,t) = corr(overlap_left, overlap_flip_left);
        
    end
end

save([outputDir '/VariantsvsFlippedVariants_spCorr_' num2str(numperms) 'permutations_5.mat'],'spCorrs')
toc

%% TO PLOT
%plot(1:1000, pseudo_maps, 'bo', 'MarkerFaceColor', 'b')
%hold on
%plot(500, true_map, 'ro', 'MarkerFaceColor', 'r')

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
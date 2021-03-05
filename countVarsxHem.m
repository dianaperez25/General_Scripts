%% Script to quantify number of variants in each hemisphere
% I'm winging this.... I'm going to write this script to count, first the number of variants in each hemisphere 
% for each participant, and then I'll find a way to compare these numbers for left vs right handers
% I think I might also want to quantify the number of variant vertices in each hemisphere for each participant

clear all
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'))
addpath(genpath('/Users/dianaperez/Box/Research/lateralizationVariants'))
%--------------------------------------------------------------------------
%% PATHS 
dataLoc = '/Users/dianaperez/Box/Research/lateralizationVariants/variantMaps/'; %location of variant maps
%% VARIABLES
threshold = 7.5;
wholeGroup = 1;
groups = {'LH', 'RH', 'MH'};
LH = [];
MH = [];
RH = [];
varcount = {};
vertcount = {};

if wholeGroup == 1 % 
    load('goodSubs752.mat')
    subs = goodSubs752;
elseif wholeGroup == 0
    load('goodSubs384.mat')
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

subsxgroup = {LH RH MH};

for t = 1:length(threshold)
    disp(['Threshold ' num2str(threshold(t)) '%'])
    fileName = ['_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
    for g = 1:length(groups)
        var_count = [];
        vert_count = [];
        subjects = subsxgroup{g};
        [files] = getfilenames(dataLoc, fileName, subjects);
        for n = 1:length(files)
            var_count(n,1:2) = subjects(n,:);
            vert_count(n,1:2) = subjects(n,:);
            %% Step 1: Load variant maps with unique ID's
            map = ft_read_cifti_mod(files{n});

            %% Step 2: Take out non-brain vertices
            brain = map.brainstructure(map.brainstructure>0);

            %% Step 3: Index left hemisphere; Count number of variants in left hemisphere
            Lhem = map.data(brain==1);
            var_count(n,3) = sum(unique(Lhem)>0);
            vert_count(n,3) = sum(Lhem>0);
            
            
            %% Step 4: Index right hemisphere; Count number of variants in right hemisphere
            Rhem = map.data(brain==2);
            var_count(n,4) = sum(unique(Rhem)>0);
            vert_count(n,4) = sum(Rhem>0);
            
        end
        
        varcount{g} = var_count;
        vertcount{g} = vert_count;
    end
end

function [filenames] = getfilenames(dataLoc, fileName, subs)
filenames = [];
    for x = 1:length(subs)
        file = [dataLoc num2str(subs(x)) fileName];
        filenames{x,1} = file;
    end
end

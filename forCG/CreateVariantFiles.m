%% CreateVariantFiles.m
%This script makes the variant maps from spatial correlation maps

clear all

%% Paths
%change paths
workbenchdir = '/Applications/workbench/bin_macosx64/';
leftsurf = '/Users/dianaperez/Box/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.L.midthickness.32k_fs_LR.surf.gii';
rightsurf = '/Users/dianaperez/Box/Dependencies/32k_ConteAtlas_v2_distribute/Conte69.R.midthickness.32k_fs_LR.surf.gii';
dataLoc = '/Users/dianaperez/Box/Quest_Backup_Pre_1_20_2020/HCP_analyses/from_Ben/HCP_variants/spCorr/'; %location of spatial correlation maps
fileName = '_allRuns_mean_vs_120_avg_corr_LR_corr.dtseries.nii'; %string after subject ID in file names
SNRpath = '/Users/dianaperez/Box/HCP_variants/bottomBrainMask.dtseries.nii'; %location of brain mask
outfilepath = '/Users/dianaperez/Box/Research/lateralizationVariants/variantMaps/'; %location of output

threshold = [5, 10];  %% Thresholds used to calculate variants (lowest % or correlation values)
SNRexclusion = 1;  %% Toggles whether to exclude variants based on SNR, 1 = exclude, 0 = don't exclude
ExcludeBySize = 1; % 1 -> will exclude variants < 50 vertices
wholeGroup = 1; % 1 for 752 subs (includes sets of relatives) 0 for 384 subs (excludes related subs)

if wholeGroup == 1
    nLH = 55;
    nMid = 43;
    nRH = 651;
    allSubs = xlsread('goodSubs_allInfo', 7); % sheet 3 for 384goodsubs, sheet 7 for 752goodsubs
elseif wholeGroup == 0
    nLH = 27;
    nMid = 26;
    nRH = 330;
    allSubs = xlsread('goodSubs_allInfo', 3); % sheet 3 for 384goodsubs, sheet 7 for 752goodsubs
end

%create variables with subs ID's
LH = allSubs(1:nLH,1);
middle = allSubs(1:nMid,2);
RH = allSubs(1:nRH,3);

%create variables with paths to subs' variant maps
[LH_files, RH_files, middle_files] = makingdatafiles(dataLoc, fileName, LH, middle, RH, 1)

for t = 1:length(threshold)
    %% Makes variant maps for left handers
    for x = 1:length(LH_files)
        subject = LH(x);
        %reads cifti files from the txt files
        cifti_LH = ft_read_cifti_mod(LH_files{x});
        % apply brain mask
        if SNRexclusion == 1
            SNRmap = ft_read_cifti_mod(SNRpath);
            SNRmap.data = SNRmap.data(1:59412,:);
            SNRexclude = find(SNRmap.data == 1);
            cifti_LH.data(SNRexclude,1) = NaN;
        end
        
        %findsvertices that are below the threshold, making those vertices = 1 
        %and all other vertices = 0
        cifti_LH_threshold = find(cifti_LH.data < prctile(cifti_LH.data,threshold(t)));
        cifti_LH_thresh_dat = zeros(size(cifti_LH.data));
        cifti_LH_thresh_dat(cifti_LH_threshold,1) = 1;
        cifti_LH_final_dat = zeros(size(cifti_LH.data));

        for w = 1:length(cifti_LH.data)
            if cifti_LH_thresh_dat(w) == 1
                cifti_LH_final_dat(w) = 1;
            end
        end
        
        % create names for output files
        outfileLH = [outfilepath '/Left/' num2str(subject) '_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii']; %final variant map
        outfilewbLH = [outfilepath '/temp/' num2str(subject) '_LH_Variant_Size_50_TempVariantMap_' num2str(threshold(t)) '.dtseries.nii']; %temp variant map

        %writes the file in cifti format
        cifti_LH.data = cifti_LH_final_dat;
        ft_write_cifti_mod(outfileLH, cifti_LH)
        
        % numbers that variants, clustering adjacent vertices that = 1
        % together and assigning each cluster a number
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfileLH ' 0 0 0 0 COLUMN ' outfilewbLH ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

        % writes file with numbered clusters in preparation for size
        % exclusion
        cifti_LH = ft_read_cifti_mod(outfilewbLH);

        if ExcludeBySize == 1             
            % exclusion criteria is set to 15 vertices (any variant less
            % than 15 vertices big will be excluded
            [cifti_LH.data] = ExcludeVariantSize(cifti_LH.data, subject, threshold(t), 50);
        end 
        
        % writes size-excluded variant masks
        ft_write_cifti_mod(outfileLH, cifti_LH)
    end

    
    % makes variant maps for right handed participants
    for x = 1:length(RH_files)
        subject = RH(x);
        
        %reads cifti files from the txt files
        cifti_RH = ft_read_cifti_mod(RH_files{x});
        
        % apply brain mask
        if SNRexclusion == 1
            SNRmap = ft_read_cifti_mod(SNRpath);
            SNRmap.data = SNRmap.data(1:59412,:);
            SNRexclude = find(SNRmap.data == 1);
            cifti_RH.data(SNRexclude,1) = NaN;
        end
        
        %finds vertices that are below the threshold, making those vertices = 1 
        %and all other vertices = 0
        cifti_RH_threshold = find(cifti_RH.data < prctile(cifti_RH.data,threshold(t)));
        cifti_RH_thresh_dat = zeros(size(cifti_RH.data));
        cifti_RH_thresh_dat(cifti_RH_threshold,1) = 1;
        cifti_RH_final_dat = zeros(size(cifti_RH.data));

        for w = 1:length(cifti_RH.data)
            if cifti_RH_thresh_dat(w) == 1
                cifti_RH_final_dat(w) = 1;
            end
        end

        outfileRH = [outfilepath '/Right/' num2str(subject) '_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
        outfilewbRH = [outfilepath '/temp/' num2str(subject) '_RH_Variant_Size_50_TempVariantMap_' num2str(threshold(t)) '.dtseries.nii'];

        %writes the file in cifti format
        cifti_RH.data = cifti_RH_final_dat;
        ft_write_cifti_mod(outfileRH, cifti_RH)
        
        % numbers that variants, clustering adjacent vertices that =1
        % together and assigning each cluster a number
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfileRH ' 0 0 0 0 COLUMN ' outfilewbRH ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

        % writes file with numbered clusters in preparation for size
        % exclusion
        cifti_RH = ft_read_cifti_mod(outfilewbRH);

        if ExcludeBySize == 1             
            % exclusion criteria is set to 15 vertices (any variant less
            % than 15 vertices big will be excluded
            [cifti_RH.data] = ExcludeVariantSize(cifti_RH.data, subject, threshold(t), 50);
        end 
        
        % writes size-excluded variant masks
        ft_write_cifti_mod(outfileRH, cifti_RH)
    end

    
    %% Make variant maps for "middle" handed participants
    for x = 1:length(middle_files)
        subject = middle(x);
        
        %reads cifti files from the txt files
        cifti_mid = ft_read_cifti_mod(middle_files{x});
        
        % apply brain mask
        if SNRexclusion == 1
            SNRmap = ft_read_cifti_mod(SNRpath);
            SNRmap.data = SNRmap.data(1:59412,:);
            SNRexclude = find(SNRmap.data == 1);
            cifti_mid.data(SNRexclude,1) = NaN;
        end
        
        %finds vertices that are below the threshold, making those vertices = 1 
        %and all other vertices = 0
        cifti_mid_threshold = find(cifti_mid.data < prctile(cifti_mid.data,threshold(t)));
        cifti_mid_thresh_dat = zeros(size(cifti_mid.data));
        cifti_mid_thresh_dat(cifti_mid_threshold,1) = 1;
        cifti_mid_final_dat = zeros(size(cifti_mid.data));

        for w = 1:length(cifti_mid.data)
            if cifti_mid_thresh_dat(w) == 1
                cifti_mid_final_dat(w) = 1;
            end
        end

        outfilemid = [outfilepath '/Middle/' num2str(subject) '_ThresholdedVariantMap_SNRExclude_' num2str(threshold(t)) '.dtseries.nii'];
        outfilewbmid = [outfilepath '/temp/' num2str(subject) '_midH_Variant_Size_50_TempVariantMap_' num2str(threshold(t)) '.dtseries.nii'];

        %writes the file in cifti format
        cifti_mid.data = cifti_mid_final_dat;
        ft_write_cifti_mod(outfilemid, cifti_mid)
        
        % numbers that variants, clustering adjacent vertices that =1
        % together and assigning each cluster a number
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfilemid ' 0 0 0 0 COLUMN ' outfilewbmid ' -left-surface ' leftsurf ' -right-surface ' rightsurf])

        % writes file with numbered clusters in preparation for size
        % exclusion
        cifti_mid = ft_read_cifti_mod(outfilewbmid);

        if ExcludeBySize == 1             
            % exclusion criteria is set to 15 vertices (any variant less
            % than 15 vertices big will be excluded
            [cifti_mid.data] = ExcludeVariantSize(cifti_mid.data, subject, threshold(t), 50);
        end 
        
        % writes size-excluded variant masks
        ft_write_cifti_mod(outfilemid, cifti_mid)
    end
end
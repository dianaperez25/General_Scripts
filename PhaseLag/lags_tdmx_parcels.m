function lags_tdmx_parcels(task)
%% Script originally from Ryan for doing lag analyses
% previously found in:
% /data/nil-bluarc/raichle/ryan/scripts/lags/cifti_lags_toolbox/
% CG edits: 5.25.2018

%% Setup
addpath('/data/nil-bluearc/ances/RELEASE/MATLAB');
addpath('/data/nil-bluearc/ances/Buckner/Lags_multi');
addpath('/data/nil-bluearc/raichle/ryan/scripts/lags/cifti_lags_toolbox');
wb_dir = '/usr/local/pkg/workbench1.2.3/bin_rh_linux64';
clear subjects patids
num_nodes = 333; %7320;  % cortex [Ryan usually does on 4k surface]

outdir = ['/data/nil-bluearc/GMT/Caterina/lags/' task '/']; % set directory for saving out images here
if ~exist(outdir)
    mkdir(outdir);
end 

%%% Set parameters
lags = -3:3;    % range of TR shifts -- (max(lags)-1)*TR must be > lag limit
lag_lim = 4;    % lag limit (in seconds)

%% MSC subjects
root_dir = ['/data/nil-bluearc/GMT/Caterina/TaskFC/FC_Parcels/' task '/'];
%subjects = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
subjects = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'};
tr = 2.2;
fd_thresh = .2; % IMPORTANT -- must match motion criteria used uring preproc; may want to change to 0.5 eventually
patids = [1:10]; % all datasets have 10 sessions [but not all sessions are good]

%% some atlas information
atlas_params = atlas_parameters('Parcels','/data/cn5/caterina/TaskConn_Methods/all_data/');
% parcel giftis - needed for creating parcel files
watershed_L = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_L_watershedmerge_0.35_tweaked.func.gii'];
watershed_R = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'];


%% Body

min_frames = max(lags)+1; % min. contiguous frames (in seconds)

% initialize group matrices
grp_lags = single(nan(num_nodes,num_nodes));    % peak lags
grp_ZL = single(nan(num_nodes,num_nodes));      % zero-lag correlation

grp_lags_nans = single(zeros(num_nodes,num_nodes));
grp_ZL_nans = single(zeros(num_nodes,num_nodes));

for s = 1:numel(subjects)
    
    subj = subjects{s};
    disp(['Processing ' subj]);
    %patids = textread([root_dir '/TESTING/' subj '/' subj '_patids.lst'],'%s');
    
    % initialize subject matrices
    subj_lags = single(nan(num_nodes,num_nodes,length(patids))); % peak lags
    subj_ZL = single(nan(num_nodes,num_nodes,length(patids)));   % zero-lag correlation
    
    %subj_lags_nans = single(zeros(num_nodes,num_nodes));
    %subj_ZL_nans = single(zeros(num_nodes,num_nodes));
    

    % load data for this subject
    subj_parcel_data = load([root_dir subj '_parcel_timecourse.mat']); % loads parcel_time and tmask_all
    
    for p = 1:length(patids)
        tic
        patid = patids(p);
        
        % load BOLD time series
        %cifti = ft_read_cifti_mod(['/data/nil-bluearc/raichle/ryan/lags/MSC/surface_lags/' subj '/cifti_timeseries_normalwall_atlas_4k/' patid '_LR_surf_subcort_333_4k_fsLR_smooth3.4_modesubcort.dtseries.nii']);
        % do whatever you want below to include/exclude certain brain regions
        %mask = cifti.brainstructure(cifti.brainstructure>0); % to make sure you're getting good signal voxels
        %mask1 = mask<3;
        %mask2 = mask<3; % cortex only
        %BOLD1 = nanmean(cifti.data(mask1,:))';
        %BOLD2 = cifti.data(mask2,:)';
        BOLD1 = subj_parcel_data.parcel_time{p}; % assuming this is a mean over vertices?
        BOLD2 = subj_parcel_data.parcel_time{p};
        
        % defined vertices (can include subject-wise SNR mask here)
        %good1 = true(size(Avg1,2),1);
        %good2 = true(size(Avg2,2),1);
        good1 = true(size(BOLD1,2),1); % assume all parcels are good
        good2 = true(size(BOLD2,2),1); % assume all parcels are good       
        
        % load motion file; can alternatively load your own temporal mask; always make sure frames masked during preproc are masked here
        %fd_file = [root_dir '/' subj '/' patid '/movement/' patid '_faln_dbnd_xr3d.FD'];      
        %fd = dlmread(fd_file); fd = fd(:,1);
        %format = true(size(fd));
        %format(fd > fd_thresh) = false;
        %format(1:2) = false; % pre-steady state
        format = logical(subj_parcel_data.tmask_all{p}); % tmask
                     
        %% Do the lagged correlation/covariance computation of TD matrices
        % CG: This has been simplified now that windowing is gone
        % So move temporal mask for each shift of the time series
        % no minimum for now... we'll see later (Anish has used minimum of 60 seconds, Ryan has used 20 seconds)
        % analysis done per session currently, then averaged over sessions
        % Only keep blocks >= min_frames
        
        % create a tmask that accounts for consecutive frames
        % [We shouldn't need this]
        consec_good_frames = 0;
        format_use = format;
        for i = 1:length(format)
            if format(i)
                if consec_good_frames > 0
                    consec_good_frames = consec_good_frames + 1;
                else
                    consec_good_frames = 1;
                    start = i;
                end
            else
                if consec_good_frames > 0
                    consec_good_frames = 0;
                    stop = i-1;
                    if numel(start:stop) < min_frames
                        format_use(start:stop) = false;
                    end               
                end
            end
        end
        
        % Compute cross-covariance function
        % produces a nodes x nodes x # of lags matrix
        Cov = lagged_cov_v3(BOLD1(:,good1),BOLD2(:,good2),format_use,max(lags));
        
        % PLOTS - for testing (CG)
        %for i = 1:7
        %    figure_corrmat_network_generic(Cov(:,:,i),atlas_params,-2,2)
        %    title(['Lag: ' num2str(lags(i))]);
        %    save_fig(gcf,[outdir 'CrossCov_' subj '_' num2str(patids(p)) '_lag' num2str(i) '.pdf']);
        %end
        %close('all');
        
        [pl,pc] = parabolic_interp_v2_TEMP(Cov,tr); % peak lag = 333 x 333 matrix (time delay)

        pl(abs(pl) > lag_lim) = nan; % Long lags excluded [means they didn't estimate well]
        pc(abs(pl) > lag_lim) = nan;
        
        % running sum for averaging across sessions
        %temp_lags = single(zeros(size(subj_lags)));
        %temp_ZL = single(zeros(size(subj_lags)));
        %temp_pZL = single(zeros(size(subj_lags)));
        %temp_lags(good1,good2) = pl; %vec_mat_flip(pl,'TD'); % CG - already in matrix format
        %temp_ZL(good1,good2) = Cov(:,:,lags==0); %vec_mat_flip(Cov(:,lags==0),'ZL'); % CG - already in matrix format
        %temp_pZL(good1,good2) = pc; % CG added this output to save peak covariance
        
        % CG: changed this to save out all of the sessions, since not huge
        %subj_lags = cat(3,subj_lags,temp_lags);
        %subj_lags = nansum(subj_lags,3);
        %subj_ZL = cat(3,subj_ZL,temp_ZL);
        %subj_ZL = nansum(subj_ZL,3);
        subj_lags(good1,good2,p) = pl;
        subj_ZL(good1,good2,p) = Cov(:,:,lags==0);
        subj_pZL(good1,good2,p) = pc;        
        
        % running sum of nans
        %subj_lags_nans = subj_lags_nans + isnan(temp_lags);
        %subj_ZL_nans = subj_ZL_nans + isnan(temp_ZL); 

        % save out upper triangulars - CG not going to bother to do this for now
        %temp_lags = vec_mat_flip(temp_lags);
        %temp_ZL = vec_mat_flip(temp_ZL);
        %save([outdir '/' subj '_' num2str(p) '_TD.mat'],'temp_lags');
        %save([outdir '/' subj '_' num2str(p) '_ZL.mat'],'temp_ZL');
        
        toc
    end
    
    % take average
    %subj_lags_mean = subj_lags ./ (numel(patids) - subj_lags_nans);
    %subj_ZL_mean = subj_ZL ./ (numel(patids) - subj_ZL_nans);
    subj_lags_mean = nanmean(subj_lags,3);
    subj_ZL_mean = nanmean(subj_ZL,3);
    subj_pZL_mean = nanmean(subj_pZL,3);
    
    subj_lags_D = nanmean(subj_lags,3)./nanstd(subj_lags,[],3);
    
    % covariance matrix to correlation matrix
    d = zeros(size(subj_ZL_mean));
    d(logical(eye(length(subj_ZL_mean)))) = sqrt(diag(subj_ZL_mean));
    subj_ZL_mean = FisherTransform(d^(-1)*subj_ZL_mean/d); % CG - added Fisher transform

    d = zeros(size(subj_pZL_mean));
    d(logical(eye(length(subj_pZL_mean)))) = sqrt(diag(subj_pZL_mean));
    subj_pZL_mean = FisherTransform(d^(-1)*subj_pZL_mean/d); % CG - added Fisher transform
    
    % weighted lag projection for subject 
    % LP: CG - lag projection (average time delay relative to the rest of the
    % brain)..  low dimensional, loses a lot of information since different
    % lag threads are exclusionary
    % BUT lets you get a quick look if things look normal
    % Ryan will change this section of the code
    % now lag projections are weighted toward things that have high
    % correlations [at zero lag]
    lag_weights = tan((pi/2)*(1-abs(subj_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
    lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
    lag_weights(isnan(subj_lags_mean)) = nan;
    subj_lags_mean_wghtd = subj_lags_mean.*lag_weights;
    subj_LP = nansum(subj_lags_mean_wghtd)./nansum(lag_weights); 

    % save out these calculations
    save([outdir subj '_lagmats_persession.mat'],'subj_lags','subj_ZL','subj_pZL');
    save([outdir subj '_lagmats_mean.mat'],'subj_lags_mean','subj_ZL_mean','subj_pZL_mean','subj_lags_D','subj_LP');

    % save lag projection as a CIFTI
    assign_data_to_parcel_cifti(subj_LP',watershed_L,watershed_R,outdir,[subj '_lagprojection']);

    % make some figures
    [bb i] = sortrows([atlas_params.mods_array; subj_LP]',[1 2]);
    atlas_params_new = atlas_params;
    atlas_params_new.sorti = i;
    figure_corrmat_network_generic(subj_lags_D,atlas_params_new,-5,5);
    save_fig(gcf,[outdir subj '_lagD_matrix.pdf']);
    figure_corrmat_network_generic(subj_lags,atlas_params_new,-2,2);
    save_fig(gcf,[outdir subj '_lag_matrix.pdf']);
    close('all');
    
    % running sum for averaging across subjects
    %grp_lags = cat(3,grp_lags,subj_lags_mean);
    %grp_lags = nansum(grp_lags,3);
    %grp_ZL = cat(3,grp_ZL,subj_ZL_mean);
    %grp_ZL = nansum(grp_ZL,3);
    grp_lags(:,:,s) = subj_lags_mean;
    grp_ZL(:,:,s) = subj_ZL_mean;
    grp_pZL(:,:,s) = subj_pZL_mean;
    
    % running sum of nans
    %grp_lags_nans = grp_lags_nans + isnan(subj_lags_mean);
    %grp_ZL_nans = grp_ZL_nans + isnan(subj_ZL_mean);

    clear pl Cov subj_* temp*
    
end

% take average
%grp_lags_mean = grp_lags ./ (numel(subjects) - grp_lags_nans);
%grp_ZL_mean = grp_ZL ./ (numel(subjects) - grp_ZL_nans);
grp_lags_mean = nanmean(grp_lags,3);
grp_ZL_mean = nanmean(grp_ZL,3);
grp_pZL_mean = nanmean(grp_pZL,3);
grp_lags_D = nanmean(grp_lags,3)./nanstd(grp_lags,[],3);

% weighted lag projection for group
lag_weights = tan((pi/2)*(1-abs(grp_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
lag_weights(isnan(grp_lags_mean)) = nan;
grp_lags_mean_wghtd = grp_lags_mean.*lag_weights;
grp_LP = nansum(grp_lags_mean_wghtd)./nansum(lag_weights);

% save this out
save([outdir 'Group_lagmats_mean.mat'],'grp_lags_mean','grp_ZL_mean','grp_pZL_mean','grp_LP','grp_lags_D');

% save lag projection as a CIFTI
assign_data_to_parcel_cifti(grp_LP',watershed_L,watershed_R,outdir,['Group_lagprojection']);


% make some figures
[bb i] = sortrows([atlas_params.mods_array; grp_LP]',[1 2]);
atlas_params_new = atlas_params;
atlas_params_new.sorti = i;
figure_corrmat_network_generic(grp_lags_D,atlas_params_new,-5,5);
save_fig(gcf,[outdir 'Group_lagD_matrix.pdf']);
figure_corrmat_network_generic(grp_lags,atlas_params_new,-2,2);
save_fig(gcf,[outdir 'Group_lag_matrix.pdf']);
figure_corrmat_network_generic(grp_ZL_mean,atlas_params_new,-1,1);
save_fig(gcf,[outdir 'Group_ZeroLagFC_matrix.pdf']);
figure_corrmat_network_generic(grp_pZL_mean,atlas_params_new,-1,1);
save_fig(gcf,[outdir 'Group_PeakLagFC_matrix.pdf']);

end
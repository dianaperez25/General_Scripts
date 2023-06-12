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

outdir = '/data/nil-bluearc/GMT/Caterina/lags/'; % set directory for saving out images here

%%% Set parameters
lags = -3:3;    % range of TR shifts -- (max(lags)-1)*TR must be > lag limit
lag_lim = 4;    % lag limit (in seconds)

%% MSC subjects
root_dir = '/data/nil-bluearc/raichle/ryan/MSC';
%subjects = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'};
subjects = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tr = 2.2;
fd_thresh = .2; % IMPORTANT -- must match motion criteria used uring preproc; may want to change to 0.5 eventually
patids = 1:10;

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
    patids = textread([root_dir '/TESTING/' subj '/' subj '_patids.lst'],'%s');
    
    % initialize subject matrices
    subj_lags = single(nan(num_nodes,num_nodes)); % peak lags
    subj_ZL = single(nan(num_nodes,num_nodes));   % zero-lag correlation
    
    subj_lags_nans = single(zeros(num_nodes,num_nodes));
    subj_ZL_nans = single(zeros(num_nodes,num_nodes));
    

    for p = 1:numel(patids)
        tic
        patid = patids{p};
        
        % load BOLD time series
        cifti = ft_read_cifti_mod(['/data/nil-bluearc/raichle/ryan/lags/MSC/surface_lags/' subj '/cifti_timeseries_normalwall_atlas_4k/' patid '_LR_surf_subcort_333_4k_fsLR_smooth3.4_modesubcort.dtseries.nii']);
        % do whatever you want below to include/exclude certain brain regions
        mask = cifti.brainstructure(cifti.brainstructure>0); % to make sure you're getting good signal voxels
        mask1 = mask<3;
        mask2 = mask<3; % cortex only
        BOLD1 = nanmean(cifti.data(mask1,:))';
        BOLD2 = cifti.data(mask2,:)';
        
        % defined vertices (can include subject-wise SNR mask here)
        good1 = true(size(Avg1,2),1);
        good2 = true(size(Avg2,2),1);
        
        % load motion file; can alternatively load your own temporal mask; always make sure frames masked during preproc are masked here
        fd_file = [root_dir '/' subj '/' patid '/movement/' patid '_faln_dbnd_xr3d.FD'];      
        fd = dlmread(fd_file); fd = fd(:,1);
        format = true(size(fd));
        format(fd > fd_thresh) = false;
        format(1:2) = false; % pre-steady state
                     
        %% Do the lagged correlation/covariance computation of TD matrices
        % CG: This has been simplified now that windowing is gone
        % So move temporal mask for each shift of the time series
        % no minimum for now... we'll see later (Anish has used minimum of 60 seconds, Ryan has used 20 seconds)
        % analysis done per session currently, then averaged over sessions
        % Only keep blocks >= min_frames
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
        
        [pl,~] = parabolic_interp_v2(Cov,tr); % peak lag = 333 x 333 matrix (time delay)

        pl(abs(pl) > lag_lim) = nan; % Long lags excluded [means they didn't estimate well]
        
        % running sum for averaging across sessions
        temp_lags = single(zeros(size(subj_lags)));
        temp_ZL = single(zeros(size(subj_lags)));
        temp_lags(good1,good2) = vec_mat_flip(pl,'TD'); % double check that this format works
        temp_ZL(good1,good2) = vec_mat_flip(Cov(:,lags==0),'ZL');
        
        subj_lags = cat(3,subj_lags,temp_lags);
        subj_lags = nansum(subj_lags,3);
        subj_ZL = cat(3,subj_ZL,temp_ZL);
        subj_ZL = nansum(subj_ZL,3);
        
        % running sum of nans
        subj_lags_nans = subj_lags_nans + isnan(temp_lags);
        subj_ZL_nans = subj_ZL_nans + isnan(temp_ZL); 

        % save out upper triangulars
        temp_lags = vec_mat_flip(temp_lags);
        temp_ZL = vec_mat_flip(temp_ZL);
%         save([outdir '/' subj '_' num2str(p) '_TD.mat'],'temp_lags');
%         save([outdir '/' subj '_' num2str(p) '_ZL.mat'],'temp_ZL');

        toc
    end
    
    % take average
    subj_lags_mean = subj_lags ./ (numel(patids) - subj_lags_nans);
    subj_ZL_mean = subj_ZL ./ (numel(patids) - subj_ZL_nans);
    
    % covariance matrix to correlation matrix
    d = zeros(size(subj_ZL_mean));
    d(logical(eye(length(subj_ZL_mean)))) = sqrt(diag(subj_ZL_mean));
    subj_ZL_mean = d^(-1)*subj_ZL_mean/d;
    
    % weighted lag projection for subject
    lag_weights = tan((pi/2)*(1-abs(subj_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
    lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
    lag_weights(isnan(subj_lags_mean)) = nan;
    subj_lags_mean_wghtd = subj_lags_mean.*lag_weights;
    subj_LP = nansum(subj_lags_mean_wghtd)./nansum(lag_weights); 

    % save lag projection
    % LP: CG - lag projection (average time delay relative to the rest of the
    % brain)..  low dimensional, loses a lot of information since different
    % lag threads are exclusionary
    % BUT lets you get a quick look if things look normal
    % Ryan will change this section of the code
    % now lag projections are weighted toward things that have high
    % correlations [at zero lag]
%     cifti.data = zeros(size(cifti.data,1),1);
%     cifti.data(mask) = subj_lags_proj';
%     cifti.hdr.dim(6) = 1;
%     ft_write_cifti_mod([outdir '/' subj '_LP_nosm.dtseries.nii'],cifti);
%     system([wb_dir '/wb_command -cifti-smoothing ' outdir '/' subj '_LP_nosm.dtseries.nii 3 3 COLUMN ' outdir '/' subj ...
%         '_LP_3mm.dtseries.nii -left-surface /data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/downsampled/' subj '.L.inflated.4k.surf.gii' ...
%         ' -right-surface /data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/downsampled/' subj '.R.inflated.4k.surf.gii']);
    
    % running sum for averaging across subjects
    grp_lags = cat(3,grp_lags,subj_lags_mean);
    grp_lags = nansum(grp_lags,3);
    grp_ZL = cat(3,grp_ZL,subj_ZL_mean);
    grp_ZL = nansum(grp_ZL,3);
    
    % running sum of nans
    grp_lags_nans = grp_lags_nans + isnan(subj_lags_mean);
    grp_ZL_nans = grp_ZL_nans + isnan(subj_ZL_mean);
    
end

% take average
grp_lags_mean = grp_lags ./ (numel(subjects) - grp_lags_nans);
grp_ZL_mean = grp_ZL ./ (numel(subjects) - grp_ZL_nans);
clear pl Cov covlags subj_* temp*

% weighted lag projection for group
lag_weights = tan((pi/2)*(1-abs(grp_ZL_mean))).^(-2);    % weighted by 1/f^2(r); f(r) = tan[(pi/2)(1-|r|)]
lag_weights(logical(eye(size(lag_weights)))) = 0;   % zero-out diagonal weights
lag_weights(isnan(grp_lags_mean)) = nan;
grp_lags_mean_wghtd = grp_lags_mean.*lag_weights;
grp_LP = nansum(grp_lags_mean_wghtd)./nansum(lag_weights);

% save lag projection
% cifti = ft_read_cifti_mod(['/data/nil-bluearc/raichle/ryan/lags/MSC/surface_lags/' subj '/cifti_timeseries_normalwall_atlas_4k/' patid '_LR_surf_subcort_333_4k_fsLR_smooth3.4_modesubcort.dtseries.nii']);
% cifti.data(mask) = grp_lags_proj';
% ft_write_cifti_mod([outdir '/MSC_LP_nosm.dtseries.nii'],cifti);
% system([wb_dir '/wb_command -cifti-smoothing ' outdir '/MSC_LP_nosm.dtseries.nii 3 3 COLUMN ' outdir '/MSC_LP_3mm.dtseries.nii ' ...
%     '-left-surface /data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/downsampled/Conte69.L.inflated.4k_fs_LR.surf.gii ' ...
%     '-right-surface /data/nil-bluearc/GMT/Evan/MSC/Analysis_V1/downsampled/Conte69.R.inflated.4k_fs_LR.surf.gii']);

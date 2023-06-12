%%% Script for playing around with lag results
topdir = '/data/nil-bluearc/GMT/Caterina/lags/';
outdir = [topdir 'comparisons/'];

% some things to set up
atlas_params = atlas_parameters('Parcels','/data/cn5/caterina/TaskConn_Methods/all_data/');
watershed_L = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_L_watershedmerge_0.35_tweaked.func.gii'];
watershed_R = ['/data/cn5/caterina/Atlases/Evan_parcellation/120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'];

subject_list = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC08','MSC09','MSC10'};
subjects = 1:10;
tasks = {'rest','motor','mixed','mem'};
maskmat = ones(333);
maskmat = logical(triu(maskmat,1));

%%% Load data (from lags_tdmx_parcels.m script)
for t = 1:length(tasks)
    cd([topdir tasks{t} '/']);
    for s = 1:10
        sub = sprintf('MSC%02d',s);
        subdata.(tasks{t})(s) = load([sub '_lagmats_mean.mat']);
        subsessdata.(tasks{t})(s) = load([sub '_lagmats_persession.mat']);
        lpmapdata.(tasks{t})(s,:) = subdata.(tasks{t})(s).subj_LP;
    end
end
cd(topdir);


%%% Also, create some z-scored differences comparing lags from each day
for t = 2:length(tasks) %excluding rest
    [h p ci stats] = ttest(lpmapdata.(tasks{t}),lpmapdata.rest);
    p_FDR = FDR(p);
    assign_data_to_parcel_cifti(stats.tstat',watershed_L,watershed_R,outdir,['Group_laprojection_' tasks{t} '_minus_rest_tstat']);
    assign_data_to_parcel_cifti(p_FDR',watershed_L,watershed_R,outdir,['Group_laprojection_' tasks{t} '_minus_rest_pFDR']);
end


%%% Create some interesting matrices
% sub x task lags
count = 1;
for s = 1:10
    for t = 1:length(tasks)
        alltasksub_lags_vec(count,:) = subdata.(tasks{t})(s).subj_lags_mean(maskmat);
        alltasksub_LP_vec(count,:) = subdata.(tasks{t})(s).subj_LP;
        alltasksub_ZL_vec(count,:) = subdata.(tasks{t})(s).subj_ZL_mean(maskmat);
        alltasksub_task(count) = t;
        alltasksub_sub(count) = s;
        count = count+1;        
    end
end

% sub x task x sess lags
count = 1;
for s = 1:10
    for t = 1:length(tasks)
        for d = 1:10 % day
             tmp = subsessdata.(tasks{t})(s).subj_lags(:,:,d);
             allsesstasksub_lags_vec(count,:) = tmp(maskmat);
             allsesstasksub_task(count) = t;
             allsesstasksub_day(count) = d;
             allsesstasksub_sub(count) = s;
             count = count+1;             
        end
    end
end

% sub x task x sess lags - SPLIT HALF
count = 1;
for s = 1:10
    for t = 1:length(tasks)
        
        % half one
        tmp = nanmean(subsessdata.(tasks{t})(s).subj_lags(:,:,1:5),3);
        allhalftasksub_lags_vec(count,:) = tmp(maskmat);
        allhalftasksub_task(count) = t;
        allhalftasksub_day(count) = 1;
        allhalftasksub_sub(count) = s;
        count = count+1;             
        
        % half one
        tmp = nanmean(subsessdata.(tasks{t})(s).subj_lags(:,:,6:10),3);
        allhalftasksub_lags_vec(count,:) = tmp(maskmat);
        allhalftasksub_task(count) = t;
        allhalftasksub_day(count) = 2;
        allhalftasksub_sub(count) = s;
        count = count+1;                     
    end
end


% Corr magnitude mask
mean_ZL = mean(alltasksub_ZL_vec,1);
mean_ZL_mat = zeros(333);
mean_ZL_mat(maskmat) = mean_ZL;
mean_ZL_mat = mean_ZL_mat + mean_ZL_mat';
%figure; imagesc(mean_ZL_mat);
maskmat2 = abs(mean_ZL) > 0.4;


%%% Make some figures
figure;
corrvals = nancov(alltasksub_lags_vec')./sqrt((nanvar(alltasksub_lags_vec')'*nanvar(alltasksub_lags_vec')));
imagesc(corrvals,[0 0.5]);
colorbar();
axis square;
title('Lags Similarity');
figure;
corrvals = nancov(alltasksub_lags_vec(:,maskmat2)')./sqrt((nanvar(alltasksub_lags_vec(:,maskmat2)')'*nanvar(alltasksub_lags_vec(:,maskmat2)')));
imagesc(corrvals,[0 0.5]);
colorbar();
axis square;
title('Lags Similarity, masked > 0.4');

figure;
corrvals = nancov(alltasksub_LP_vec')./sqrt((nanvar(alltasksub_LP_vec')'*nanvar(alltasksub_LP_vec')));
imagesc(corrvals,[0 0.5]);
colorbar();
axis square;
title('Lag Projection Similarity');

figure;
corrvals = nancov(allsesstasksub_lags_vec')./sqrt((nanvar(allsesstasksub_lags_vec')'*nanvar(allsesstasksub_lags_vec')));
imagesc(corrvals,[0 0.2]);
colorbar();
axis square;
title('Lags Similarity, SubxTaskxSess');

figure;
corrvals = nancov(allhalftasksub_lags_vec')./sqrt((nanvar(allhalftasksub_lags_vec')'*nanvar(allhalftasksub_lags_vec')));
imagesc(corrvals,[0 0.2]);
colorbar();
axis square;
title('Lags Similarity, SubxTaskxHalf');
figure;
corrvals = nancov(allhalftasksub_lags_vec(:,maskmat2)')./sqrt((nanvar(allhalftasksub_lags_vec(:,maskmat2)')'*nanvar(allhalftasksub_lags_vec(:,maskmat2)')));
imagesc(corrvals,[0 0.5]);
colorbar();
axis square;
title('Lags Similarity, SubxTaskxHalf, masked>0.2');

[aa ind] = sort(allhalftasksub_task);
figure;
corrvals = nancov(allhalftasksub_lags_vec(ind,maskmat2)')./sqrt((nanvar(allhalftasksub_lags_vec(ind,maskmat2)')'*nanvar(allhalftasksub_lags_vec(ind,maskmat2)')));
imagesc(corrvals,[0 0.5]);
colorbar();
axis square;
title('Lags Similarity, SubxTaskxHalf, by task, masked>0.2');
[simvals.raw subvals.raw condvals.raw] = matrix_similarity_calcs(corrvals,allhalftasksub_task,allhalftasksub_sub,allhalftasksub_day);
matrix_similarity_quantification_fig(subvals,'orig','corr',subjects,0,0);

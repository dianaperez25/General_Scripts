%% make fake timeseries for figure


% load some subjects timeseries
%timeseries = ft_read_cifti_mod('/Volumes/fsmresfiles/PBS/Gratton_Lab/iNetworks/BIDS/Nifti/derivatives/postFCproc_CIFTI_20.2.0/sub-INET001/ses-1/cifti_timeseries_normalwall/sub-INET001_ses-1_task-rest_run-1_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii')
timeseries = load('/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/postFCproc_CIFTI/FC_Parcels_333/sub-INET001_rest_ses-1_parcel_timecourse.mat');
% apply t-mask?
thresh = 0.05;
bins = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, thresh];
mod_threshold = 0.05;
fake_signals = {};

%% generate fake data for each subject, calculate mod and fc 
for i = 1:length(bins)
    disp(i) 
    count = 0;
    

            bin = bins(i);
%             load(functionals_path)
%             parcel_timeseries = parcel_time{session};
%             tmask = tmask_all{session};
            %timeseries = parcel_timeseries(tmask==1, :)';

%             if(size(timeseries,2) <= 333)
%                 fc_sim(i,j,k) = NaN;
%                 mod_sim(i,j,k) = NaN;
%                 continue;
%             end

            real_signal = timeseries.parcel_time;
            real_fc = paircorr_mod(real_signal);
            rand_signal = randn(size(real_signal));
            
            
            % make sure that the matrix is positive-determinant
            try
               R = chol(real_fc);
            catch
               fc_sim(i) = NaN;
               mod_sim(i) = NaN;
               %continue;
            end  
            
            % get eigenvectors of real FC 
            [V D] = eig(real_fc);
            V_use = V*sqrt(D);

            % apply them to our random signal calculate coflux, modularity
            % etc
            fake_signal = V_use*rand_signal;
            fake_signals{j,k} = fake_signal;
           

end
load('/Volumes/fsmresfiles/PBS/Gratton_Lab/Lifespan/Post-COVID/BIDS/derivatives/postFCproc_CIFTI/FC_Parcels_333/sub-INET001_rest_ses-1_parcel_timecourse.mat')
real_signal = parcel_time;
addpath '/Users/dianaperez/Documents/GitHub/Lifespan-Analysis/utilities'
real_fc = paircorr_mod(real_signal);
rand_signal = randn(size(real_signal));
try
    R = chol(real_fc);
catch
    fc_sim(i) = NaN;
    mod_sim(i) = NaN;
end
[V D] = eig(real_fc);
V_use = V*sqrt(D);
fake_signal = V_use*rand_signal';
fake_fc = corr(fake_signal');
% find correlation values that are both low and high
smooth_sig_1 = smoothdata(fake_signal(15, 1:400));
smooth_sig_1 = smoothdata(smooth_sig_1);
smooth_sig_1 = smoothdata(smooth_sig_1);
smooth_sig_1 = smoothdata(smooth_sig_1);

smooth_sig_2 = smoothdata(fake_signal(137, 1:400));
smooth_sig_2 = smoothdata(fake_signal(8, 1:400));
smooth_sig_2 = smoothdata(smooth_sig_2);
smooth_sig_2 = smoothdata(smooth_sig_2);
smooth_sig_2 = smoothdata(smooth_sig_2);

% this is negative corr with sig 1 -- don't use
smooth_sig_3 = smoothdata(fake_signal(147, 1:400));
smooth_sig_3 = smoothdata(smooth_sig_3);
smooth_sig_3 = smoothdata(smooth_sig_3);
smooth_sig_3 = smoothdata(smooth_sig_3);

% positive low corr with sig 1
smooth_sig_4 = smoothdata(fake_signal(14, 1:400));
smooth_sig_4 = smoothdata(smooth_sig_4);
smooth_sig_4 = smoothdata(smooth_sig_4);
smooth_sig_4 = smoothdata(smooth_sig_4);

% figure 1, high correlation
figure;
% orange
plot(smooth_sig_1, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980])
hold on
% blue
plot(smooth_sig_2, 'LineWidth', 3, 'Color', [0 0.4470 0.7410])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

% % figure 2; low (negative) correlation
% figure;
% % orange
% plot(smooth_sig_1, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980])
% hold on
% % green
% plot(smooth_sig_3, 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880])
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);


% figure 2; low correlation
figure;
% orange
plot(smooth_sig_1, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980])
hold on
% green
plot(smooth_sig_4, 'LineWidth', 3, 'Color', [0.4660 0.6740 0.1880])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
function make_spCorr_Map(cortexOnly, dconn_dat, outputname, groupAvg_path)

if cortexOnly
    voxnum = 59412;
end

if ischar(dconn_dat)
    dconn_dat = ft_read_cifti_mod(dconn_dat);
    cifti_corrmap = dconn_dat.data(1:voxnum,1:voxnum);
else
    cifti_corrmap = dconn_dat(1:voxnum, 1:voxnum);
end

clear dconn_dat

tic;
disp('making spatial correlations...')

cifti_corrmap(isnan(cifti_corrmap)) = 0;
cifti_corrmap1 = cifti_corrmap(1:voxnum,1:(voxnum/4));
disp(['First quarter of Cifti corrmap is ' num2str(size(cifti_corrmap1,1)) ' by ' num2str(size(cifti_corrmap1,2)) ': ' datestr(now)])
% Remove NaNs (produced if vertices have no data)


for dconnrow = 1:size(cifti_corrmap1,1)            
    cifti_corrmap1(dconnrow,:) = single(FisherTransform(cifti_corrmap1(dconnrow,:)));
end

disp(sprintf('First quarter of Fisher Transform Finished: %s', datestr(now)));
                    
disp(sprintf('Loading Template: %s', datestr(now)));

% Load group-average corrmat (assumes it's a dconn)
group = ft_read_cifti_mod(groupAvg_path);
group = group.data(1:voxnum,1:voxnum);
group1 = group(1:voxnum,1:(voxnum/4));
sizedata = size(group1,2);

if cortexOnly==1
    % Load template variants file
    template = ft_read_cifti_mod('/projects/b1081/member_directories/bkraus/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
    template.data = [];
else
    template = ft_read_cifti_mod('/projects/b1081/MSC/TaskFC/FCProc_MSC05_motor_pass2/cifti_timeseries_normalwall_native_freesurf/vc39006_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii');
    template.data = [];
    template.hdr(6).dim = 1;
end

disp(['First quarter of Cifti group data is ' num2str(size(group1,1)) ' by ' num2str(size(group1,2)) ': ' datestr(now)])

disp(sprintf('Template Loaded: %s', datestr(now)));
            
            % Compare to group average
for i=1:sizedata
    template.data(i,1) = paircorr_mod(group1(:,i),cifti_corrmap1(:,i));
end
            
clear cifti_corrmap1
clear group1

disp(sprintf('First quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap2 = cifti_corrmap(1:voxnum,(voxnum/4)+1:(voxnum/2));

disp(['Second quarter of Cifti corrmap is ' num2str(size(cifti_corrmap2,1)) ' by ' num2str(size(cifti_corrmap2,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap2,1)
    cifti_corrmap2(dconnrow,:) = single(FisherTransform(cifti_corrmap2(dconnrow,:)));
end

disp(sprintf('Second quarter of Fisher Transform Finished: %s', datestr(now)));

group2 = group(1:voxnum,(voxnum/4)+1:(voxnum/2));
sizedata = size(group2,2);

disp(['Second quarter of Cifti group data is ' num2str(size(group2,1)) ' by ' num2str(size(group2,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data((i+(voxnum/4)),1) = paircorr_mod(group2(1:voxnum,i),cifti_corrmap2(1:voxnum,i));
end

clear cifti_corrmap2
clear group2

disp(sprintf('Second quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap3 = cifti_corrmap(:,(voxnum/2)+1:(voxnum*.75));

disp(['Third quarter of Cifti corrmap is ' num2str(size(cifti_corrmap3,1)) ' by ' num2str(size(cifti_corrmap3,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap3,1)
    cifti_corrmap3(dconnrow,:) = single(FisherTransform(cifti_corrmap3(dconnrow,:)));
end

disp(sprintf('Third quarter of Fisher Transform Finished: %s', datestr(now)));

group3 = group(:,(voxnum/2)+1:(voxnum*.75));
sizedata = size(group3,2);

disp(['Third quarter of Cifti group data is ' num2str(size(group3,1)) ' by ' num2str(size(group3,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data(i+(voxnum/2),1) = paircorr_mod(group3(:,i),cifti_corrmap3(:,i));
end

clear cifti_corrmap3
clear group3

disp(sprintf('Third quarter of Correlation Finished: %s', datestr(now)));

cifti_corrmap4 = cifti_corrmap(:,(voxnum*.75)+1:end);

disp(['Fourth quarter of Cifti corrmap is ' num2str(size(cifti_corrmap4,1)) ' by ' num2str(size(cifti_corrmap4,2)) ': ' datestr(now)])

for dconnrow = 1:size(cifti_corrmap4,1)
    cifti_corrmap4(dconnrow,:) = single(FisherTransform(cifti_corrmap4(dconnrow,:)));
end

disp(sprintf('Fourth quarter of Fisher Transform Finished: %s', datestr(now)));

group4 = group(:,(voxnum*.75)+1:end);
sizedata = size(group4,2);

disp(['Fourth quarter of Cifti group data is ' num2str(size(group4,1)) ' by ' num2str(size(group4,2)) ': ' datestr(now)])

for i=1:sizedata
    template.data(i+(voxnum*.75),1) = paircorr_mod(group4(:,i),cifti_corrmap4(:,i));
end

clear cifti_corrmap4
clear group4

disp(sprintf('Fourth quarter of Correlation Finished: %s', datestr(now)));


disp(['Final size of spatial correlation map is ' num2str(size(template.data,1)) ' by ' num2str(size(template.data,2)) ': ' datestr(now)])


% Write out the variants

if cortexOnly == 1

    ft_write_cifti_mod([outputdir '/' outputname '_vs_' groupAvgName '_cortex_corr'],template)
    template.data = [];

elseif cortexOnly == 0

    ft_write_cifti_mod([outputdir '/' outputname '_vs_' groupAvgName '_subcort_corr'],template)
    template.data = [];

end
 
disp('done')
toc
end
 

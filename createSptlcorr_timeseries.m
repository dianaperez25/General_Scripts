function createSptlcorr_MSCdconns_timeseries(groupAvgLoc,groupAvgName,cortexOnly,outputdir,dconnData, outputname)

%addpath(genpath('/data/cn/data1/scripts/CIFTI_RELATED/Resources/'))

%clear global template  %% Removes global dconn correlation matrix from memory
%clear global template2  

if ~exist(outputdir , 'dir')
mkdir(outputdir)
end


   
    
    % Compute spatial correlations (split into quarters)
    
            cifti_timeseries = dconnData;
            
            %cifti_corrmap = dconnData;
            clear dconnData
            
            cifti_corrmap = paircorr_mod(cifti_timeseries');
            
            cifti_corrmap = cifti_corrmap(:,1:(59412/4));

            disp(['First quarter of Cifti corrmap is ' num2str(size(cifti_corrmap,1)) ' by ' num2str(size(cifti_corrmap,2)) ': ' datestr(now)])

            % Remove NaNs (produced if vertices have no data)
            cifti_corrmap(isnan(cifti_corrmap)) = 0;

            % Apply the Fisher tranformation
            
            %cifti_corrmap = single(FisherTransform(cifti_corrmap));    %%Changed to row-wise transformation for memory allocation
            
            
            for dconnrow = 1:size(cifti_corrmap,1)
            
                cifti_corrmap(dconnrow,:) = single(FisherTransform(cifti_corrmap(dconnrow,:)));
                
            end
            
          

            disp(sprintf('First quarter of Fisher Transform Finished: %s', datestr(now)));
            
            
disp(sprintf('Loading Template: %s', datestr(now)));

% Load group-average corrmat (assumes it's a dconn)
group = ft_read_cifti_mod([groupAvgLoc '/' groupAvgName '.dconn.nii']);
if cortexOnly==1
    group = group.data(1:59412,1:59412);
    group = group(:,1:(59412/4));
    sizedata = size(group,2);
    
    % Load template variants file
    template = ft_read_cifti_mod('/projects/b1081/member_directories/bkraus/Brian_MSC/dconn_scripts/templates/MSC01_allses_mean_native_freesurf_vs_120sub_corr.dtseries.nii');
    template.data = [];
else
    group = group.data(1:65625,1:65625);
    sizedata = size(group,1);
    
    %%%%%%% NEED TO MAKE YOUR A TEMPLATE FILE SO DIMENSIONS MATCH %%%%%%%%
    template = ft_read_cifti_mod('/projects/b1081/MSC/TaskFC/FCProc_MSC05_motor_pass2/cifti_timeseries_normalwall_native_freesurf/vc39006_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii');
    template.data = [];
    template.hdr(6).dim = 1;
    
end

disp(['First quarter of Cifti group data is ' num2str(size(group,1)) ' by ' num2str(size(group,2)) ': ' datestr(now)])

disp(sprintf('Template Loaded: %s', datestr(now)));
            
            % Compare to group average
            for i=1:sizedata
                template.data(i,1) = paircorr_mod(group(:,i),cifti_corrmap(:,i));
            end
            
            clear cifti_corrmap
            clear group
            
            disp(sprintf('First quarter of Correlation Finished: %s', datestr(now)));
            
            cifti_corrmap = paircorr_mod(cifti_timeseries');
            
            %clear cifti_timeseries
            
            cifti_corrmap = cifti_corrmap(:,(59412/4)+1:(59412/2));

            disp(['Second quarter of Cifti corrmap is ' num2str(size(cifti_corrmap,1)) ' by ' num2str(size(cifti_corrmap,2)) ': ' datestr(now)])
            
            for dconnrow = 1:size(cifti_corrmap,1)
            
                cifti_corrmap(dconnrow,:) = single(FisherTransform(cifti_corrmap(dconnrow,:)));
                
            end
            
            disp(sprintf('Second quarter of Fisher Transform Finished: %s', datestr(now)));
            
            group = ft_read_cifti_mod([groupAvgLoc '/' groupAvgName '.dconn.nii']);
            
            group = group.data(1:59412,1:59412);
            group = group(:,(59412/4)+1:(59412/2));
            sizedata = size(group,2);
            
            disp(['Second quarter of Cifti group data is ' num2str(size(group,1)) ' by ' num2str(size(group,2)) ': ' datestr(now)])

            for i=1:sizedata
                template.data(i+(59412/4),1) = paircorr_mod(group(:,i),cifti_corrmap(:,i));
            end
            
            clear cifti_corrmap
            clear group
            
            disp(sprintf('Second quarter of Correlation Finished: %s', datestr(now)));
            
            cifti_corrmap = paircorr_mod(cifti_timeseries');
            
            cifti_corrmap = cifti_corrmap(:,(59412/2)+1:(59412*.75));

            disp(['Third quarter of Cifti corrmap is ' num2str(size(cifti_corrmap,1)) ' by ' num2str(size(cifti_corrmap,2)) ': ' datestr(now)])
            
            for dconnrow = 1:size(cifti_corrmap,1)
            
                cifti_corrmap(dconnrow,:) = single(FisherTransform(cifti_corrmap(dconnrow,:)));
                
            end
            
            disp(sprintf('Third quarter of Fisher Transform Finished: %s', datestr(now)));
            
            group = ft_read_cifti_mod([groupAvgLoc '/' groupAvgName '.dconn.nii']);
            
            group = group.data(1:59412,1:59412);
            group = group(:,(59412/2)+1:(59412*.75));
            sizedata = size(group,2);
            
            disp(['Third quarter of Cifti group data is ' num2str(size(group,1)) ' by ' num2str(size(group,2)) ': ' datestr(now)])

            for i=1:sizedata
                template.data(i+(59412/2),1) = paircorr_mod(group(:,i),cifti_corrmap(:,i));
            end
            
            clear cifti_corrmap
            clear group
            
            disp(sprintf('Third quarter of Correlation Finished: %s', datestr(now)));
            
            cifti_corrmap = paircorr_mod(cifti_timeseries');
            
            clear cifti_timeseries
            
            cifti_corrmap = cifti_corrmap(:,(59412*.75)+1:end);

            disp(['Fourth quarter of Cifti corrmap is ' num2str(size(cifti_corrmap,1)) ' by ' num2str(size(cifti_corrmap,2)) ': ' datestr(now)])
            
            for dconnrow = 1:size(cifti_corrmap,1)
            
                cifti_corrmap(dconnrow,:) = single(FisherTransform(cifti_corrmap(dconnrow,:)));
                
            end
            
            disp(sprintf('Fourth quarter of Fisher Transform Finished: %s', datestr(now)));
            
            group = ft_read_cifti_mod([groupAvgLoc '/' groupAvgName '.dconn.nii']);
            
            group = group.data(1:59412,1:59412);
            group = group(:,(59412*.75)+1:end);
            sizedata = size(group,2);
            
            disp(['Fourth quarter of Cifti group data is ' num2str(size(group,1)) ' by ' num2str(size(group,2)) ': ' datestr(now)])

            for i=1:sizedata
                template.data(i+(59412*.75),1) = paircorr_mod(group(:,i),cifti_corrmap(:,i));
            end
            
            clear cifti_corrmap
            clear group
            
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
 


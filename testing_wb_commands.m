%% TRYING OUT SOME WB COMMANDS

workbenchdir = '/projects/b1081/Scripts/workbench2/bin_linux64/';
%surface = '/projects/p31161/lateralizationVariants/113619_ThresholdedVariantMap_SNRExclude_7.5.dtseries.nii';
surface = '/projects/p31161/lateralizationVariants/140117.R.midthickness.32k_fs_LR.surf.gii';
%surface = '/projects/b1081/HCP_analyses/from_Ben/HCP_variants/spCorr/100206_allRuns_mean_vs_120_avg_corr_LR_corr.dtseries.nii'; 
%outfile_L = '/projects/p31161/lateralizationVariants/100206_Testing_cifti_separate_left.func.gii';
%outfile_R = '/projects/p31161/lateralizationVariants/100206_Testing_cifti_separate_right.func.gii';
outfile = '/projects/p31161/lateralizationVariants/140117_Testing_surface_to_metric_R.func.gii';
system([workbenchdir 'wb_command -surface-coordinates-to-metric ' surface ' ' outfile])
%system([workbenchdir 'wb_command -surface-flip-lr ' surface ' ' outfile])
%system([workbenchdir 'wb_command -cifti-separate ' surface ' COLUMN -metric CORTEX ' outfile])

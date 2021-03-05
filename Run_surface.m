clear all

PostFreeSurferPipeline_fsavg2fslr_long_GrattonLab('LS05','/projects/b1081/Lifespan/derivatives/freesurfer-6.0.1/', '/projects/b1081/Lifespan/derivatives/preproc_fmriprep-20.0.6/fmriprep/sub-LS05/anat/');

PostFreeSurferPipeline_fsavg2fslr_long_GrattonLab('LS07','/projects/b1081/Lifespan/derivatives/freesurfer-6.0.1/', '/projects/b1081/Lifespan/derivatives/preproc_fmriprep-20.0.6/fmriprep/sub-LS07/anat/');

post_fc_processing_batch_GrattonLab('post_fc_processing_batch_params_iNetworks.m');
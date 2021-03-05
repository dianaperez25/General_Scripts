%% Script to run post fc surface processing
addpath(genpath('/projects/b1081/Scripts'))
addpath(genpath('/projects/p31161/Scripts'))

cd /projects/p31161/Repositories/GrattonLab-General-Repo/SurfacePipeline

%post_fc_processing_batch_GrattonLab('post_fc_processing_batch_params_iNetworks.m')
make_dconn_session_wrapper()
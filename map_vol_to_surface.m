function map_vol_to_surface_MSCspecific(volfile,subject,task)
% volfile should be without extension

%volfile = '1_L_Hand+2_R_Hand-3_L_Leg-4_R_Leg_zstat_333_t88';

% path to goodvoxels - use the one from the first day. For this purpose
% (just mapping), just use the first one in the folder. Could make a union
% mask eventually to be more accurate
goodvoxdir = ['/data/nil-bluearc/GMT/Caterina/SurfTask_analysis/' subject '/' task '_vol/goodvoxels/'];
goodvox_files = dir([goodvoxdir '*goodvoxels.nii.gz']);
goodvoxmask = [goodvoxdir goodvox_files(1).name];

surfdir = '/projects/b1081/Lifespan/derivatives/freesurfer-6.0.1/FREESURFER_fs_LR/';
timdir = '/data/nil-bluearc/GMT/Laumann/MSC';
maskdir = [timdir '/' subject '/subcortical_mask_native_freesurf'];
workbenchdir = '/data/cn5/caterina/workbench/bin_linux64/';



% first, make a nifti version fo the file
% Note: annoyingly have to adjust some FIDL output files to have the center
% coordinate information
prep_nifti_file(volfile)

% resample to the surface using the subject's own pial/white/etc. data
HEMS = {'L','R'};
for hem = 1:2

    midsurf = [surfdir '/sub-' subject '/NativeVol/Native/sub-' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
    midsurf_LR32k = [surfdir '/sub-' subject '/NativeVol/fsaverage_LR32k/sub-' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    whitesurf = [surfdir '/sub-' subject '/NativeVol/Native/sub-' subject '.' HEMS{hem} '.white.native.surf.gii'];
    pialsurf = [surfdir '/sub-' subject '/NativeVol/Native/sub-' subject '.' HEMS{hem} '.pial.native.surf.gii'];
    nativedefsphere = [surfdir '/sub-' subject '/NativeVol/Native/sub-' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
    outsphere = [surfdir '/sub-' subject '/NativeVol/fsaverage_LR32k/sub-' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];

    system([workbenchdir '/wb_command -volume-to-surface-mapping ' volfile '.nii.gz ' midsurf ' ' volfile '_' HEMS{hem} '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' goodvoxmask]);
    system([workbenchdir '/wb_command -metric-dilate ' volfile '_' HEMS{hem} '.func.gii ' midsurf ' 10 ' volfile '_' HEMS{hem} '_dil10.func.gii']);
    system([workbenchdir '/wb_command -metric-resample ' volfile '_' HEMS{hem} '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' volfile '_' HEMS{hem} '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
    
    % remove intermediate files
    delete([volfile '_' HEMS{hem} '.func.gii']);
    delete([volfile '_' HEMS{hem} '_dil10.func.gii']);

end

% put together into a cifti
maskdir = [timdir '/' subject '/subcortical_mask_native_freesurf'];
left_mask = ['/data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii'];
right_mask = ['/data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii'];

timename_L = [volfile '_L_dil10_32k_fs_LR.func.gii'];
timename_R = [volfile '_R_dil10_32k_fs_LR.func.gii'];
system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_L]);
system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_R]);
system([workbenchdir '/wb_command -cifti-create-dense-scalar ' volfile '_LR_surf_subcort_333_32k_fsLR.dscalar.nii -volume ' volfile '.nii.gz ' maskdir '/subcortical_mask_LR_333.nii -left-metric ' timename_L ' -roi-left ' left_mask ' -right-metric ' timename_R ' -roi-right ' right_mask ]);

% remove intermediate files
delete([volfile '_L_dil10_32k_fs_LR.func.gii']);
delete([volfile '_R_dil10_32k_fs_LR.func.gii']);

end

function prep_nifti_file(volfile)

% check to see if the center coordinates are present
% if not, add them along with mmpix (assumed to come together)
[s result] = system(['grep "center" ' volfile '.4dfp.ifh']); 
if isempty(result) %missing center information
    system(['cp ' volfile '.4dfp.ifh ' volfile '_ORIG.4dfp.ifh']); %copy over original for safekeeping
    fid = fopen([volfile '.4dfp.ifh'],'a');
    fprintf(fid,'center                            := 73.500000 -87.000000 -84.000000\n');
    fprintf(fid,'mmppix                            := 3.000000 -3.000000 -3.000000\n');
    fclose(fid);
end

system(['niftigz_4dfp -n ' volfile ' ' volfile]);

end
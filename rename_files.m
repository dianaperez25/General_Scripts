%script to rename fmriprep output files

d = '/projects/b1081/Lifespan/derivatives/preproc_fmriprep-20.2.0/fmriprep/sub-LS03/ses-4/func/FD_outputs/'
dirinfo = dir(d);
oldnames = {dirinfo.name};
oldnames = oldnames(15:end);
newnames = {};

for n = 1:length(oldnames)
    newname = strrep(oldnames(n),'run-0','run-');
    newnames = [newnames newname];
end

for f = 1:length(oldnames)
movefile(fullfile(d,oldnames{f}),fullfile(d,newnames{f}))
end

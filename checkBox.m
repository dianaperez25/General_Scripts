clear all
%change directory to RDSS folder
cd /Volumes/GRATTONLAB
addpath ('/Volumes/GRATTONLAB')
RDSSdirs = dir('HCP');
RDSSdirNames = RDSSdirs(:).name;

%change directory to Box folder
cd /Users/dianaperez/Box/DATA
%add folder to path
addpath ('/Users/dianaperez/Box/DATA')
%get the names of all the directories in that folder in a structure
BoxdirNames = dir('HCP');
%add all of those directories to the path
for x = 1:length(BoxdirNames)
 if ~ismember(BoxdirNames(x).name, 
 end

%setting up variables
subjID = {};
output = {};
subj = {};
count = 1;
count2 = 1;
dirs = {'rfMRI_REST1_RL', 'rfMRI_REST1_LR', 'rfMRI_REST2_RL', 'rfMRI_REST2_LR'};

for x = 1:length(dirnames)
     subj{x} = num2str(extractBefore(dirnames(x).name, 7));
    if ~ismember(subj{x}, subjID)
        subjID{count} = subj{x};
        count = count + 1;
    end
end

subjID = subjID';

for x = 1:length(subjID)
    for y = length(dirs)
        if ~isfolder([subjID{x} '_' dirs{y}])
            output{count2} = ['subject ' subjID{x} ' is missing directory ' dirs{y}];
            count2 = count2 + 1;
        end
    end
end
    
%these are the files that should be in every subject's directory
HPCfiles = {'fc.fcparams', 'Movement_Regressors.txt', [dirs{s} '.4dfp.hdr'], [dirs{s} '.4dfp.ifh'],[dirs{s} '.4dfp.img'],...
    [dirs{s} '.4dfp.img.rec'], [dirs{s} '.ddat'], [dirs{s} '_dsd0.4dfp.hdr'], [dirs{s} '_dsd0.4dfp.ifh'], [dirs{s} '_dsd0.4dfp.img'],...
    [dirs{s} '_dsd0.4dfp.img.rec'], [dirs{s} '_func_vols.conc.rec'], [dirs{s} '_func_vols.crit'], [dirs{s} '_func_vols.dat'],...
    [dirs{s} '_func_vols.lst'], [dirs{s} '_func_vols.vals'], [dirs{s} '_func_vols.xmgr'], [dirs{s} '_func_vols_ave.4dfp.hdr'], [dirs{s} '_func_vols_ave.4dfp.ifh'],...
    [dirs{s} '_func_vols_ave.4dfp.img'], [dirs{s} '_func_vols_ave.4dfp.img.rec'], 'T1w_restore_brain.4dfp.hrd', 'T1w_restore_brain.4dfp.ifh', 'T1w_restore_brain.4dfp.img',...
    'T1w_restore_brain.4dfp.img.rec', 'T1w_restore_brain.nii.gz', 'T1w_restore_brain_222.4dfp.hdr', 'T1w_restore_brain_222.4dfp.ifh', 'T1w_restore_brain_222.4dfp.img',...
     'T1w_restore_brain_222.4dfp.img.rec'};
    
for i = 1:length(dirnames)
    subj{i} = num2str(extractBefore(dirnames(i).name, 7));
    cd(dirnames(i).name)
    for j = 1:length(HPCfiles)
        
    end
    cd ..
        
        
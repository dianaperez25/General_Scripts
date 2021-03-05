numperms=1000;

% this script shuffles border-shift/ectopic labels within a subject's
% variants
% with each permutation, the labels are shuffled across 384 subjects, then
% overlap maps of all border-shift variants and all ectopic variants are
% created
% q - are distribution of border-shift and ectopic variants sig different?
% --> calculate spatial corr of each randomized bordershift overlap map and
% randomized ectopic overlap map


% load in subjects ('sublist')
load /data/cn5/seitzmanb/HCP/fcprocessing/topVariants/HCPfull/spCorr/thresholded/reassigned/correlationSplitHalf.mat
sublist= [splitHalf1(:,1); splitHalf2(:,1)]; clear split*;

% initialize mat of values i'm permuting
spCorrs=zeros(numperms,1);
diffmaps=zeros(59412,numperms);

% loop thru permutations;
for i=1:numperms
    
    % initializing overlap map that will result from this permutation
    % across all 384 subjects
    allSubsBordMap= zeros(59412,1);
    allSubsEcMap= zeros(59412,1);
    
    for s=1:length(sublist)
        
        subID=num2str(sublist(s));
        
        % read in the subjects variants labeled by UNIQUE ID (1 thru # variants)
        varIDs= ft_read_cifti_mod(['/data/cn6/allyd/BorderEctopic/HCP384/splitVars_uniqueIDs/reassigned/' subID '_uniqueIDs_afterReassign.dtseries.nii']);
        varIDs= varIDs.data;
        % read in the file which labels border-shift and ectopic variants differently
        bordEc= ft_read_cifti_mod(['/data/cn6/allyd/BorderEctopic/HCP384/splitVars_uniqueIDs/reassigned/borderEctopic/' subID '_border1ectopic2.dtseries.nii']);
        bordEc= bordEc.data;
        
        numvars= length(unique(varIDs))-1;
        
               
        % what are the subject's border/ectopic labels to shuffle?
        % (want to keep proportion of each type, but will shuffle ID labels)
        bordEcLabels= zeros(numvars,1);
        for vv=1:numvars
            varInds=find(varIDs==vv);
            bordEcLabels(vv,1)= bordEc(varInds(1));
        end
        clear varInds vv
        
        rng('shuffle'); % ensures different random seed for each permutation
        randVarIDs= randperm(numvars); % we will look through the shuffled IDs and assign them to fake border/ectopic labels

        
        % shuffle by random ID
        randBordEcLabels= zeros(numvars,1);
        randBordMap= zeros(59412,1);
        randEcMap= zeros(59412,1);
        for v=1:numvars
            randBordEcLabels(v,1)= bordEcLabels(randVarIDs(v));
            
            if randBordEcLabels(v,1)==1 % 1 means it is a border-shift variant
                randBordMap(varIDs==v) = 1; % subject's border-shift variants map gets a +1
            elseif randBordEcLabels(v,1)==2 % 2 means it is an ectopic variant
                randEcMap(varIDs==v) = 1; %subject's ectopic variants map gets a +1
            end
        end
        
        allSubsBordMap= allSubsBordMap + randBordMap; %add in this subject to the all-subjects overlap map
        allSubsEcMap = allSubsEcMap + randEcMap; %again for ectopic variants
        
        clear randBordMap randEcMap v 
    end
    
    % save the spatial corr value of overap maps for this permutation
    permDiffmat= allSubsEcMan - allSubsBordMap;
    diffmats(:,i)= permDiffmat;
    spCorrs(i,1)= corr(allSubsBordMap, allSubsEcMap);
    clear allSubsBordMap allSubsEcMap permDiffmat
    
end


save(['/data/cn6/allyd/BorderEctopic/permute/borderEctopic_spCorr_' num2str(numperms) 'permutations.mat'],'spCorrs')
save(['/data/cn6/allyd/BorderEctopic/permute/borderEctopic_diffMap_' num2str(numperms) 'permutations.mat'],'diffmats')


function [dice_summary_path, network_assignment_path,winners_summary_path]  = assignvars(bold_data_path, variants_file, outfilepath, outfilename)
%add back later
%network_assignment_path, 
%bold_data_path=timeseries data

%load templates- matrix of 16 networks x XX,XXX vertices
network_names = {'DMN'	'Vis'	'FP'	'Unassigned'	'DAN'	'Unassigned2'	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'	'PMN'	'PON'};
wb_colors = [1 2 3 5 7 8 9 10 11 12 13 14 15 16];

templates=load('/projects/b1081/Scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat');
%templates=load('/Users/Alexis/Box/Quest_Backup/Scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat');

templates=templates.templates(1:59412,:)';
wb_colors= [1 2 3 5 7 8 9 10 11 12 13 14 15 16];
cifti_geo = load('/Users/dianaperez/Box/Quest_Backup/member_directories/zladwig/dependencies/Important_Files/Cifti_geo_distances_uint8.mat');
%cifti_geo = load('/Users/Alexis/Box/Quest_Backup/member_directories/zladwig/Important_Files/Cifti_geo_distances_uint8.mat');
cifti_geo.distances = cifti_geo.distances(1:59412,1:59412);
%consensus = ft_read_cifti_mod('/projects/b1081/member_directories/zladwig/dependencies/Important_Files/120_colorassn_minsize400_manualconsensus.dtseries.nii');
consensus = ft_read_cifti_mod('/Users/dianaperez/Box/Quest_Backup/member_directories/zladwig/Important_Files/120_colorassn_minsize400_manualconsensus.dtseries.nii');
getmaxmatch = 0;    % Sets whether to get best match value for each variant
minsize = 50;       % Minimum size in vertices for each variant

if getmaxmatch == 1
    min_match_dice = 0;
    maxmatch = [];
else 
    min_match_dice=.014619; %this is number I will set
    maxmatch=[];
end


    %0.1432; % brian set (ALEXIS: this was set by the lowest 5% of matches, they get assigned to the 'Unassigned' network. You'll have to recalculate this


%threshold templates at top 5%
temp= sort(templates(:),'descend');
templateThresh= temp(round(0.05 * numel(temp)));
cifti_template_mat_full= templates > templateThresh;
clear temp

% get original bold data 
cifti = ft_read_cifti_mod(bold_data_path);
cifti.data=cifti.data(1:59412,:);
bold_data = cifti.data;

% get the variants for the person 
variants_map = ft_read_cifti_mod(variants_file);
unique_vars = unique(variants_map.data);
unique_vars(1) = [];    %% Gets rid of 0 values in cifti

% make some templates for where we will store reassigned variants and their dice 

varIDsreassign=variants_map; 
varIDsreassign.data=zeros(59412,1);

diceOfWinner= variants_map; 
diceOfWinner.data=zeros(59412,1);

% for each variant, find all of its vertices. Then, calculate the seedmap for each vertex 
dice_summary = [];
winners_summary = [];
for var=1:length(unique_vars)

    % find vertices
    vertex_idx = find(variants_map.data == unique_vars(var));
    
    if length(vertex_idx) >= minsize    %% If variant size > minsize

        % calculate timeseries for each vertex
        cifti_subset = bold_data(vertex_idx, :);
        
        % get mean timeseries
        mean_cifti_subset = mean(cifti_subset);
        
        % calculate a seedmap from all of those points
        varRmat = paircorr_mod(cifti_subset',bold_data');
        
        % old - calculate a mean seedmap from all of those points
        %varmean = paircorr_mod(mean_cifti_subset', bold_data');
        
        % get mean seedmap of all vertices in a variant
        varmean = mean(varRmat, 1);
        
        %% NEW remove all vertices which were within 30mm of ANY vertex in the variant
        
        excludedverts = [];
        
        for i=1:length(vertex_idx)
            excludedverts = [excludedverts find(cifti_geo.distances(vertex_idx(i),:) < 30)];
        end
        
        finalexclude = unique([vertex_idx' excludedverts]);
        
        varmean(1, finalexclude) = NaN;
        
        % only look at the top 5% most correlated vertices in the brain (ab 3000 vertices)
        temp= sort(varmean,'descend', 'MissingPlacement', 'last');
        varthresh= temp(round(0.05 * numel(temp)));
        varmean_thresh= varmean > varthresh;
        clear temp;
        
        %disp(['overlap was ' num2str(length(intersect(find(varmean_thresh ==1),finalexclude))) 'out of ' num2str(length(find(varmean_thresh ==1)))]);
        % loop through each network and compare how similar the network map is
        % to the thresholded mean seed map
        
        for net=1:14
            % dice calculation b/w variant & each network %
            netTemplate = cifti_template_mat_full(net,:);
            andImage = netTemplate & varmean_thresh;
            orImage = netTemplate | varmean_thresh;
            var_dice_to_template(net) = 2*sum(andImage)/(sum(orImage)+ sum(andImage));
        end
        
        dice_summary = [dice_summary; var_dice_to_template];
        
        %pick highest eta value, assign variant to that network, write the dice coefficient %%%
        [max_dice, winner] = max(var_dice_to_template',[],1);
        
        %% NEW - assign network as unassigned if the max dice is low
        
        if max_dice > min_match_dice
            winner = wb_colors(winner);
        else
            winner = 4;
        end
        
        if length(intersect(vertex_idx, find(consensus.data == winner))) > 0.5*length(vertex_idx)        %% If variant overlaps at least 50% with same canonical network
            winner = 99;
            max_dice = 0;
            
        else
            
            maxmatch = [maxmatch max_dice];
            
        end
        
        disp(['winner is ' num2str(winner)]);
        disp(['max dice is ' num2str(max_dice)]);
        
        varIDsreassign.data(vertex_idx) = winner;
        diceOfWinner.data(vertex_idx)=max_dice;
        
        %winners_summary(var, 1) = unique_vars(var);
        %winners_summary(var, 2) = winner;
        %winners_summary(var, 3)=max_dice;
        winners_summary=[winners_summary; unique_vars(var) winner max_dice];

    end
end

% write out a summary of the dice for everybody
dice_summary_path = [outfilepath outfilename '_dice_summary'];
dlmwrite(dice_summary_path, dice_summary);

% write out a summary of the winners 
winners_summary = winners_summary';
winners_summary_path = [outfilepath outfilename '_winners_summary'];
dlmwrite(winners_summary_path, winners_summary);

%write out ciftis- uniqueIDs, reassign %
network_assignment_path = [outfilepath outfilename '_assigned_networks.dtseries.nii']; 
ft_write_cifti_mod(network_assignment_path, varIDsreassign);

disp('done assigning networks');

end 

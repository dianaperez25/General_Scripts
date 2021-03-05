%% load stuff

% load difference maps
diffmat_file= '/Users/dianaperez/Desktop/Research/Lateralization_Variants/Permutation_Tests/VariantsvsFlippedVariants_LH_diffMap_1000_thresh_7.5_permutations_2.mat';
load(diffmat_file)
min_subject_threshold= .03;
diffmat_bin= logical(abs(diffmats)>= min_subject_threshold);

% load cifti vertex neighbors
neigh = load('/Users/dianaperez/Box/Quest_Backup/Scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_surf_neighbors_LR_normalwall.mat');
neigh = neigh.neighbors;

all_permutations_clusters= zeros(29696,size(diffmats,2));

%% label all clusters of N>=(min_subject_threshold) subjects by diff map
for c=1:size(all_permutations_clusters,2)
    %disp(num2str(c))
    
    %clusterized = zeros(59412,1);
    clusterized = zeros(29696,1);
    count=1;
    for i=1:29696

        % If the vertex is a potential variant, label it with a unique ID
        if diffmat_bin(i,c)

            % If the vertex is already part of a variant, continue using that ID
            if clusterized(i)>0
                id = clusterized(i);
                
            % Otherwise, use a new ID
            else
                id = count;
                count=count+1;
            end

            % Label it
            clusterized(i)=id;

            % Check that vertex's neighbors
            neighVerts = neigh(i,2:7);
            neighVerts(isnan(neighVerts))=[]; % Ignore NaNs

            % All nonzero neighbors are given the same ID
            for j=1:length(neighVerts)
                % Second part of if statement prevents overwriting previously assigned variants
                if clusterized(neighVerts(j))>0
                    clusterized(clusterized==clusterized(neighVerts(j)))=id;
                elseif diffmat_bin(neighVerts(j),c) && clusterized(neighVerts(j))==0
                    clusterized(neighVerts(j))=id;
                end
            end
        end
    end
    
    % make the IDs sequential
    ids = unique(clusterized); ids(1) = []; % Get rid of 0
    for i=1:length(ids)
        clusterized(logical(clusterized==ids(i)))=i;
    end
    
      
    all_permutations_clusters(:,c)=clusterized;
    clear clusterized
    
end


% get all cluster sizes (in # vertices)
all_cluster_sizes = [];
for i = 1:size(all_permutations_clusters,2)
    numClusters = length(unique(all_permutations_clusters(:,i))) - 1;
    for c = 1:numClusters
        all_cluster_sizes = [all_cluster_sizes; sum(all_permutations_clusters(:,i)==c)]; 
    end
end

all_cluster_sizes = sort(all_cluster_sizes,'descend');
top5p = round(length(all_cluster_sizes)*0.05);

cluster_threshold = all_cluster_sizes(top5p)

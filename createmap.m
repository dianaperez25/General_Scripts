function [overlap_map] = createmap(files, template)
allsubs = [];
groupmap = [];
nfiles = length(files);
disp(['Calculating overlap for ' num2str(nfiles) ' subjects in LH group...'])

for x = 1:nfiles
    % load subject data
    cifti = ft_read_cifti_mod(files{x});

    for q = 1:length(cifti.data)
        if cifti.data(q) > 0
            allsubs(q,x) = 1;
        elseif cifti.data(q) == 0
            allsubs(q,x) = 0;
        end
    end
    clear cifti
end

groupsum = sum(allsubs,2);
disp(['The maximum number of subjects that overlap in this group is ' num2str(max(groupsum))]) 
groupmap = sum(allsubs,2)/nfiles;
disp(['The maximum proportion of subjects that overlap in this group is ' num2str(max(groupmap))]) 
overlap_map = template;
overlap_map.data = groupmap;
end
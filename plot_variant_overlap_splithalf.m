%%overlap of variants split half analysis

outfile = '/Users/dianaperez/Desktop/sub-LS05_split-half_variant-overlap.dtseries.nii';
varmap_first = ft_read_cifti_mod('/Users/dianaperez/Box/Research/Lifespan/LS05_variants_first-half_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii');
varmap_second = ft_read_cifti_mod('/Users/dianaperez/Box/Research/Lifespan/LS05_variants_second-half_sizeExcluded_thresh-5_smooth_2.55.dtseries.nii');
data_first = varmap_first.data;
data_second = varmap_second.data;
new_cifti = zeros(size(data_first));
for v = 1:59412
    if data_first(v,1) == 1 && data_second(v,1) == 1
        new_cifti(v,1) = 3;
    elseif data_first(v,1) == 1 && data_second(v,1) == 0
        new_cifti(v,1) = 1;
    elseif data_first(v,1) == 0 && data_second(v,1) == 1
        new_cifti(v,1) = 2;
    elseif data_first(v,1) == 0 && data_second(v,1) == 0
        new_cifti(v,1) = 0;
    else
        disp('Check Variant IDs');
    end
end

varmap_first.data = new_cifti;
ft_write_cifti_mod(outfile, varmap_first);
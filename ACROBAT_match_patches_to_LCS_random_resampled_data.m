%% calculate percent match for randomly resampled ACROBAT data
% based on methods from "ACROBAT_match_patches_to_LCS.m"

% load in resampled and matched random data
load '/home/jmv208/ACROBAT/random_resample_all.mat';

%% how many patches have LCS

for i = 1:max(random_resample_all.patchID)
    
    ind_patch = find(random_resample_all.patchID ==i);
    random_resample_all.len_patch(i) = length(ind_patch);
    
    rpd_patch = random_resample_all.rpd_match(ind_patch);
    ind_rpd = find(rpd_patch > 0);
    random_resample_all.rpd_patch_positive(i) = length(ind_rpd);
    random_resample_all.rpd_patch_avg(i) = nanmean(rpd_patch);
    random_resample_all.rpd_patch_per_pos(i) = (length(ind_rpd)/length(ind_patch))*100;
    
    ftle_patch = random_resample_all.ftle_match(ind_patch);
    ind_ftle = find(ftle_patch > 0.3);
    random_resample_all.ftle_patch_positive(i) = length(ind_ftle);
    random_resample_all.ftle_patch_avg(i) = nanmean(ftle_patch);
    random_resample_all.ftle_patch_per_pos(i) = (length(ind_ftle)/length(ind_patch))*100;
    
    fprintf('patch complete\n');
end

random_resample_all_percentage = random_resample_all;

save('random_resample_all_percentage' , 'random_resample_all_percentage');

exit
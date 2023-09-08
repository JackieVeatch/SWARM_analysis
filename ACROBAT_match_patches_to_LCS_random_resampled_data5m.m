%% calculate percent match for randomly resampled ACROBAT data
% based on methods from "ACROBAT_match_patches_to_LCS.m"

% load in resampled and matched random data
load '/home/jmv208/ACROBAT/random_resample_all5m.mat';
load '/home/jmv208/ACROBAT/patchID_Adelie5m_matched.mat';

%% how many patches have LCS

for i = 1:max(random_resample_all5m.patchID)
    
    ind_patch = find(random_resample_all5m.patchID ==i);
    random_resample_all5m.len_patch(i) = length(ind_patch); % there are 100 randomized versions of this patch
    
    rpd_patch = random_resample_all5m.rpd_match(ind_patch);
    ind_rpd = find(rpd_patch > 0);
    random_resample_all5m.rpd_patch_positive(i) = length(ind_rpd);
    random_resample_all5m.rpd_patch_avg(i) = nanmean(rpd_patch);
    random_resample_all5m.rpd_patch_per_pos(i) = (length(ind_rpd)/length(ind_patch))*100;
    
    ftle_patch = random_resample_all5m.ftle_match(ind_patch);
    ind_ftle = find(ftle_patch > 0.3);
    random_resample_all5m.ftle_patch_positive(i) = length(ind_ftle);
    random_resample_all5m.ftle_patch_avg(i) = nanmean(ftle_patch);
    random_resample_all5m.ftle_patch_per_pos(i) = (length(ind_ftle)/length(ind_patch))*100;
    
    fprintf('patch complete\n');
end

random_resample_all_percentage5m = random_resample_all5m;

save('random_resample_all_percentage5m' , 'random_resample_all_percentage5m');


%% calculate len_patch variable for observed data

for i = 1:max(patchID_Adelie5m.patchID)
    
    ind_patch = find(patchID_Adelie5m.patchID ==i);
    patchID_Adelie5m.len_patch(i) = length(ind_patch);
    
    rpd_patch = patchID_Adelie5m.rpd(ind_patch);
    ind_rpd = find(rpd_patch > 0);
    patchID_Adelie5m.rpd_patch_positive(i) = length(ind_rpd);
    patchID_Adelie5m.rpd_patch_avg(i) = nanmean(rpd_patch);
    patchID_Adelie5m.rpd_patch_per_pos(i) = (length(ind_rpd)/length(ind_patch))*100;
    
    ftle_patch = patchID_Adelie5m.ftle(ind_patch);
    ind_ftle = find(ftle_patch > 0.3);
    patchID_Adelie5m.ftle_patch_positive(i) = length(ind_ftle);
    patchID_Adelie5m.ftle_patch_avg(i) = nanmean(ftle_patch);
    patchID_Adelie5m.ftle_patch_per_pos(i) = (length(ind_ftle)/length(ind_patch))*100;
    
    fprintf('patch complete\n');
end

random_resample_all_percentage5m = random_resample_all5m;

save('patchID_Adelie5m_matched' , 'patchID_Adelie5m');


exit
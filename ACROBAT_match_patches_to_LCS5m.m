%% Calculate how many phyto patches contain LCS
%jveatch 8June2023

% calculate patch length and assign patchID to patch definintions, adapted
% from 'ACROBAT_match_patches_to_LCS.m' --> this is the same analysis on the
% patch definitions integrating to only 5m instead of MLD

% calculate % match of RPD & FTLE to ACROBAT profiles (patch/no patch) with 'ACROBAT_match_variables.m'
% see "ACROBAT_match_nonpatches_to_LCS.m" for data on when LCS exist and
% patches do not, test the 'null hypothesis'

% Adelie transect definition of patches integrated to 5m
load '/Volumes/T7_shield/jmv208/scales/patchID_Adelie5m';

%% CODE FOR THE ADELIE TRANSECT DEFINITION OF PATCHES INTEGRATED TO 5M %%

% replace NaNs in patchID data surrounded by 1s with 1, and surrounded by
% zeroes with 0:
    % Example: 0 0 NaN NaN 0 0 becomes 0 0 0 0 0 0
    %           1 1 NaN NaN 1 1 becomes 1 1 1 1 1
    %           0 0 NaN NaN 1 1 stays 0 0 NaN NaN 1 1
patchID_Adelie5m.patch_noNaN = nan(size(patchID_Adelie5m.patch));

for i = 2:length(patchID_Adelie5m.patch)-2
    if isnan(patchID_Adelie5m.patch(i))

        validIndices_forward = find(~isnan(patchID_Adelie5m.patch(i:end)));
        validIndices_forward = validIndices_forward + (i-1);
        validIndices_backward = find(~isnan(patchID_Adelie5m.patch(1:i)));

        if patchID_Adelie5m.patch(validIndices_forward(1)) == patchID_Adelie5m.patch(validIndices_backward(end))
            patchID_Adelie5m.patch_noNaN(i) = patchID_Adelie5m.patch(validIndices_forward(1));
        else
            patchID_Adelie5m.patch_noNaN(i) = patchID_Adelie5m.patch(i);
        end

    else
        patchID_Adelie5m.patch_noNaN(i) = patchID_Adelie5m.patch(i);
    end
end

%%
patchID_Adelie5m.patchID = nan(size(patchID_Adelie5m.patch));
previous = 1;
id = [];
ID_BLOOM = 0;

for i = 2: length(patchID_Adelie5m.patch_noNaN) % loop through profs for one survey
    current = i;

    if patchID_Adelie5m.patch_noNaN(previous) == false && patchID_Adelie5m.patch_noNaN(current) == false
        % NOBLOOM
        IND_START = NaN;
        IND_END = NaN;
    end

    if patchID_Adelie5m.patch_noNaN(previous) == false && patchID_Adelie5m.patch_noNaN(current) == true
        % BLOOM START
        IND_START = current;
    end
    
    if isnan(patchID_Adelie5m.patch_noNaN(previous)) && patchID_Adelie5m.patch_noNaN(current) == true
        % BLOOM START after NaN value
        IND_START = current;
    end

    if patchID_Adelie5m.patch_noNaN(previous) == true && patchID_Adelie5m.patch_noNaN(current) == false && ~isnan(IND_START)
        % BLOOM END
        IND_END = previous;

        ID_BLOOM = ID_BLOOM + 1;
        id = [id ID_BLOOM];
        if IND_START > IND_END
            BLOOM_DURATION = IND_END;
        else
            BLOOM_DURATION = IND_END - IND_START;
        end
        patchID_Adelie5m.patchID(IND_START:IND_END) = ID_BLOOM;
        duration = [duration BLOOM_DURATION]; % bloom duration in units of # profiles
        
        % reset
        IND_START = NaN;
        IND_END = NaN;
        BLOOM_DURATION = NaN;

    end
    
    if patchID_Adelie5m.patch_noNaN(previous) == true && isnan(patchID_Adelie5m.patch_noNaN(current)) && ~isnan(IND_START)
        % BLOOM END
        IND_END = previous;

        ID_BLOOM = ID_BLOOM + 1;
        id = [id ID_BLOOM];
        if IND_START > IND_END
            BLOOM_DURATION = IND_END;
        else
            BLOOM_DURATION = IND_END - IND_START;
        end
        patchID_Adelie5m.patchID(IND_START:IND_END) = ID_BLOOM;
        duration = [duration BLOOM_DURATION]; % bloom duration in units of # profiles
        
        % reset
        IND_START = NaN;
        IND_END = NaN;
        BLOOM_DURATION = NaN;

    end

    if patchID_Adelie5m.patch_noNaN(previous) == true && patchID_Adelie5m.patch_noNaN(current) == true
        % BLOOM
    end

    if current == length(patchID_Adelie5m.patch_noNaN) && patchID_Adelie5m.patch_noNaN(current) == true 
        % day ends on a bloom
        BLOOM_DURATION = current - IND_START;
        duration = [duration BLOOM_DURATION];
        ID_BLOOM = ID_BLOOM +1;
        patchID_Adelie5m.patchID(IND_START:IND_END) = ID_BLOOM;
    end

    previous = current;

end % END OF THE LOOP


%% BELOW IS SOME PRELIMINARY ANALYSIS --> 
% better and more conclusive anaylsis done in 'compare_real_to_random.m'


%% how many patches have LCS

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
    
end

%% categorize patches into three size categories
% small patches = below the median
% median patches = median: one std above median
% large patches = larger than one std above median

med = median(patchID_Adelie5m.len_patch);
st_dev = std(patchID_Adelie5m.len_patch);
med_std = med+st_dev;

ind_small_patch = find(patchID_Adelie5m.len_patch <= med);
ind_med_patch = find(patchID_Adelie5m.len_patch > med & patchID_Adelie5m.len_patch <= med_std);
ind_large_patch = find(patchID_Adelie5m.len_patch > med_std);

small_patch_per_pos = patchID_Adelie5m.rpd_patch_per_pos(ind_small_patch);
med_patch_per_pos = patchID_Adelie5m.rpd_patch_per_pos(ind_med_patch);
large_patch_per_pos = patchID_Adelie5m.rpd_patch_per_pos(ind_large_patch);

%% number of patches that have any positive RPD
small_patch_pos = nnz(patchID_Adelie5m.rpd_patch_positive(ind_small_patch));
med_patch_pos = nnz(patchID_Adelie5m.rpd_patch_positive(ind_med_patch));
large_patch_pos = nnz(patchID_Adelie5m.rpd_patch_positive(ind_large_patch));

small_patch_pos_ftle = nnz(patchID_Adelie5m.ftle_patch_positive(ind_small_patch));
med_patch_pos_ftle = nnz(patchID_Adelie5m.ftle_patch_positive(ind_med_patch));
large_patch_pos_ftle = nnz(patchID_Adelie5m.ftle_patch_positive(ind_large_patch));

%% percent of patches that have more than half of their profiles with LCS

small_patch_50match = length(find(patchID_Adelie5m.rpd_patch_per_pos(ind_small_patch) > 50))/length(ind_small_patch);
med_patch_50match = length(find(patchID_Adelie5m.rpd_patch_per_pos(ind_med_patch) > 50))/length(ind_med_patch);
large_patch_50match = length(find(patchID_Adelie5m.rpd_patch_per_pos(ind_large_patch) > 50))/length(ind_large_patch);

small_patch_50match_ftle = length(find(patchID_Adelie5m.ftle_patch_per_pos(ind_small_patch) > 50))/length(ind_small_patch);
med_patch_50match_ftle = length(find(patchID_Adelie5m.ftle_patch_per_pos(ind_med_patch) > 50))/length(ind_med_patch);
large_patch_50match_ftle = length(find(patchID_Adelie5m.ftle_patch_per_pos(ind_large_patch) > 50))/length(ind_large_patch);


%%
% THINGS LEFT TO TRY %
% percent pos of non-patch prof and non-patch sections
% break patches into days and see which days were the best
% percentages of edges of large patches that have FTLE


% Felipe's code to determine if homogenous or not from RU32
    % returns 'Q' index of how mixed it is
% or you could do deep MLD vs shallow MLD from the ACROBAT data
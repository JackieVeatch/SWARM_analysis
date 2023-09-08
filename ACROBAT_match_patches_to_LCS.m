%% Calculate how many phyto patches contain LCS
% calculate % match of RPD & FTLE to ACROBAT profiles (patch/no patch) with 'ACROBAT_match_variables.m'
% see "ACROBAT_match_nonpatches_to_LCS.m" for data on when LCS exist and
% patches do not, test the 'null hypothesis'
% see "ACROBAT_match_patches_to_LCS_5m.m" for the same analysis on the
% patch definitions integrating to only 5m instead of MLD

% line by line definition of patches
load '/Volumes/home/jmv208/scales/patchID_edges';
% Adelie transect definition of patches
load '/Volumes/home/jmv208/scales/patchID_Adelie';

%% CODE FOR THE LINE BY LINE DEFINITION OF PATCHES %%

% replace NaNs in patchID data surrounded by 1s with 1, and surrounded by
% zeroes with 0:
    % Example: 0 0 NaN NaN 0 0 becomes 0 0 0 0 0 0
    %           1 1 NaN NaN 1 1 becomes 1 1 1 1 1
    %           0 0 NaN NaN 1 1 stays 0 0 NaN NaN 1 1
patchID.patch_noNaN = nan(size(patchID.patch));

for i = 2:length(patchID.patch)-1
    if isnan(patchID.patch(i))

        validIndices_forward = find(~isnan(patchID.patch(i:end)));
        validIndices_forward = validIndices_forward + (i-1);
        validIndices_backward = find(~isnan(patchID.patch(1:i)));

        if patchID.patch(validIndices_forward(1)) == patchID.patch(validIndices_backward(end))
            patchID.patch_noNaN(i) = patchID.patch(validIndices_forward(1));
        else
            patchID.patch_noNaN(i) = patchID.patch(i);
        end

    else
        patchID.patch_noNaN(i) = patchID.patch(i);
    end
end

%%
patchID.patchID = nan(size(patchID.patch));
previous = 1;
id = [];
ID_BLOOM = 0;

for i = 2: length(patchID.patch_noNaN) % loop through profs for one survey
    current = i;

    if patchID.patch_noNaN(previous) == false && patchID.patch_noNaN(current) == false
        % NOBLOOM
        IND_START = NaN;
        IND_END = NaN;
    end

    if patchID.patch_noNaN(previous) == false && patchID.patch_noNaN(current) == true
        % BLOOM START
        IND_START = current;
    end
    
    if isnan(patchID.patch_noNaN(previous)) && patchID.patch_noNaN(current) == true
        % BLOOM START
        IND_START = current;
    end

    if patchID.patch_noNaN(previous) == true && patchID.patch_noNaN(current) == false && ~isnan(IND_START)
        % BLOOM END
        IND_END = previous;

        ID_BLOOM = ID_BLOOM + 1;
        id = [id ID_BLOOM];
        if IND_START > IND_END
            BLOOM_DURATION = IND_END;
        else
            BLOOM_DURATION = IND_END - IND_START;
        end
        patchID.patchID(IND_START:IND_END) = ID_BLOOM;
        duration = [duration BLOOM_DURATION]; % bloom duration in units of # profiles
        
        % reset
        IND_START = NaN;
        IND_END = NaN;
        BLOOM_DURATION = NaN;

    end
    
    if patchID.patch_noNaN(previous) == true && isnan(patchID.patch_noNaN(current)) && ~isnan(IND_START)
        % BLOOM END
        IND_END = previous;

        ID_BLOOM = ID_BLOOM + 1;
        id = [id ID_BLOOM];
        if IND_START > IND_END
            BLOOM_DURATION = IND_END;
        else
            BLOOM_DURATION = IND_END - IND_START;
        end
        patchID.patchID(IND_START:IND_END) = ID_BLOOM;
        duration = [duration BLOOM_DURATION]; % bloom duration in units of # profiles
        
        % reset
        IND_START = NaN;
        IND_END = NaN;
        BLOOM_DURATION = NaN;

    end

    if patchID.patch_noNaN(previous) == true && patchID.patch_noNaN(current) == true
        % BLOOM
    end

    if current == length(patchID.patch_noNaN) && patchID.patch_noNaN(current) == true 
        % day ends on a bloom
        BLOOM_DURATION = current - IND_START;
        duration = [duration BLOOM_DURATION];
        ID_BLOOM = ID_BLOOM +1;
        patchID.patchID(IND_START:IND_END) = ID_BLOOM;
    end

    previous = current;

end % END OF THE LOOP

%% how many patches have LCS

for i = 1:max(patchID.patchID)
    
    ind_patch = find(patchID.patchID ==i);
    patchID.len_patch(i) = length(ind_patch);
    
    rpd_patch = patchID.rpd(ind_patch);
    ind_rpd = find(rpd_patch > 0);
    patchID.rpd_patch_positive(i) = length(ind_rpd);
    patchID.rpd_patch_avg(i) = nanmean(rpd_patch);
    patchID.rpd_patch_per_pos(i) = (length(ind_rpd)/length(ind_patch))*100;
    
    ftle_patch = patchID.ftle(ind_patch);
    ind_ftle = find(ftle_patch > 0.3);
    patchID.ftle_patch_positive(i) = length(ind_ftle);
    patchID.ftle_patch_avg(i) = nanmean(ftle_patch);
    patchID.ftle_patch_per_pos(i) = (length(ind_ftle)/length(ind_patch))*100;
    
end


%% CODE FOR THE ADELIE TRANSECT DEFINITION OF PATCHES %%

% replace NaNs in patchID data surrounded by 1s with 1, and surrounded by
% zeroes with 0:
    % Example: 0 0 NaN NaN 0 0 becomes 0 0 0 0 0 0
    %           1 1 NaN NaN 1 1 becomes 1 1 1 1 1
    %           0 0 NaN NaN 1 1 stays 0 0 NaN NaN 1 1
patchID_Adelie.patch_noNaN = nan(size(patchID_Adelie.patch));

for i = 2:length(patchID_Adelie.patch)-2
    if isnan(patchID_Adelie.patch(i))

        validIndices_forward = find(~isnan(patchID_Adelie.patch(i:end)));
        validIndices_forward = validIndices_forward + (i-1);
        validIndices_backward = find(~isnan(patchID_Adelie.patch(1:i)));

        if patchID_Adelie.patch(validIndices_forward(1)) == patchID_Adelie.patch(validIndices_backward(end))
            patchID_Adelie.patch_noNaN(i) = patchID_Adelie.patch(validIndices_forward(1));
        else
            patchID_Adelie.patch_noNaN(i) = patchID_Adelie.patch(i);
        end

    else
        patchID_Adelie.patch_noNaN(i) = patchID_Adelie.patch(i);
    end
end

%%
patchID_Adelie.patchID = nan(size(patchID_Adelie.patch));
previous = 1;
id = [];
ID_BLOOM = 0;

for i = 2: length(patchID_Adelie.patch_noNaN) % loop through profs for one survey
    current = i;

    if patchID_Adelie.patch_noNaN(previous) == false && patchID_Adelie.patch_noNaN(current) == false
        % NOBLOOM
        IND_START = NaN;
        IND_END = NaN;
    end

    if patchID_Adelie.patch_noNaN(previous) == false && patchID_Adelie.patch_noNaN(current) == true
        % BLOOM START
        IND_START = current;
    end
    
    if isnan(patchID_Adelie.patch_noNaN(previous)) && patchID_Adelie.patch_noNaN(current) == true
        % BLOOM START after NaN value
        IND_START = current;
    end

    if patchID_Adelie.patch_noNaN(previous) == true && patchID_Adelie.patch_noNaN(current) == false && ~isnan(IND_START)
        % BLOOM END
        IND_END = previous;

        ID_BLOOM = ID_BLOOM + 1;
        id = [id ID_BLOOM];
        if IND_START > IND_END
            BLOOM_DURATION = IND_END;
        else
            BLOOM_DURATION = IND_END - IND_START;
        end
        patchID_Adelie.patchID(IND_START:IND_END) = ID_BLOOM;
        duration = [duration BLOOM_DURATION]; % bloom duration in units of # profiles
        
        % reset
        IND_START = NaN;
        IND_END = NaN;
        BLOOM_DURATION = NaN;

    end
    
    if patchID_Adelie.patch_noNaN(previous) == true && isnan(patchID_Adelie.patch_noNaN(current)) && ~isnan(IND_START)
        % BLOOM END
        IND_END = previous;

        ID_BLOOM = ID_BLOOM + 1;
        id = [id ID_BLOOM];
        if IND_START > IND_END
            BLOOM_DURATION = IND_END;
        else
            BLOOM_DURATION = IND_END - IND_START;
        end
        patchID_Adelie.patchID(IND_START:IND_END) = ID_BLOOM;
        duration = [duration BLOOM_DURATION]; % bloom duration in units of # profiles
        
        % reset
        IND_START = NaN;
        IND_END = NaN;
        BLOOM_DURATION = NaN;

    end

    if patchID_Adelie.patch_noNaN(previous) == true && patchID_Adelie.patch_noNaN(current) == true
        % BLOOM
    end

    if current == length(patchID_Adelie.patch_noNaN) && patchID_Adelie.patch_noNaN(current) == true 
        % day ends on a bloom
        BLOOM_DURATION = current - IND_START;
        duration = [duration BLOOM_DURATION];
        ID_BLOOM = ID_BLOOM +1;
        patchID_Adelie.patchID(IND_START:IND_END) = ID_BLOOM;
    end

    previous = current;

end % END OF THE LOOP

%% how many patches have LCS

for i = 1:max(patchID_Adelie.patchID)
    
    ind_patch = find(patchID_Adelie.patchID ==i);
    patchID_Adelie.len_patch(i) = length(ind_patch);
    
    rpd_patch = patchID_Adelie.rpd(ind_patch);
    ind_rpd = find(rpd_patch > 0);
    patchID_Adelie.rpd_patch_positive(i) = length(ind_rpd);
    patchID_Adelie.rpd_patch_avg(i) = nanmean(rpd_patch);
    patchID_Adelie.rpd_patch_per_pos(i) = (length(ind_rpd)/length(ind_patch))*100;
    
    ftle_patch = patchID_Adelie.ftle(ind_patch);
    ind_ftle = find(ftle_patch > 0.3);
    patchID_Adelie.ftle_patch_positive(i) = length(ind_ftle);
    patchID_Adelie.ftle_patch_avg(i) = nanmean(ftle_patch);
    patchID_Adelie.ftle_patch_per_pos(i) = (length(ind_ftle)/length(ind_patch))*100;
    
end

%% categorize patches into three size categories
% small patches = below the median
% median patches = median: one std above median
% large patches = larger than one std above median

med = median(patchID_Adelie.len_patch);
st_dev = std(patchID_Adelie.len_patch);
med_std = med+st_dev;

ind_small_patch = find(patchID_Adelie.len_patch <= med);
ind_med_patch = find(patchID_Adelie.len_patch > med & patchID_Adelie.len_patch <= med_std);
ind_large_patch = find(patchID_Adelie.len_patch > med_std);

small_patch_per_pos = patchID_Adelie.rpd_patch_per_pos(ind_small_patch);
med_patch_per_pos = patchID_Adelie.rpd_patch_per_pos(ind_med_patch);
large_patch_per_pos = patchID_Adelie.rpd_patch_per_pos(ind_large_patch);

%% number of patches that have any positive RPD
small_patch_pos = nnz(patchID_Adelie.rpd_patch_positive(ind_small_patch));
med_patch_pos = nnz(patchID_Adelie.rpd_patch_positive(ind_med_patch));
large_patch_pos = nnz(patchID_Adelie.rpd_patch_positive(ind_large_patch));

small_patch_pos_ftle = nnz(patchID_Adelie.ftle_patch_positive(ind_small_patch));
med_patch_pos_ftle = nnz(patchID_Adelie.ftle_patch_positive(ind_med_patch));
large_patch_pos_ftle = nnz(patchID_Adelie.ftle_patch_positive(ind_large_patch));

%% percent of patches that have more than half of their profiles with LCS

small_patch_50match = length(find(patchID_Adelie.rpd_patch_per_pos(ind_small_patch) > 50))/length(ind_small_patch);
med_patch_50match = length(find(patchID_Adelie.rpd_patch_per_pos(ind_med_patch) > 50))/length(ind_med_patch);
large_patch_50match = length(find(patchID_Adelie.rpd_patch_per_pos(ind_large_patch) > 50))/length(ind_large_patch);

small_patch_50match_ftle = length(find(patchID_Adelie.ftle_patch_per_pos(ind_small_patch) > 50))/length(ind_small_patch);
med_patch_50match_ftle = length(find(patchID_Adelie.ftle_patch_per_pos(ind_med_patch) > 50))/length(ind_med_patch);
large_patch_50match_ftle = length(find(patchID_Adelie.ftle_patch_per_pos(ind_large_patch) > 50))/length(ind_large_patch);


%%
% THINGS LEFT TO TRY %
% percent pos of non-patch prof and non-patch sections
% break patches into days and see which days were the best
% percentages of edges of large patches that have FTLE


% Felipe's code to determine if homogenous or not from RU32
    % returns 'Q' index of how mixed it is
% or you could do deep MLD vs shallow MLD from the ACROBAT data
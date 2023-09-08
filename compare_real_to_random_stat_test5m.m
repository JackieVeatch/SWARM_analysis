%% How different are the Observed and Null LCS values for FTLE and RPD? 
% adapted from 'compare_real_to_random_stat_test.m' for 5m version of analysis

% first question: is the data normal?
    % no --> Wilcoxon Rank Sum test
        % Whereas the null hypothesis of the two-sample t test is equal means, 
        % the null hypothesis of the Wilcoxon test is usually taken as equal medians. 
        % Another way to think of the null is that the two populations have the same 
        % distribution with the same median. If we reject the null, that means we 
        % have evidence that one distribution is shifted to the left or right of the other. 
        % https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
    % yes --> ANOVA
        % is a statistical formula used to compare variances across the means 
        % (or average) of different groups. A range of scenarios use it to 
        % determine if there is any difference between the means of different groups.
        % https://www.tibco.com/reference-center/what-is-analysis-of-variance-anova#:~:text=Analysis%20of%20Variance%20(ANOVA)%20is,the%20means%20of%20different%20groups.
        
        
%% load in observed and null data
% created in last lines of 'compare_real_to_random.m'
load '/Volumes/T7_Shield/jmv208/ACROBAT/LCS_obs_null5m.mat'

%% is it normal? Shapiro-Wilk test for normality

addpath '/Volumes/T7_Shield/jmv208/scales/stats_functions/swtest'

fprintf('FTLE patch\n');

% Data to be tested
data = LCS_obs_null5M.rand_ftle_patch_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Null patch FTLE values are normally distributed.\n');
else
    fprintf('Hypothesis: Null patch FTLE values are not normally distributed.\n');
end

% Data to be tested
data = LCS_obs_null5M.obs_ftle_patch_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Observed patch FTLE values are normally distributed.\n');
else
    fprintf('Hypothesis: Observed patch FTLE values are not normally distributed.\n');
end

%% use Wilcoxon Rank Sum test because some data are not normal

% Data for two groups
group1 = LCS_obs_null5M.rand_ftle_patch_avg;
group2 = LCS_obs_null5M.obs_ftle_patch_avg;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2, 'tail', 'left');
% left ==> group 1 < group 1

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end

%% is it normal? Shapiro-Wilk test for normality

fprintf('************\n');
fprintf('FTLE edge\n');

% Data to be tested
data = LCS_obs_null5M.rand_ftle_edge_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Null edge FTLE values are normally distributed.\n');
else
    fprintf('Hypothesis: Null edge FTLE values are not normally distributed.\n');
end

% Data to be tested
data = LCS_obs_null5M.obs_ftle_edge_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Observed patch FTLE values are normally distributed.\n');
else
    fprintf('Hypothesis: Observed patch FTLE values are not normally distributed.\n');
end

%% use Wilcoxon Rank Sum test because some data are not normal

% Data for two groups
group1 = LCS_obs_null5M.rand_ftle_edge_avg;
group2 = LCS_obs_null5M.obs_ftle_edge_avg;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2, 'tail', 'left');
% left ==> group 1 < group 2

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end

%% is it normal? Shapiro-Wilk test for normality

fprintf('************\n');
fprintf('RPD patch\n');

% Data to be tested
data = LCS_obs_null5M.rand_rpd_patch_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Null patch RPD values are normally distributed.\n');
else
    fprintf('Hypothesis: Null patch RPD values are not normally distributed.\n');
end

% Data to be tested
data = LCS_obs_null5M.obs_rpd_patch_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Observed patch RPD values are normally distributed.\n');
else
    fprintf('Hypothesis: Observed patch RPD values are not normally distributed.\n');
end

%% use Wilcoxon Rank Sum test because some data are not normal

% Data for two groups
group1 = LCS_obs_null5M.rand_rpd_patch_avg;
group2 = LCS_obs_null5M.obs_rpd_patch_avg;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2, 'tail', 'left');

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end

%% is it normal? Shapiro-Wilk test for normality

fprintf('************\n');
fprintf('RPD edge\n');

% Data to be tested
data = LCS_obs_null5M.rand_rpd_edge_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Null edge RPD values are normally distributed.\n');
else
    fprintf('Hypothesis: Null edge RPD values are not normally distributed.\n');
end

% Data to be tested
data = LCS_obs_null5M.obs_rpd_edge_avg;

% Perform Shapiro-Wilk test
[h, p_value, ~] = swtest(data);

% Display results
fprintf('Shapiro-Wilk Test\n');
fprintf('-----------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: Observed edge RPD values are normally distributed.\n');
else
    fprintf('Hypothesis: Observed edge RPD values are not normally distributed.\n');
end

%% use Wilcoxon Rank Sum test because some data are not normal

% Data for two groups
group1 = LCS_obs_null5M.rand_rpd_edge_avg;
group2 = LCS_obs_null5M.obs_rpd_edge_avg;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2, 'tail', 'left');
%left ==> group 1 < group 2

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end

%% Same methods for data of obs-rand
fprintf('***STARTING OBS-RAND HISTOGRAMS***\n');

fprintf('Diff FTLE patch and edge\n');
% Data for two groups
group1 = LCS_obs_null5M.diff_ftle_patch;
group2 = LCS_obs_null5M.diff_ftle_edge;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2);

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end


%%
fprintf('Diff FTLE patch and edge strat days\n');
% Data for two groups
group1 = LCS_obs_null5M.diff_ftle_strat_edge;
group2 = LCS_obs_null5M.diff_ftle_strat_patch;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2);

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end



%%
fprintf('Diff FTLE patch and edge mixed days\n');
% Data for two groups
group1 = LCS_obs_null5M.diff_ftle_mixed_patch;
group2 = LCS_obs_null5M.diff_ftle_mixed_edge;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2);

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end

%%
fprintf('Diff FTLE edge mixed vs strat days\n');
% Data for two groups
group1 = LCS_obs_null5M.diff_ftle_mixed_edge;
group2 = LCS_obs_null5M.diff_ftle_strat_edge;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2);

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end

%%
fprintf('Diff FTLE patch mixed vs strat days\n');
% Data for two groups
group1 = LCS_obs_null5M.diff_ftle_mixed_patch;
group2 = LCS_obs_null5M.diff_ftle_strat_patch;

% Perform Wilcoxon rank sum test
[p_value, h] = ranksum(group1, group2);

% Display results
fprintf('Wilcoxon Rank Sum Test\n');
fprintf('----------------------\n');
fprintf('P-value: %.4f\n', p_value);
if h == 0
    fprintf('Hypothesis: The distributions are not significantly different.\n');
else
    fprintf('Hypothesis: The distributions are significantly different.\n');
end




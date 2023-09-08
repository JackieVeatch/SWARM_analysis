%% compare percent matches of random background phytoplankton patches and real patches

% created 14June2023 based off of "compare_real_to_random.m"
% using 5m integration definition of patch
% created to build context on LCS predicting phytoplankton patch locations

%%% DOUBLE CHECK THAT YOU ASSIGNED EDGE THE CORRECT WAY %%%

load '/Volumes/T7_Shield/jmv208/ACROBAT/random_resample_all_percentage5m.mat'
load '/Volumes/T7_Shield/jmv208/ACROBAT/patchID_Adelie5m_matched.mat'

% index for patch size, randomly resampled data
% NOTE: random resampled patches have 100 repeats, so length of 4 profiles = 400
r_ind_small_patch = find(random_resample_all_percentage5m.len_patch < 400);
r_ind_med_patch = find(random_resample_all_percentage5m.len_patch >= 400 & random_resample_all_percentage5m.len_patch <= 700);
r_ind_large_patch = find(random_resample_all_percentage5m.len_patch > 700);


% index for patch size, real data
a_ind_small_patch = find(patchID_Adelie5m.len_patch < 4);
a_ind_med_patch = find(patchID_Adelie5m.len_patch >= 4 & patchID_Adelie5m.len_patch <= 7);
a_ind_large_patch = find(patchID_Adelie5m.len_patch > 7);

patch_num = 1:length(random_resample_all_percentage5m.len_patch);

%% calculate std for each randomly resampled patch RPD

rpd_per_pos_all = [];
patch_num_all = [];

for i = 1:length(random_resample_all_percentage5m.len_patch)
    ind_patch = find(random_resample_all_percentage5m.patchID == i);
    rpd_match = random_resample_all_percentage5m.rpd_match(ind_patch);
    last_ind = 0;
    len = random_resample_all_percentage5m.len_patch(i)/100;
    for j = 1:100
        ind_patch = (last_ind+1):(last_ind + len);
        rpd_patch = rpd_match(ind_patch);
        rpd_per_pos(j) = (length(find(rpd_patch > 0))/len)*100;
        last_ind = ind_patch(end);
        
    end
  
    std_rpd_per_pos(i) = std(rpd_per_pos);
    rpd_per_pos_all = [rpd_per_pos_all, rpd_per_pos];
    patch_num_all = [patch_num_all, i.* ones(1,100)];
end

%% create plot to compare real with random based on patch size RPD

a_large_patch_rpd = patchID_Adelie5m.rpd_patch_per_pos(a_ind_large_patch);
a_med_patch_rpd = patchID_Adelie5m.rpd_patch_per_pos(a_ind_med_patch);
a_small_patch_rpd = patchID_Adelie5m.rpd_patch_per_pos(a_ind_small_patch);

r_large_patch_rpd = random_resample_all_percentage5m.rpd_patch_per_pos(r_ind_large_patch);
r_med_patch_rpd = random_resample_all_percentage5m.rpd_patch_per_pos(r_ind_med_patch);
r_small_patch_rpd = random_resample_all_percentage5m.rpd_patch_per_pos(r_ind_small_patch);

figure(1)
errorbar(patch_num(r_ind_large_patch), r_large_patch_rpd, std_rpd_per_pos(r_ind_large_patch));
hold on;
scatter(patch_num(a_ind_large_patch), a_large_patch_rpd, 'filled');
xlabel('patch number');
ylabel('percent of patch with positive RPD');
title('LARGE PATCHES: percent of patch with positive RPD, compare observed data to random');
ylim([-10, 120]);
legend('randomly resampled', 'observed');


figure(2)
errorbar(patch_num(r_ind_med_patch), r_med_patch_rpd, std_rpd_per_pos(r_ind_med_patch));
hold on;
scatter(patch_num(a_ind_med_patch), a_med_patch_rpd, 'filled');
xlabel('patch number');
ylabel('percent of patch with positive RPD');
title('MED PATCHES: percent of patch with positive RPD, compare observed data to random');
ylim([-10, 120]);
legend('randomly resampled', 'observed');


figure(3)
errorbar(patch_num(r_ind_small_patch), r_small_patch_rpd, std_rpd_per_pos(r_ind_small_patch));
hold on;
scatter(patch_num(a_ind_small_patch), a_small_patch_rpd, 'filled');
xlabel('patch number');
ylabel('percent of patch with positive RPD');
title('SMALL PATCHES: percent of patch with positive RPD, compare observed data to random');
ylim([-10, 120]);
legend('randomly resampled', 'observed');

%% large and medium patches together

figure(6)
errorbar(patch_num(r_ind_large_patch), r_large_patch_rpd, std_rpd_per_pos(r_ind_large_patch));
hold on;
scatter(patch_num(a_ind_large_patch), a_large_patch_rpd, 'filled');
errorbar(patch_num(r_ind_med_patch), r_med_patch_rpd, std_rpd_per_pos(r_ind_med_patch));
hold on;
scatter(patch_num(a_ind_med_patch), a_med_patch_rpd, 'filled');

xlabel('patch number');
ylabel('percent of patch with positive RPD');
title('LARGE + MEDIUM PATCHES: percent of patch with positive RPD, compare observed data to random');
ylim([-10, 120]);
legend('large randomly resampled', 'large observed', 'medium randomly resampled', 'medium observed');

%% calculate STD of average patch FTLE for randomly resampled data


for i = 1:length(random_resample_all_percentage5m.len_patch)
    ind_patch = find(random_resample_all_percentage5m.patchID == i);
    ftle_match = random_resample_all_percentage5m.ftle_match(ind_patch);
    last_ind = 0;
    len = random_resample_all_percentage5m.len_patch(i)/100;
    for j = 1:100
        ind_patch = (last_ind+1):(last_ind + len);
        ftle_patch = ftle_match(ind_patch);
        ftle_patch_avg(j) = mean(rpd_patch);
        last_ind = ind_patch(end);

    end
  
    std_ftle_patch_avg(i) = std(ftle_patch_avg);
end

%% subtract random background from obersved LCS average strength of patch

% assign timestamp to each patch
for i=1:length(patchID_Adelie5m.len_patch)
    
    ind_patch = find(patchID_Adelie5m.patchID == i);
    time = patchID_Adelie5m.timestamp(ind_patch);
    patchID_Adelie5m.patch_time(i) = mean(time);
    
end

% index patches into survey days
survey_days = ['15-Jan-2020'; '18-Jan-2020'; '21-Jan-2020'; '24-Jan-2020'; 
    '28-Jan-2020'; '01-Feb-2020';'05-Feb-2020';'07-Feb-2020';
    '12-Feb-2020'; '14-Feb-2020'; '18-Feb-2020'; '21-Feb-2020';'22-Feb-2020'; 
    '25-Feb-2020'; '28-Feb-2020'; '03-Mar-2020'];
% took out March 6th becuase there is no HFR data there

for i = 1:length(survey_days)
    
    day = survey_days(i,:);
    ind = find(patchID_Adelie5m.patch_time >= datenum(day) & patchID_Adelie5m.patch_time <= datenum(day)+1);
    indices_surveys{i} = ind;
    
end



%% Calculate difference between observed data and random background

diff_rpd = patchID_Adelie5m.rpd_patch_avg - random_resample_all_percentage5m.rpd_patch_avg;
diff_ftle = patchID_Adelie5m.ftle_patch_avg - random_resample_all_percentage5m.ftle_patch_avg;


%% calculate average difference between observed and random for each survey day


for i = 1:length(indices_surveys)
     rpd = diff_rpd(indices_surveys{i});
     diff_rpd_days(i) = mean(rpd);
     diff_rpd_days_std(i) = std(rpd);
     ftle = diff_ftle(indices_surveys{i});
     diff_ftle_days(i) = nanmean(ftle);
     diff_ftle_days_std(i) = nanstd(ftle);
end

%% calculate average difference between observed and random for just medium and large patches

a_ind_med_large = sort([a_ind_large_patch, a_ind_med_patch]);

diff_rpd_med_large = patchID_Adelie5m.rpd_patch_avg(a_ind_med_large) - random_resample_all_percentage5m.rpd_patch_avg(a_ind_med_large);
diff_ftle_med_large = patchID_Adelie5m.ftle_patch_avg(a_ind_med_large) - random_resample_all_percentage5m.ftle_patch_avg(a_ind_med_large);
timestamp_med_large = patchID_Adelie5m.patch_time(a_ind_med_large);

for i = 1:length(survey_days)
    day = survey_days(i,:);
    ind = find(timestamp_med_large >= datenum(day) & timestamp_med_large <= datenum(day)+1);
    diff_rpd_day = diff_rpd_med_large(ind);
    diff_ftle_day = diff_ftle_med_large(ind);
    diff_rpd_days_med_large(i) = nanmean(diff_rpd_day);
    diff_ftle_days_med_large(i) = nanmean(diff_ftle_day);
    
end

%% calculate difference between observed and random edges, (only M + L patches)

% index edges as the first profile of patch +/- 1 and the last profile of patch +/- 1
ind_edges_all = [];
edges = NaN(length(patchID_Adelie5m.patchID),1);
c = 1;

for i = 1:length(patchID_Adelie5m.len_patch)
    
    if patchID_Adelie5m.len_patch(i) >= 4
        ind_patch = find(patchID_Adelie5m.patchID == i);
        ind_edge1 = [(ind_patch(1) - 1), ind_patch(1), ind_patch(2)];
        ind_edge2 = [(ind_patch(end)-1), ind_patch(end), (ind_patch(end)+1)];
        edges(ind_edge1) = c;
        edges(ind_edge2) = c+1;
        c = c+2;
        ind_edges_all = [ind_edges_all, ind_edge1, ind_edge2];
    else
    end

end

patchID_Adelie5m.edgeID = edges;

%% calculate average ftle and rpd for each edge

for i = 1:max(patchID_Adelie5m.edgeID)
    
    ind_edge = find(patchID_Adelie5m.edgeID == i);
    rpd = nanmean(patchID_Adelie5m.rpd(ind_edge));
    ftle = nanmean(patchID_Adelie5m.ftle(ind_edge));
    patchID_Adelie5m.rpd_edge_avg(i) = rpd;
    patchID_Adelie5m.ftle_edge_avg(i) = ftle;
    
end

%take out first profile from each survey day, this is how you did the random resample u dummy

ind_edge_edited_all = [];

for i = 1:length(survey_days)
    day = survey_days(i,:);
    ind = find(patchID_Adelie5m.timestamp >= datenum(day) & patchID_Adelie5m.timestamp <= datenum(day)+1);
    ind_edge_survey = edges(ind);
    ind_edge_edited = ind_edge_survey(2:end)';
    ind_edge_edited_all = [ind_edge_edited_all, ind_edge_edited];
end
    
% repeat indices 100 times to match random resampled data
random_resample_all_percentage5m.edgeID = repmat(ind_edge_edited_all, 1,100);

% average LCS values for radomly resampled edges

for i = 1:max(random_resample_all_percentage5m.edgeID)
    
    ind_edge = find(random_resample_all_percentage5m.edgeID == i);
    rpd = nanmean(random_resample_all_percentage5m.rpd_match(ind_edge));
    ftle = nanmean(random_resample_all_percentage5m.ftle_match(ind_edge));
    random_resample_all_percentage5m.rpd_edge_avg(i) = rpd;
    random_resample_all_percentage5m.ftle_edge_avg(i) = ftle;
end

% calculate standard deviation for each randomly resampled edge
for i = 1:length(random_resample_all_percentage5m.rpd_edge_avg)
    ind_edge = find(random_resample_all_percentage5m.edgeID == i);
    ftle_match = random_resample_all_percentage5m.ftle_match(ind_edge);
    rpd_match = random_resample_all_percentage5m.ftle_match(ind_edge);
    last_ind = 0;
    len = length(ftle_match)/100;
    for j = 1:100
        ind_one_edge = (last_ind+1):(last_ind + len);
        ftle_edge = ftle_match(ind_one_edge);
        rpd_edge = rpd_match(ind_one_edge);
        ftle_edge_avg(j) = mean(ftle_edge);
        rpd_edge_avg(j) = mean(rpd_edge);
        last_ind = ind_one_edge(end);

    end
  
    std_ftle_edge_avg(i) = std(ftle_edge_avg);
    std_rpd_edge_avg(i) = std(rpd_edge_avg);
end

%% Calculate difference between randomly resampled edges and observed edges and plot

diff_edge_rpd = patchID_Adelie5m.rpd_edge_avg - random_resample_all_percentage5m.rpd_edge_avg;
diff_edge_ftle = patchID_Adelie5m.ftle_edge_avg - random_resample_all_percentage5m.ftle_edge_avg;

%assign timestamp to edges
for i = 1:max(patchID_Adelie5m.edgeID)
    ind = find(patchID_Adelie5m.edgeID == i);
    time = nanmean(patchID_Adelie5m.timestamp(ind));
    patchID_Adelie5m.edge_timestamp(i) = time;
end

%% read in ACROBAT data to look for stratified days

load('/Volumes/T7_Shield/jmv208/SWARM_data/acro_data_reprocessed/ACRO_reprocessed.mat');

% Small detour to calculate density difference
for i = 1:length(ACRO.mld)
    if isnan(ACRO.mld(i))
        ACRO.dens_diff(i) = NaN;
    else
        ind_no_nan = ~isnan(ACRO.dens(:,i));
        dens = ACRO.dens(:,i);
        dens_no_nan = dens(ind_no_nan);
        
        if dens_no_nan(end) < 20 % if the profile did not go below 20 meters
            ACRO.dens_diff(i) = NaN;
        end
        
        ACRO.dens_diff(i) = real(dens_no_nan(end) - dens_no_nan(2));
        
    end
end

% index ACRO data for just the ADELIE transect

ad_poly_lon = [-64.08, -64.21, -64.3, -64.1, -64.08];
ad_poly_lat = [-64.83, -64.79, -64.86, -64.89, -64.83];

adelie_all_inds = inpolygon(ACRO.lon, ACRO.lat, ad_poly_lon, ad_poly_lat);

%% find average density difference for each survey day
for i = 1:length(survey_days)
    
    day = survey_days(i,:);
    ind = find(ACRO.mtime >= datenum(day) & ACRO.mtime <= datenum(day)+1);
    dens_diff_day = ACRO.dens_diff(ind);
    ind_bad = find(dens_diff_day > 5); %some high values that I don't trust
    dens_diff_day(ind_bad) = NaN;
    dens_diff_daily_avg(i) = nanmean(dens_diff_day);
    dens_diff_daily_std(i) = nanstd(dens_diff_day);

end

%% stratified and unstratified days
t = datetime(survey_days);

median_dens_diff = median(dens_diff_daily_avg);
ind_strat = find(dens_diff_daily_avg >= median_dens_diff);
ind_mixed = find(dens_diff_daily_avg < median_dens_diff);

strat_days = t(ind_strat);
mixed_days = t(ind_mixed);

ftle_patch_ML = patchID_Adelie5m.ftle_patch_avg(a_ind_med_large);
ftle_rand_ML = random_resample_all_percentage5m.ftle_patch_avg(a_ind_med_large);
rpd_patch_ML = patchID_Adelie5m.rpd_patch_avg(a_ind_med_large);
rpd_rand_ML = random_resample_all_percentage5m.rpd_patch_avg(a_ind_med_large);

ftle_strat_patch = [];
ftle_strat_randpatch = [];
rpd_strat_patch = [];
rpd_strat_randpatch = [];
ftle_strat_edge = [];
ftle_strat_randedge = [];
rpd_strat_edge = [];
rpd_strat_randedge = [];

ftle_strat_patch_diff = [];
ftle_strat_edge_diff = [];


for i = 1:length(strat_days)
    
    day = datenum(strat_days(i));
    
    ind = find(timestamp_med_large > day & timestamp_med_large < day+1);
    ftle_strat_patch = [ftle_strat_patch, ftle_patch_ML(ind)];
    ftle_strat_randpatch = [ftle_strat_randpatch, ftle_rand_ML(ind)];
    rpd_strat_patch = [rpd_strat_patch, rpd_patch_ML(ind)];
    rpd_strat_randpatch = [rpd_strat_randpatch, rpd_rand_ML(ind)];
    ftle_strat_patch_diff = [ftle_strat_patch_diff, diff_ftle_med_large(ind)];
    
    ind = find(patchID_Adelie5m.edge_timestamp > day & patchID_Adelie5m.edge_timestamp < day+1);
    ftle_strat_edge = [ftle_strat_edge, patchID_Adelie5m.ftle_edge_avg(ind)];
    ftle_strat_randedge = [ftle_strat_randedge, random_resample_all_percentage5m.ftle_edge_avg(ind)];
    rpd_strat_edge = [rpd_strat_edge, patchID_Adelie5m.rpd_edge_avg(ind)];
    rpd_strat_randedge = [rpd_strat_randedge, random_resample_all_percentage5m.rpd_edge_avg(ind)];
    ftle_strat_edge_diff = [ftle_strat_edge_diff, diff_edge_ftle(ind)];
    
end

% and again for mixed days

ftle_mixed_patch = [];
ftle_mixed_randpatch = [];
rpd_mixed_patch = [];
rpd_mixed_randpatch = [];
ftle_mixed_edge = [];
ftle_mixed_randedge = [];
rpd_mixed_edge = [];
rpd_mixed_randedge = [];

ftle_mixed_patch_diff = [];
ftle_mixed_edge_diff = [];

for i = 1:length(mixed_days)
    
    day = datenum(mixed_days(i));
    
    ind = find(timestamp_med_large > day & timestamp_med_large < day+1);
    ftle_mixed_patch = [ftle_mixed_patch, ftle_patch_ML(ind)];
    ftle_mixed_randpatch = [ftle_mixed_randpatch, ftle_rand_ML(ind)];
    rpd_mixed_patch = [rpd_mixed_patch, rpd_patch_ML(ind)];
    rpd_mixed_randpatch = [rpd_mixed_randpatch, rpd_rand_ML(ind)];
    ftle_mixed_patch_diff = [ftle_mixed_patch_diff, diff_ftle_med_large(ind)];
    
    ind = find(patchID_Adelie5m.edge_timestamp > day & patchID_Adelie5m.edge_timestamp < day+1);
    ftle_mixed_edge = [ftle_mixed_edge, patchID_Adelie5m.ftle_edge_avg(ind)];
    ftle_mixed_randedge = [ftle_mixed_randedge, random_resample_all_percentage5m.ftle_edge_avg(ind)];
    rpd_mixed_edge = [rpd_mixed_edge, patchID_Adelie5m.rpd_edge_avg(ind)];
    rpd_mixed_randedge = [rpd_mixed_randedge, random_resample_all_percentage5m.rpd_edge_avg(ind)];
    ftle_mixed_edge_diff = [ftle_mixed_edge_diff, diff_edge_ftle(ind)];
    
end


%% Lets make some box and whisker plots already!!

figure(51)
A = ftle_strat_randpatch';
B = ftle_strat_patch';
C = ftle_mixed_randpatch';
D = ftle_mixed_patch';
group = [    ones(size(A));
         2 * ones(size(B))
         3 * ones(size(C))
         4 * ones(size(D))];
boxplot([A; B; C; D], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
ylim([0, 0.45]);
h = findobj(gca,'Tag','Box');
colors = [0 0.6 0.6; 0.8 0.4 0 ; 0 0.6 0.6; 0.8 0.4 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1.5 3.5]);
xticklabels({'stratified surveys', 'mixed surveys'});
legend('observed','randomly generated', 'Location', 'NorthWest');
title('FTLE values of patch centers on stratified and mixed surveys');

%%
figure(52)
A = random_resample_all_percentage5m.ftle_patch_avg(a_ind_med_large)';
B = patchID_Adelie5m.ftle_patch_avg(a_ind_med_large)';
C = random_resample_all_percentage5m.ftle_edge_avg';
D = patchID_Adelie5m.ftle_edge_avg';
group = [    ones(size(A));
         2 * ones(size(B))
         3 * ones(size(C))
         4 * ones(size(D))];
boxplot([A; B; C; D], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
ylim([0, 0.45]);
h = findobj(gca,'Tag','Box');
colors = [0 0.6 0.6; 0.8 0.4 0 ; 0 0.6 0.6; 0.8 0.4 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1.5 3.5]);
xticklabels({'patch centers', 'patch edge'});
legend('observed','randomly generated', 'Location', 'NorthWest');
title('FTLE values of patch centers and patch edges');
%%
figure(53)
A = random_resample_all_percentage5m.rpd_patch_avg(a_ind_med_large)';
B = patchID_Adelie5m.rpd_patch_avg(a_ind_med_large)';
C = random_resample_all_percentage5m.rpd_edge_avg';
D = patchID_Adelie5m.rpd_edge_avg';
group = [    ones(size(A));
         2 * ones(size(B))
         3 * ones(size(C))
         4 * ones(size(D))];
boxplot([A; B; C; D], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
ylim([-40, 100]);
h = findobj(gca,'Tag','Box');
colors = [0 0.6 0.6; 0.8 0.4 0 ; 0 0.6 0.6; 0.8 0.4 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1.5 3.5]);
xticklabels({'patch centers', 'patch edge'});
legend('observed','randomly generated', 'Location', 'NorthWest');
title('RPD values of patch centers and patch edges');

%% BW of difference between observed and random

figure(54)
A = diff_ftle_med_large';
B = diff_edge_ftle';
group = [ ones(size(A));
        2 * ones(size(B))];
boxplot([A; B], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
h = findobj(gca,'Tag','Box');
colors = [0 0.6 0; 0 0.6 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1 2]);
xticklabels({'patch centers', 'patch edge'});
title('FTLE difference between observed and random of patch centers and patch edges');

% again for RPD
figure(74)
A = diff_rpd_med_large';
B = diff_edge_rpd';
group = [ ones(size(A));
        2 * ones(size(B))];
boxplot([A; B], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
h = findobj(gca,'Tag','Box');
colors = [0.8 0.33 0; 0.8 0.33 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1 2]);
xticklabels({'patch centers', 'patch edge'});
title('RPD difference between observed and random of patch centers and patch edges');

figure(55)
A = ftle_strat_patch_diff';
B = ftle_strat_edge_diff';
C = ftle_mixed_patch_diff';
D = ftle_mixed_edge_diff';
group = [ ones(size(A));
        2 * ones(size(B))
        3* ones(size(C))
        4 * ones(size(D))];
boxplot([A; B; C; D], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
h = findobj(gca,'Tag','Box');
colors = [ 0.5 0.5 0.5; 0.5 0.5 0.5; 0.8 0 0; 0.8 0 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1 2 3 4]);
xticklabels({'patch centers', 'patch edge','patch centers', 'patch edge'});
legend('stratified surveys', 'mixed surveys', 'Location', 'SouthEast');
title('FTLE difference between observed and random of patch centers and patch edges, stratified surveys');

%% centers and edges for mixed days
figure(56)
A = ftle_strat_patch_diff';
B = ftle_mixed_patch_diff';
group = [ ones(size(A));
        2 * ones(size(B))];
boxplot([A; B], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
h = findobj(gca,'Tag','Box');
colors = [ 0.5 0.5 0.5; 0.8 0 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1 2]);
xticklabels({'stratified', 'mixed'});
hold on;
yline(0);
legend('mixed surveys','stratified surveys', 'zero line', 'Location', 'SouthEast');
title('FTLE difference between observed and random of patch centers');

%%
figure(57)
A = ftle_strat_edge_diff';
B = ftle_mixed_edge_diff';
group = [ ones(size(A));
        2 * ones(size(B))];
boxplot([A; B], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
h = findobj(gca,'Tag','Box');
colors = [ 0.5 0.5 0.5; 0.8 0 0];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1 2]);
xticklabels({'stratified', 'mixed'});
hold on;
yline(0);
legend('mixed surveys','stratified surveys', 'zero line', 'Location', 'SouthEast');
title('FTLE difference between observed and random of patch edges');
%%
figure(75)
A = diff_ftle_med_large';
B = diff_edge_ftle';
C = ftle_strat_patch_diff';
D = ftle_strat_edge_diff';
E = ftle_mixed_patch_diff';
F = ftle_mixed_edge_diff';
group = [ ones(size(A));
        2 * ones(size(B))
        3* ones(size(C))
        4 * ones(size(D))
        5 * ones(size(E))
        6 * ones(size(F))];
boxplot([A; B; C; D; E; F], group, 'notch', 'on', 'boxstyle', 'outline', 'symbol', 'k+', 'medianstyle', 'line');
h = findobj(gca,'Tag','Box');
colors = [ 0.5 0.5 0.5; 0.5 0.5 0.5; 0.8 0 0; 0.8 0 0 ; 0 0.6 0; 0 0.6 0 ];
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
xticks([1 2 3 4 5 6]);
xticklabels({'patch centers', 'patch edges', 'patch centers', 'patch edges', 'patch centers', 'patch edges'});
title('Difference Between FTLE Values of Observed and Null Model 5m integration');
yline(0)
legend('all surveys','stratified surveys', 'mixed surveys', 'Location', 'eastoutside');



%%
LCS_obs_null5M.rand_ftle_patch_avg = random_resample_all_percentage5m.ftle_patch_avg(a_ind_med_large)';
LCS_obs_null5M.obs_ftle_patch_avg = patchID_Adelie5m.ftle_patch_avg(a_ind_med_large)';
LCS_obs_null5M.rand_ftle_edge_avg = random_resample_all_percentage5m.ftle_edge_avg';
LCS_obs_null5M.obs_ftle_edge_avg = patchID_Adelie5m.ftle_edge_avg';
LCS_obs_null5M.rand_rpd_patch_avg = random_resample_all_percentage5m.rpd_patch_avg(a_ind_med_large)';
LCS_obs_null5M.obs_rpd_patch_avg = patchID_Adelie5m.rpd_patch_avg(a_ind_med_large)';
LCS_obs_null5M.rand_rpd_edge_avg = random_resample_all_percentage5m.rpd_edge_avg';
LCS_obs_null5M.obs_rpd_edge_avg = patchID_Adelie5m.rpd_edge_avg';
LCS_obs_null5M.diff_ftle_patch = diff_ftle_med_large';
LCS_obs_null5M.diff_ftle_edge = diff_edge_ftle';
LCS_obs_null5M.diff_ftle_strat_patch = ftle_strat_patch_diff';
LCS_obs_null5M.diff_ftle_strat_edge = ftle_strat_edge_diff';
LCS_obs_null5M.diff_ftle_mixed_patch = ftle_mixed_patch_diff';
LCS_obs_null5M.diff_ftle_mixed_edge = ftle_mixed_edge_diff';


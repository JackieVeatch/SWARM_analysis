%% compare percent matches of random background phytoplankton patches and real patches

% created 30Jan2023
% created to build context on LCS predicting phytoplankton patch locations

load '/Volumes/T7_Shield/jmv208/ACROBAT/random_resample_all_percentage.mat'
load '/Volumes/T7_Shield/jmv208/ACROBAT/patchID_Adelie_matched.mat'

% index for patch size, randomly resampled data
% NOTE: random resampled patches have 100 repeats, so length of 4 profiles = 400
r_ind_small_patch = find(random_resample_all_percentage.len_patch < 400);
r_ind_med_patch = find(random_resample_all_percentage.len_patch >= 400 & random_resample_all_percentage.len_patch <= 700);
r_ind_large_patch = find(random_resample_all_percentage.len_patch > 700);


% index for patch size, real data
a_ind_small_patch = find(patchID_Adelie.len_patch < 4);
a_ind_med_patch = find(patchID_Adelie.len_patch >= 4 & patchID_Adelie.len_patch <= 7);
a_ind_large_patch = find(patchID_Adelie.len_patch > 7);

patch_num = 1:length(random_resample_all_percentage.len_patch);

%% calculate std for each randomly resampled patch RPD

rpd_per_pos_all = [];
patch_num_all = [];

for i = 1:length(random_resample_all_percentage.len_patch)
    ind_patch = find(random_resample_all_percentage.patchID == i);
    rpd_match = random_resample_all_percentage.rpd_match(ind_patch);
    last_ind = 0;
    len = random_resample_all_percentage.len_patch(i)/100;
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

%% scatter all data
figure(7)
scatter(patch_num_all, rpd_per_pos_all);
hold on;
scatter(patch_num, patchID_Adelie.rpd_patch_per_pos, 'filled');
xlabel('patch number');
ylabel('percent of patch with positive RPD');
title('percent of patch with RPD, compare observed data to all random data');
legend('randomly resampled', 'observed');
ylim([-10, 120]);


%% create plot to compare real with random with error bars RPD
figure(4);
errorbar(patch_num, random_resample_all_percentage.rpd_patch_per_pos, std_rpd_per_pos);
hold on
scatter(patch_num, patchID_Adelie.rpd_patch_per_pos, 'filled');
xlabel('patch number');
ylabel('percent of patch with positive RPD');
title('percent of patch with RPD, compare observed data to random');
ylim([-10, 120]);
legend('randomly resapled', 'observed');


%% calculate std for each randomly resampled patch FTLE
for i = 1:length(random_resample_all_percentage.len_patch)
    ind_patch = find(random_resample_all_percentage.patchID == i);
    ftle_match = random_resample_all_percentage.ftle_match(ind_patch);
    last_ind = 0;
    len = random_resample_all_percentage.len_patch(i)/100;
    for j = 1:100
        ind_patch = (last_ind+1):(last_ind + len);
        ftle_patch = ftle_match(ind_patch);
        ftle_per_pos(j) = (length(find(ftle_patch > 0.3))/len)*100;
        last_ind = ind_patch(end);

    end
  
    std_ftle_per_pos(i) = std(ftle_per_pos);
end

patch_num = 1:length(random_resample_all_percentage.len_patch);

%% create plot to compare real with random FTLE
figure(5);
errorbar(patch_num, random_resample_all_percentage.ftle_patch_per_pos, std_ftle_per_pos);
hold on
scatter(patch_num, patchID_Adelie.ftle_patch_per_pos, 'filled');
xlabel('patch number');
ylabel('percent of patch with FTLE');
title('percent of patch with FTLE, compare observed data to random');
ylim([-10, 120]);
legend('randomly resampled', 'observed');

%% create plot to compare real with random based on patch size RPD

a_large_patch_rpd = patchID_Adelie.rpd_patch_per_pos(a_ind_large_patch);
a_med_patch_rpd = patchID_Adelie.rpd_patch_per_pos(a_ind_med_patch);
a_small_patch_rpd = patchID_Adelie.rpd_patch_per_pos(a_ind_small_patch);

r_large_patch_rpd = random_resample_all_percentage.rpd_patch_per_pos(r_ind_large_patch);
r_med_patch_rpd = random_resample_all_percentage.rpd_patch_per_pos(r_ind_med_patch);
r_small_patch_rpd = random_resample_all_percentage.rpd_patch_per_pos(r_ind_small_patch);

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

%% calulate percent of small/medium/large patches that are half RPD, comparing observed to random
% DOUBLE CHECK THIS IT DOESN'T QUITE MAKE SENSE

r_per_pos_rpd_small = random_resample_all_percentage.rpd_patch_per_pos(r_ind_small_patch);
r_per_pos_rpd_med = random_resample_all_percentage.rpd_patch_per_pos(r_ind_med_patch);
r_per_pos_rpd_large = random_resample_all_percentage.rpd_patch_per_pos(r_ind_large_patch);

r_per_half_rpd_small = (length(find(r_per_pos_rpd_small > 55))/length(r_ind_small_patch))*100;
r_per_half_rpd_med = (length(find(r_per_pos_rpd_med > 55))/length(r_ind_med_patch))*100;
r_per_half_rpd_large = (length(find(r_per_pos_rpd_large > 55))/length(r_ind_large_patch))*100;

a_per_pos_rpd_small = patchID_Adelie.rpd_patch_per_pos(a_ind_small_patch);
a_per_pos_rpd_med = patchID_Adelie.rpd_patch_per_pos(a_ind_med_patch);
a_per_pos_rpd_large = patchID_Adelie.rpd_patch_per_pos(a_ind_large_patch);

a_per_half_rpd_small = (length(find(a_per_pos_rpd_small > 55))/length(a_ind_small_patch))*100;
a_per_half_rpd_med = (length(find(a_per_pos_rpd_med > 55))/length(a_ind_med_patch))*100;
a_per_half_rpd_large = (length(find(a_per_pos_rpd_large > 55))/length(a_ind_large_patch))*100;

%% calculate median of randomly resampled data

for i = 1:length(patchID_Adelie.len_patch)
    ind_patch = find(patch_num_all == i);
    patch = rpd_per_pos_all(ind_patch);
    med = median(patch);
    rpd_per_pos_median(i) = med;
end


figure(8)
scatter(patch_num, rpd_per_pos_median);
hold on;
scatter(patch_num, patchID_Adelie.rpd_patch_per_pos);
hold on;
ylim([-10, 110]);
legend('random', 'observed');
ylim([-1, 120]);
title('percent of patch with positive RPD, compare observed data to median random');
xlabel('patch number');
ylabel('percent of patch with positive RPD');

%% compare random with observed average patch RPD

figure(9)
scatter(patch_num,random_resample_all_percentage.rpd_patch_avg);
hold on;
scatter(patch_num, patchID_Adelie.rpd_patch_avg, 'filled');
hold on;
legend('random', 'observed');
ylim([-1, 110]);
title('patch average RPD, compare observed data to mean random');
xlabel('patch number');
ylabel('percent of patch with positive RPD');

%% calculate STD of average patch RPD for randomly resampled data


for i = 1:length(random_resample_all_percentage.len_patch)
    ind_patch = find(random_resample_all_percentage.patchID == i);
    rpd_match = random_resample_all_percentage.rpd_match(ind_patch);
    last_ind = 0;
    len = random_resample_all_percentage.len_patch(i)/100;
    for j = 1:100
        ind_patch = (last_ind+1):(last_ind + len);
        rpd_patch = rpd_match(ind_patch);
        rpd_patch_avg(j) = mean(rpd_patch);
        last_ind = ind_patch(end);

    end
  
    std_rpd_patch_avg(i) = std(rpd_patch_avg);
end


figure(10)
errorbar(patch_num,random_resample_all_percentage.rpd_patch_avg, std_rpd_patch_avg);
hold on;
scatter(patch_num, patchID_Adelie.rpd_patch_avg, 'filled');
hold on;
legend('random', 'observed');
ylim([-1, 110]);
title('patch average RPD, compare observed data to mean random');
xlabel('patch number');
ylabel('average RPD value');

%%% compare random with observed average patch RPD

figure(9)
scatter(patch_num,random_resample_all_percentage.rpd_patch_avg);
hold on;
scatter(patch_num, patchID_Adelie.rpd_patch_avg, 'filled');
hold on;
legend('random', 'observed');
ylim([-1, 110]);
title('patch average RPD, compare observed data to mean random');
xlabel('patch number');
ylabel('percent of patch with positive RPD');

%% calculate STD of average patch FTLE for randomly resampled data


for i = 1:length(random_resample_all_percentage.len_patch)
    ind_patch = find(random_resample_all_percentage.patchID == i);
    ftle_match = random_resample_all_percentage.ftle_match(ind_patch);
    last_ind = 0;
    len = random_resample_all_percentage.len_patch(i)/100;
    for j = 1:100
        ind_patch = (last_ind+1):(last_ind + len);
        ftle_patch = ftle_match(ind_patch);
        ftle_patch_avg(j) = mean(rpd_patch);
        last_ind = ind_patch(end);

    end
  
    std_ftle_patch_avg(i) = std(ftle_patch_avg);
end


figure(10)
errorbar(patch_num,random_resample_all_percentage.ftle_patch_avg, std_ftle_patch_avg);
hold on;
scatter(patch_num, patchID_Adelie.ftle_patch_avg, 'filled');
hold on;
legend('random', 'observed');
ylim([-0.1, 0.6]);
title('patch average FTLE, compare observed data to mean random');
xlabel('patch number');
ylabel('average RPD value');


%% average LCS of patch compare by patch size

figure(11)
errorbar(patch_num,random_resample_all_percentage.rpd_patch_avg, std_rpd_patch_avg);
hold on;
scatter(patch_num(a_ind_small_patch), patchID_Adelie.rpd_patch_avg(a_ind_small_patch), 'filled', 'r');
scatter(patch_num(a_ind_med_patch), patchID_Adelie.rpd_patch_avg(a_ind_med_patch), 'filled', 'b');
scatter(patch_num(a_ind_large_patch), patchID_Adelie.rpd_patch_avg(a_ind_large_patch), 'filled', 'g');
hold on;
legend('random', 'observed small', 'observed medium', 'observed large');
ylim([-1, 110]);
title('patch average RPD, compare observed data to mean random');
xlabel('patch number');
ylabel('average RPD of patch');


figure(12)
errorbar(patch_num,random_resample_all_percentage.ftle_patch_avg, std_ftle_patch_avg);
hold on;
scatter(patch_num(r_ind_small_patch), patchID_Adelie.ftle_patch_avg(r_ind_small_patch), 'filled', 'r');
scatter(patch_num(r_ind_med_patch), patchID_Adelie.ftle_patch_avg(r_ind_med_patch), 'filled', 'b');
scatter(patch_num(r_ind_large_patch), patchID_Adelie.ftle_patch_avg(r_ind_large_patch), 'filled', 'g');
hold on;
legend('random', 'observed small', 'observed medium', 'observed large');
ylim([-0.1, 0.6]);
title('patch average FTLE, compare observed data to mean random');
xlabel('patch number');
ylabel('average FTLE of patch');

%% subtract random background from obersved LCS average strength of patch

% assign timestamp to each patch
for i=1:length(patchID_Adelie.len_patch)
    
    ind_patch = find(patchID_Adelie.patchID == i);
    time = patchID_Adelie.timestamp(ind_patch);
    patchID_Adelie.patch_time(i) = mean(time);
    
end

% index patches into survey days
survey_days = ['15-Jan-2020'; '18-Jan-2020'; '21-Jan-2020'; '24-Jan-2020'; 
    '28-Jan-2020'; '01-Feb-2020';'05-Feb-2020';'07-Feb-2020';
    '12-Feb-2020'; '14-Feb-2020'; '18-Feb-2020'; '21-Feb-2020';'22-Feb-2020'; 
    '25-Feb-2020'; '28-Feb-2020'; '03-Mar-2020'];
% took out March 6th becuase there is no HFR data there

for i = 1:length(survey_days)
    
    day = survey_days(i,:);
    ind = find(patchID_Adelie.patch_time >= datenum(day) & patchID_Adelie.patch_time <= datenum(day)+1);
    indices_surveys{i} = ind;
    
end



%% Calculate difference between observed data and random background

diff_rpd = patchID_Adelie.rpd_patch_avg - random_resample_all_percentage.rpd_patch_avg;
diff_ftle = patchID_Adelie.ftle_patch_avg - random_resample_all_percentage.ftle_patch_avg;


figure(11)
fill([indices_surveys{1}(1), indices_surveys{1}(1), indices_surveys{1}(end), indices_surveys{1}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
hold on;
fill([indices_surveys{3}(1), indices_surveys{3}(1), indices_surveys{3}(end), indices_surveys{3}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{5}(1), indices_surveys{5}(1), indices_surveys{5}(end), indices_surveys{5}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{7}(1), indices_surveys{7}(1), indices_surveys{7}(end), indices_surveys{7}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{9}(1), indices_surveys{9}(1), indices_surveys{9}(end), indices_surveys{9}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{11}(1), indices_surveys{11}(1), indices_surveys{11}(end), indices_surveys{11}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{13}(1), indices_surveys{13}(1), indices_surveys{13}(end), indices_surveys{13}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{15}(1), indices_surveys{15}(1), indices_surveys{15}(end), indices_surveys{15}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
scatter(patch_num, diff_rpd, 'filled');
yline(0);

figure(12)
fill([indices_surveys{1}(1), indices_surveys{1}(1), indices_surveys{1}(end), indices_surveys{1}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
hold on;
fill([indices_surveys{3}(1), indices_surveys{3}(1), indices_surveys{3}(end), indices_surveys{3}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{5}(1), indices_surveys{5}(1), indices_surveys{5}(end), indices_surveys{5}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{7}(1), indices_surveys{7}(1), indices_surveys{7}(end), indices_surveys{7}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{9}(1), indices_surveys{9}(1), indices_surveys{9}(end), indices_surveys{9}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{11}(1), indices_surveys{11}(1), indices_surveys{11}(end), indices_surveys{11}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{13}(1), indices_surveys{13}(1), indices_surveys{13}(end), indices_surveys{13}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{15}(1), indices_surveys{15}(1), indices_surveys{15}(end), indices_surveys{15}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
scatter(patch_num, diff_ftle, 'filled');
yline(0);

%% plot difference between observed data and random background by patch size

figure(13)
fill([indices_surveys{1}(1), indices_surveys{1}(1), indices_surveys{1}(end), indices_surveys{1}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
hold on;
fill([indices_surveys{3}(1), indices_surveys{3}(1), indices_surveys{3}(end), indices_surveys{3}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{5}(1), indices_surveys{5}(1), indices_surveys{5}(end), indices_surveys{5}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{7}(1), indices_surveys{7}(1), indices_surveys{7}(end), indices_surveys{7}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{9}(1), indices_surveys{9}(1), indices_surveys{9}(end), indices_surveys{9}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{11}(1), indices_surveys{11}(1), indices_surveys{11}(end), indices_surveys{11}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{13}(1), indices_surveys{13}(1), indices_surveys{13}(end), indices_surveys{13}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
fill([indices_surveys{15}(1), indices_surveys{15}(1), indices_surveys{15}(end), indices_surveys{15}(end)],[-0.3,0.4,0.4,-0.3], [0.75 0.75 0.75]);
scatter(patch_num(a_ind_small_patch), diff_ftle(a_ind_small_patch), 200, 'r.');
scatter(patch_num(a_ind_med_patch), diff_ftle(a_ind_med_patch), 200, 'g.');
scatter(patch_num(a_ind_large_patch), diff_ftle(a_ind_large_patch), 200, 'b.');
xlabel('patch number shaded by survey day');
ylabel('difference between average observed and average random FTLE');
title('Difference between Observed and Random average FTLE value');

yline(0);

figure(14)
fill([indices_surveys{1}(1), indices_surveys{1}(1), indices_surveys{1}(end), indices_surveys{1}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
hold on;
fill([indices_surveys{3}(1), indices_surveys{3}(1), indices_surveys{3}(end), indices_surveys{3}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{5}(1), indices_surveys{5}(1), indices_surveys{5}(end), indices_surveys{5}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{7}(1), indices_surveys{7}(1), indices_surveys{7}(end), indices_surveys{7}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{9}(1), indices_surveys{9}(1), indices_surveys{9}(end), indices_surveys{9}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{11}(1), indices_surveys{11}(1), indices_surveys{11}(end), indices_surveys{11}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{13}(1), indices_surveys{13}(1), indices_surveys{13}(end), indices_surveys{13}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
fill([indices_surveys{15}(1), indices_surveys{15}(1), indices_surveys{15}(end), indices_surveys{15}(end)],[-50,450,450,-50], [0.75 0.75 0.75]);
scatter(patch_num(a_ind_small_patch), diff_rpd(a_ind_small_patch), 200, 'r.');
scatter(patch_num(a_ind_med_patch), diff_rpd(a_ind_med_patch), 200, 'g.');
scatter(patch_num(a_ind_large_patch), diff_rpd(a_ind_large_patch), 200, 'b.');
ylim([-50, 150]);
yline(0);
xlabel('patch number shaded by survey day');
ylabel('difference between average observed and average random RPD');
title('Difference between Observed and Random average RPD value');

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

diff_rpd_med_large = patchID_Adelie.rpd_patch_avg(a_ind_med_large) - random_resample_all_percentage.rpd_patch_avg(a_ind_med_large);
diff_ftle_med_large = patchID_Adelie.ftle_patch_avg(a_ind_med_large) - random_resample_all_percentage.ftle_patch_avg(a_ind_med_large);
timestamp_med_large = patchID_Adelie.patch_time(a_ind_med_large);

for i = 1:length(survey_days)
    day = survey_days(i,:);
    ind = find(timestamp_med_large >= datenum(day) & timestamp_med_large <= datenum(day)+1);
    diff_rpd_day = diff_rpd_med_large(ind);
    diff_ftle_day = diff_ftle_med_large(ind);
    diff_rpd_days_med_large(i) = nanmean(diff_rpd_day);
    diff_ftle_days_med_large(i) = nanmean(diff_ftle_day);
    
end

% plot the difference between observed and random for RPD and FTLE of med + large patches
figure(15)
scatter(datetime(timestamp_med_large, 'ConvertFrom', 'datenum'), diff_ftle_med_large, 300, 'k.');
hold on;
yline(0);
xlabel('Survey Day');
ylabel('difference between observed and random');
title('Average FTLE of patch difference between observed and random, medium and large patches');

figure(16)
scatter(datetime(timestamp_med_large, 'ConvertFrom', 'datenum'), diff_rpd_med_large, 300, 'k.');
hold on;
yline(0);
xlabel('Survey Day');
ylabel('difference between observed and random');
title('Average RPD of patch difference between observed and random, medium and large patches');

%% calculate difference between observed and random edges, (only M + L patches)

% index edges as the first profile of patch +/- 1 and the last profile of patch +/- 1
ind_edges_all = [];
edges = NaN(length(patchID_Adelie.patchID),1);
c = 1;

for i = 1:length(patchID_Adelie.len_patch)
    
    if patchID_Adelie.len_patch(i) >= 4
        ind_patch = find(patchID_Adelie.patchID == i);
        ind_edge1 = [(ind_patch(1) - 1), ind_patch(1), ind_patch(2)];
        ind_edge2 = [(ind_patch(end)-1), ind_patch(end), (ind_patch(end)+1)];
        edges(ind_edge1) = c;
        edges(ind_edge2) = c+1;
        c = c+2;
        ind_edges_all = [ind_edges_all, ind_edge1, ind_edge2];
    else
    end

end

patchID_Adelie.edgeID = edges;

%% calculate average ftle and rpd for each edge

for i = 1:max(patchID_Adelie.edgeID)
    
    ind_edge = find(patchID_Adelie.edgeID == i);
    rpd = nanmean(patchID_Adelie.rpd(ind_edge));
    ftle = nanmean(patchID_Adelie.ftle(ind_edge));
    patchID_Adelie.rpd_edge_avg(i) = rpd;
    patchID_Adelie.ftle_edge_avg(i) = ftle;
    
end

%take out first profile from each survey day, this is how you did the random resample u dummy

ind_edge_edited_all = [];

for i = 1:length(survey_days)
    day = survey_days(i,:);
    ind = find(patchID_Adelie.timestamp >= datenum(day) & patchID_Adelie.timestamp <= datenum(day)+1);
    ind_edge_survey = edges(ind);
    ind_edge_edited = ind_edge_survey(2:end)';
    ind_edge_edited_all = [ind_edge_edited_all, ind_edge_edited];
end
    
% repeat indices 100 times to match random resampled data
random_resample_all_percentage.edgeID = repmat(ind_edge_edited_all, 1,100);

% average LCS values for radomly resampled edges

for i = 1:max(random_resample_all_percentage.edgeID)
    
    ind_edge = find(random_resample_all_percentage.edgeID == i);
    rpd = nanmean(random_resample_all_percentage.rpd_match(ind_edge));
    ftle = nanmean(random_resample_all_percentage.ftle_match(ind_edge));
    random_resample_all_percentage.rpd_edge_avg(i) = rpd;
    random_resample_all_percentage.ftle_edge_avg(i) = ftle;
end

% calculate standard deviation for each randomly resampled edge
for i = 1:length(random_resample_all_percentage.rpd_edge_avg)
    ind_edge = find(random_resample_all_percentage.edgeID == i);
    ftle_match = random_resample_all_percentage.ftle_match(ind_edge);
    rpd_match = random_resample_all_percentage.ftle_match(ind_edge);
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

diff_edge_rpd = patchID_Adelie.rpd_edge_avg - random_resample_all_percentage.rpd_edge_avg;
diff_edge_ftle = patchID_Adelie.ftle_edge_avg - random_resample_all_percentage.ftle_edge_avg;

%assign timestamp to edges
for i = 1:max(patchID_Adelie.edgeID)
    ind = find(patchID_Adelie.edgeID == i);
    time = nanmean(patchID_Adelie.timestamp(ind));
    patchID_Adelie.edge_timestamp(i) = time;
end

figure(17)
scatter(datetime(patchID_Adelie.edge_timestamp, 'ConvertFrom', 'datenum'), diff_edge_rpd, 300, 'k.');
hold on;
yline(0);
xlabel('Survey Day');
ylabel('difference between observed and random');
title('Average RPD of edge difference between observed and random, medium and large patches');

figure(18)
scatter(datetime(patchID_Adelie.edge_timestamp, 'ConvertFrom', 'datenum'), diff_edge_ftle, 300, 'k.');
hold on;
yline(0);
xlabel('Survey Day');
ylabel('difference between observed and random');
title('Average FTLE of edge difference between observed and random, medium and large patches');

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

%% master plot of FTLE info

fig = figure(19);
left_color = [.5 .5 .5];
right_color = [1 1 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

scatter(datetime(patchID_Adelie.edge_timestamp, 'ConvertFrom', 'datenum'), diff_edge_ftle, 300, 'k.');
hold on;
scatter(datetime(timestamp_med_large, 'ConvertFrom', 'datenum'), diff_ftle_med_large, 300, 'g.');
yline(0);
yyaxis left
ylabel('difference in FTLE value, observed minus random (1/day)');
ylim([-0.2, 0.8]);
t = datetime(survey_days);
yyaxis right
scatter(t, dens_diff_daily_avg);
ylabel('daily averaged density difference 50m minus surface (kg/m3)');
ylim([-0.2, 1.3]);
title('Difference between observed and random FTLE values compared to daily stratification');
legend('edge difference', 'patch difference', 'zero difference', 'stratification', 'Location', 'NorthWest');


%% read in wind data
load('/Volumes/T7_Shield/jmv208/SWARM_data/Winds/2020/Joubins_2020_all.mat');

Joub.matlab_time = matlab_time;
Joub.GMT = matlab_time + 3/24;
Joub.winddir = winddir_all;
Joub.windgust = windgust_all;
Joub.windgustdir = windgustdir_all;
Joub.windspeed = windspeed_all;
Joub.createdBy = 'Data cleaned by Josh, struct created by Jackie';
Joub.doc = "'BuildMetTimeSeries.m' by Josh and 'compare_real_random.m' by Jackie";
Joub.rawdatadir = '/Volumes/T7_Shield/jmv208/SWARM_data/SWARM_Pal_Weather';

load('/Volumes/T7_Shield/jmv208/SWARM_data/Winds/2020/Wauwermans_2020_all.mat');

Wauw.matlab_time = matlab_time;
Wauw.GMT = matlab_time + 3/24;
Wauw.winddir = winddir_all;
Wauw.windgust = windgust_all;
Wauw.windgustdir = windgustdir_all;
Wauw.windspeed = windspeed_all;
Wauw.createdBy = 'Data cleaned by Josh, struct created by Jackie';
Wauw.doc = "'BuildMetTimeSeries.m' by Josh and 'compare_real_random.m' by Jackie";
Wauw.rawdatadir = '/Volumes/T7_Shield/jmv208/SWARM_data/SWARM_Pal_Weather';

% average winds for the day before and the day of survey

for i = 1:length(survey_days)
    day = datenum(survey_days(i,:));
    ind1 = find(Joub.GMT > day-1 & Joub.GMT < day);
    ind2 = find(Joub.GMT > day & Joub.GMT < day+1);
    joub_daybefore_avg_windsp(i) = nanmean(Joub.windspeed(ind1));
    joub_daily_avg_windsp(i) = nanmean(Joub.windspeed(ind2));
    
    ind1 = find(Wauw.GMT > day-1 & Wauw.GMT < day);
    ind2 = find(Wauw.GMT > day & Wauw.GMT < day+1);
    wauw_daybefore_avg_windsp(i) = nanmean(Wauw.windspeed(ind1));
    wauw_daily_avg_windsp(i) = nanmean(Wauw.windspeed(ind2));
    
end

%% master plot of RPD info with wind and stratification

figure(20)
tcl = tiledlayout(3,1);

nexttile
scatter(datetime(patchID_Adelie.edge_timestamp, 'ConvertFrom', 'datenum'), diff_edge_rpd, 300, 'k.');
hold on;
scatter(datetime(timestamp_med_large, 'ConvertFrom', 'datenum'), diff_rpd_med_large, 300, 'g.');
legend('edge difference', 'patch difference', 'zero difference');
yline(0);
ylabel('RPD difference');
ylim([-50, 150]);

nexttile
t = datetime(survey_days);
scatter(t, dens_diff_daily_avg, 300, 'b.');
ylabel('density difference (kg/m3)');

nexttile
scatter(t, joub_daily_avg_windsp, '+'); hold on;
scatter(t, joub_daybefore_avg_windsp, 'o');
ylabel('wind speed');
legend('day of survey', 'day before survey', 'Location', 'NorthWest');

title(tcl,'Relative Particle Density, compare observed to random with environment');


figure(21)
tcl = tiledlayout(3,1);

nexttile
scatter(datetime(patchID_Adelie.edge_timestamp, 'ConvertFrom', 'datenum'), diff_edge_ftle, 300, 'k.');
hold on;
scatter(datetime(timestamp_med_large, 'ConvertFrom', 'datenum'), diff_ftle_med_large, 300, 'g.');
legend('edge difference', 'patch difference', 'zero difference');
yline(0);
ylabel('FTLE difference');
ylim([-0.2, 0.5]);

nexttile
t = datetime(survey_days);
scatter(t, dens_diff_daily_avg, 300, 'b.');
ylabel('density difference (kg/m3)');

nexttile
scatter(t, joub_daily_avg_windsp, '+'); hold on;
scatter(t, joub_daybefore_avg_windsp, 'o');
ylabel('wind speed');
legend('day of survey', 'day before survey', 'Location', 'NorthWest');

title(tcl,'FTLE, compare observed to random with environment');

%% box and whisker plots of random and observed 'populations'

% figure(21)
% boxplot(patchID_Adelie.ftle_patch_avg(a_ind_med_large),'Notch','on', 'Labels', 'FTLE observed patch');
% hold on; ylim([0, 0.5]);
% 
% figure(22)
% boxplot(random_resample_all_percentage.ftle_patch_avg(a_ind_med_large), 'Notch','on','Labels', 'FTLE random patch');
% hold on; ylim([0, 0.5]);
% 
% figure(23)
% boxplot(patchID_Adelie.rpd_patch_avg(a_ind_med_large), 'Notch','on','Labels', 'RPD observed patch');
% hold on; ylim([-25, 80]);
% 
% figure(24)
% boxplot(random_resample_all_percentage.rpd_patch_avg(a_ind_med_large),'Notch','on', 'Labels', 'RPD random patch');
% hold on; ylim([-25, 80]);
% 
% figure(25)
% boxplot(patchID_Adelie.ftle_edge_avg, 'Notch','on','Labels', 'FTLE observed edge');
% hold on; ylim([0, 0.5]);
% 
% figure(26)
% boxplot(random_resample_all_percentage.ftle_edge_avg, 'Notch','on','Labels', 'FTLE random edge');
% hold on; ylim([0, 0.5]);
% 
% figure(27)
% boxplot(patchID_Adelie.rpd_edge_avg,'Notch','on', 'Labels', 'RPD observed edge');
% hold on; ylim([-25, 80]);
% 
% figure(28)
% boxplot(random_resample_all_percentage.rpd_edge_avg,'Notch','on', 'Labels', 'RPD random edge');
% hold on; ylim([-25, 80]);

%% box and whisker plots for stratified and unstratified days

median_dens_diff = median(dens_diff_daily_avg);
ind_strat = find(dens_diff_daily_avg >= median_dens_diff);
ind_mixed = find(dens_diff_daily_avg < median_dens_diff);

strat_days = t(ind_strat);
mixed_days = t(ind_mixed);

ftle_patch_ML = patchID_Adelie.ftle_patch_avg(a_ind_med_large);
ftle_rand_ML = random_resample_all_percentage.ftle_patch_avg(a_ind_med_large);
rpd_patch_ML = patchID_Adelie.rpd_patch_avg(a_ind_med_large);
rpd_rand_ML = random_resample_all_percentage.rpd_patch_avg(a_ind_med_large);

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
    
    ind = find(patchID_Adelie.edge_timestamp > day & patchID_Adelie.edge_timestamp < day+1);
    ftle_strat_edge = [ftle_strat_edge, patchID_Adelie.ftle_edge_avg(ind)];
    ftle_strat_randedge = [ftle_strat_randedge, random_resample_all_percentage.ftle_edge_avg(ind)];
    rpd_strat_edge = [rpd_strat_edge, patchID_Adelie.rpd_edge_avg(ind)];
    rpd_strat_randedge = [rpd_strat_randedge, random_resample_all_percentage.rpd_edge_avg(ind)];
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
    
    ind = find(patchID_Adelie.edge_timestamp > day & patchID_Adelie.edge_timestamp < day+1);
    ftle_mixed_edge = [ftle_mixed_edge, patchID_Adelie.ftle_edge_avg(ind)];
    ftle_mixed_randedge = [ftle_mixed_randedge, random_resample_all_percentage.ftle_edge_avg(ind)];
    rpd_mixed_edge = [rpd_mixed_edge, patchID_Adelie.rpd_edge_avg(ind)];
    rpd_mixed_randedge = [rpd_mixed_randedge, random_resample_all_percentage.rpd_edge_avg(ind)];
    ftle_mixed_edge_diff = [ftle_mixed_edge_diff, diff_edge_ftle(ind)];
    
end

%% plotty plot plot plot

figure(29)
boxplot(ftle_strat_patch,'Notch','on', 'Labels', 'FTLE observed patch stratified days');
hold on; ylim([0.1, 0.45]);

figure(30)
boxplot(ftle_strat_randpatch,'Notch','on', 'Labels', 'FTLE random patch stratified days');
hold on; ylim([0.1, 0.45]);


figure(31)
boxplot(ftle_strat_edge,'Notch','on', 'Labels', 'FTLE observed edge stratified days');
hold on; ylim([0.1, 0.45]);

figure(32)
boxplot(ftle_strat_randedge,'Notch','on', 'Labels', 'FTLE random edge stratified days');
hold on; ylim([0.1, 0.45]);

%% plotty plot mcplotface

figure(45)
boxplot(ftle_mixed_patch,'Notch','on', 'Labels', 'FTLE observed patch mixed days');
hold on; ylim([0.1, 0.45]);

figure(46)
boxplot(ftle_mixed_randpatch,'Notch','on', 'Labels', 'FTLE random patch mixed days');
hold on; ylim([0.1, 0.45]);


figure(47)
boxplot(ftle_mixed_edge,'Notch','on', 'Labels', 'FTLE observed edge mixed days');
hold on; ylim([0.1, 0.45]);

figure(48)
boxplot(ftle_mixed_randedge,'Notch','on', 'Labels', 'FTLE random edge mixed days');
hold on; ylim([0.1, 0.45]);

%% sum big plots
clear patch
figure(50)
A = ftle_strat_randedge';
B = ftle_strat_edge';
C = ftle_mixed_randedge';
D = ftle_mixed_edge';
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
title('FTLE values of patch edges on stratified and mixed surveys');
%%
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
A = random_resample_all_percentage.ftle_patch_avg(a_ind_med_large)';
B = patchID_Adelie.ftle_patch_avg(a_ind_med_large)';
C = random_resample_all_percentage.ftle_edge_avg';
D = patchID_Adelie.ftle_edge_avg';
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
A = random_resample_all_percentage.rpd_patch_avg(a_ind_med_large)';
B = patchID_Adelie.rpd_patch_avg(a_ind_med_large)';
C = random_resample_all_percentage.rpd_edge_avg';
D = patchID_Adelie.rpd_edge_avg';
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





%% BW of difference between observed and random, stratified and mixed days

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

%% combine 54 and 55 in some plot (oh man here we go)

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
title('Difference Between FTLE Values of Observed and Null Model');
yline(0)
legend('all surveys','stratified surveys', 'mixed surveys', 'Location', 'eastoutside');

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

%% does stratification correlate with how well FTLE predict patches

% assign average stratification to each patch

for i = 1:length(patchID_Adelie.len_patch)
    ind = find(patchID_Adelie.patchID == i);
    dens_diff = patchID_Adelie.dens_diff(ind);
    patchID_Adelie.patch_dens_diff(i) = nanmean(dens_diff);
end

for i = 1:length(patchID_Adelie.edge_timestamp)
    ind = find(patchID_Adelie.edgeID == i);
    dens_diff = patchID_Adelie.dens_diff(ind);
    patchID_Adelie.edge_dens_diff(i) = nanmean(dens_diff); 
end
    
    
patch_dens_diff_ML = patchID_Adelie.patch_dens_diff(a_ind_med_large);

figure(33)
scatter(patch_dens_diff_ML, diff_ftle_med_large);
hold on;
xlim([0 2]);
yline(0);
xline(median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));
ylabel('FTLE patch observed - random');
xlabel('density difference of patch');

figure(34)
scatter(patchID_Adelie.edge_dens_diff, diff_edge_ftle);
hold on;
xlim([0 2]);
yline(0);
xline(median(patchID_Adelie.edge_dens_diff));
ylabel('FTLE edge observed - random');
xlabel('density difference of edge');

%% scatter survey averages of density difference and LCS observed - random
figure(44)
errorbar(dens_diff_daily_avg, diff_ftle_days,diff_ftle_days_std./2, diff_ftle_days_std./2, dens_diff_daily_std./2, dens_diff_daily_std./2, 'LineStyle', 'none');
hold on;
c = NaN(length(dens_diff_daily_avg), 3);
for i = 1:length(dens_diff_daily_avg)
    c(i,:,:) = rand(1,3);
    scatter(dens_diff_daily_avg(i), diff_ftle_days(i),100, c(i,:,:), 'filled');
end   
xlabel('density difference survey average');
ylabel('FTLE observed - random survey average');
title('Comparison of Stratification and FTLE performance, Survey Averages');
yline(0);
xline(median(dens_diff_daily_avg));

figure(45)
errorbar(dens_diff_daily_avg, diff_rpd_days,diff_rpd_days_std./2, diff_rpd_days_std./2, dens_diff_daily_std./2, dens_diff_daily_std./2, 'LineStyle', 'none');
hold on;
for i = 1:length(dens_diff_daily_avg)
    scatter(dens_diff_daily_avg(i), diff_rpd_days(i),100, c(i,:,:), 'filled');
end   
xlabel('density difference survey average');
ylabel('RPD observed - random survey average');
title('Comparison of Stratification and RPD performance, Survey Averages');
ylim([ -20 , 50]);
yline(0);
xline(median(dens_diff_daily_avg));

%% calculate daily difference of edges

diff_ftle_days_edges = [];
diff_rpd_days_edges = [];
diff_ftle_days_edges_std = [];
diff_rpd_days_edges_std = [];

for i = 1:length(survey_days)
    
    day = datenum(survey_days(i,:));
    ind = find(patchID_Adelie.edge_timestamp > day & patchID_Adelie.edge_timestamp < day+1);
    
    diff_ftle_days_edges = [diff_ftle_days_edges, nanmean(diff_edge_ftle(ind))];
    diff_rpd_days_edges = [diff_rpd_days_edges, nanmean(diff_edge_rpd(ind))];
    
    diff_ftle_days_edges_std = [diff_ftle_days_edges_std, nanstd(diff_edge_ftle(ind))];
    diff_rpd_days_edges_std = [diff_rpd_days_edges_std, nanstd(diff_edge_rpd(ind))];
    
end

% errorbar(X,Y,YNEG,YPOS,XNEG,XPOS)

figure(46)
errorbar(dens_diff_daily_avg, diff_ftle_days_edges,diff_ftle_days_edges_std./2, diff_ftle_days_edges_std./2, dens_diff_daily_std./2, dens_diff_daily_std./2, 'LineStyle', 'none');
hold on;
for i = 1:length(dens_diff_daily_avg)
    scatter(dens_diff_daily_avg(i), diff_ftle_days_edges(i),100, c(i,:,:), 'filled');
end
xlabel('density difference survey average');
ylabel('FTLE edges observed - random survey average');
title('Comparison of Stratification and FTLE edge performance, Survey Averages');
yline(0);
xline(median(dens_diff_daily_avg));

figure(47)
errorbar(dens_diff_daily_avg, diff_rpd_days_edges,diff_rpd_days_edges_std./2, diff_rpd_days_edges_std./2, dens_diff_daily_std./2, dens_diff_daily_std./2, 'LineStyle', 'none');
hold on;
for i = 1:length(dens_diff_daily_avg)
    scatter(dens_diff_daily_avg(i), diff_rpd_days_edges(i), 100, c(i,:,:), 'filled');
end
xlabel('density difference survey average');
ylabel('RPD edges observed - random survey average');
title('Comparison of Stratification and RPD edge performance, Survey Averages');
yline(0);
xline(median(dens_diff_daily_avg));

%% Comparison of Stratification and FTLE performance, Survey Averages for four case studies

case_study_days = ['03-Mar-2020'; '28-Jan-2020'; '21-Feb-2020'; '24-Jan-2020'];
dnum = datenum(survey_days);
% color = ('b', 

for i = 1:4
    day = datenum(case_study_days(i,:));
    ind_day = find(dnum == day);
    s(i) = errorbar(dens_diff_daily_avg(ind_day), diff_ftle_days_edges(ind_day), diff_ftle_days_edges_std(ind_day)./2,diff_ftle_days_edges_std(ind_day)./2, dens_diff_daily_std(ind_day)./2, dens_diff_daily_std(ind_day)./2);
    hold on;
    scatter(dens_diff_daily_avg(ind_day), diff_ftle_days_edges(ind_day), 'ko');
end
legend( s, '03-Mar-2020', '28-Jan-2020', '21-Feb-2020', '24-Jan-2020', 'Location', 'SouthEast');
xlabel('density difference survey average');
ylabel('FTLE edges observed - random survey average');
title('Comparison of Stratification and FTLE edge performance, Survey Averages Case Study Days');
yline(0);
xline(median(dens_diff_daily_avg));

%%
for i = 1:4
    day = datenum(case_study_days(i,:));
    ind_day = find(dnum == day);
    s(i) = errorbar(dens_diff_daily_avg(ind_day), diff_ftle_days(ind_day), diff_ftle_days_std(ind_day)./2,diff_ftle_days_std(ind_day)./2, dens_diff_daily_std(ind_day)./2, dens_diff_daily_std(ind_day)./2);
    hold on;
    scatter(dens_diff_daily_avg(ind_day), diff_ftle_days(ind_day), 'ko');
end
legend( s, '03-Mar-2020', '28-Jan-2020', '21-Feb-2020', '24-Jan-2020', 'Location', 'SouthEast');
xlabel('density difference survey average');
ylabel('FTLE patches observed - random survey average');
title('Comparison of Stratification and FTLE patch performance, Survey Averages Case Study Days');
yline(0);
xline(median(dens_diff_daily_avg));

%% how many patches have high strat and high difference (observed-random LCS)

med_dens_diff = median(patchID_Adelie.patch_dens_diff(a_ind_med_large));

figure(35)
histogram(patchID_Adelie.patch_dens_diff(a_ind_med_large), 40); hold on;
title('density difference of medium and large patches');

ind_strat_lcs_patch = find(diff_ftle_med_large > 0 & patchID_Adelie.patch_dens_diff(a_ind_med_large) > median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));
ind_strat_nolcs_patch = find(diff_ftle_med_large <= 0 & patchID_Adelie.patch_dens_diff(a_ind_med_large) > median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));
ind_nostrat_lcs_patch = find(diff_ftle_med_large > 0 & patchID_Adelie.patch_dens_diff(a_ind_med_large) <= median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));
ind_nostrat_nolcs_patch = find(diff_ftle_med_large <= 0 & patchID_Adelie.patch_dens_diff(a_ind_med_large) <= median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));

%%

ind_strat_lcs_edge = find(diff_edge_ftle >= 0 & patchID_Adelie.edge_dens_diff >= median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));
ind_strat_nolcs_edge = find(diff_edge_ftle <= 0 & patchID_Adelie.edge_dens_diff >= median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));
ind_nostrat_lcs_edge = find(diff_edge_ftle >= 0 & patchID_Adelie.edge_dens_diff <= median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));
ind_nostrat_nolcs_edge = find(diff_edge_ftle <= 0 & patchID_Adelie.edge_dens_diff <= median(patchID_Adelie.patch_dens_diff(a_ind_med_large)));


%% box and whisker plots for high wind and low wind days

daily_wind_speed = mean(horzcat(joub_daily_avg_windsp', wauw_daily_avg_windsp'), 2);

median_wind_sp = median(daily_wind_speed);
ind_windy = find(daily_wind_speed >= median_wind_sp);
ind_calm = find(daily_wind_speed < median_wind_sp);

windy_days = t(ind_windy);
calm_days = t(ind_calm);

ftle_patch_ML = patchID_Adelie.ftle_patch_avg(a_ind_med_large);
ftle_rand_ML = random_resample_all_percentage.ftle_patch_avg(a_ind_med_large);
rpd_patch_ML = patchID_Adelie.rpd_patch_avg(a_ind_med_large);
rpd_rand_ML = random_resample_all_percentage.rpd_patch_avg(a_ind_med_large);

ftle_windy_patch = [];
ftle_windy_randpatch = [];
rpd_windy_patch = [];
rpd_windy_randpatch = [];
ftle_windy_edge = [];
ftle_windy_randedge = [];
rpd_windy_edge = [];
rpd_windy_randedge = [];


ftle_calm_patch = [];
ftle_calm_randpatch = [];
rpd_calm_patch = [];
rpd_calm_randpatch = [];
ftle_calm_edge = [];
ftle_calm_randedge = [];
rpd_calm_edge = [];
rpd_calm_randedge = [];


for i = 1:length(windy_days)
    
    day = datenum(windy_days(i));
    
    ind = find(timestamp_med_large > day & timestamp_med_large < day+1);
    ftle_windy_patch = [ftle_windy_patch, ftle_patch_ML(ind)];
    ftle_windy_randpatch = [ftle_windy_randpatch, ftle_rand_ML(ind)];
    rpd_windy_patch = [rpd_windy_patch, rpd_patch_ML(ind)];
    rpd_windy_randpatch = [rpd_windy_randpatch, rpd_rand_ML(ind)];
    
    ind = find(patchID_Adelie.edge_timestamp > day & patchID_Adelie.edge_timestamp < day+1);
    ftle_windy_edge = [ftle_windy_edge, patchID_Adelie.ftle_edge_avg(ind)];
    ftle_windy_randedge = [ftle_windy_randedge, random_resample_all_percentage.ftle_edge_avg(ind)];
    rpd_windy_edge = [rpd_windy_edge, patchID_Adelie.rpd_edge_avg(ind)];
    rpd_windy_randedge = [rpd_windy_randedge, random_resample_all_percentage.rpd_edge_avg(ind)];
    
end

for i = 1:length(calm_days)
    
    day = datenum(calm_days(i));
    
    ind = find(timestamp_med_large > day & timestamp_med_large < day+1);
    ftle_calm_patch = [ftle_calm_patch, ftle_patch_ML(ind)];
    ftle_calm_randpatch = [ftle_calm_randpatch, ftle_rand_ML(ind)];
    rpd_calm_patch = [rpd_calm_patch, rpd_patch_ML(ind)];
    rpd_calm_randpatch = [rpd_calm_randpatch, rpd_rand_ML(ind)];
    
    ind = find(patchID_Adelie.edge_timestamp > day & patchID_Adelie.edge_timestamp < day+1);
    ftle_calm_edge = [ftle_calm_edge, patchID_Adelie.ftle_edge_avg(ind)];
    ftle_calm_randedge = [ftle_calm_randedge, random_resample_all_percentage.ftle_edge_avg(ind)];
    rpd_calm_edge = [rpd_calm_edge, patchID_Adelie.rpd_edge_avg(ind)];
    rpd_calm_randedge = [rpd_calm_randedge, random_resample_all_percentage.rpd_edge_avg(ind)];
    
end


%% plotty plot plot plot

figure(36)
boxplot(ftle_windy_patch,'Notch','on', 'Labels', 'FTLE observed patch windy days');
hold on; ylim([0.1, 0.45]);

figure(37)
boxplot(ftle_windy_randpatch,'Notch','on', 'Labels', 'FTLE random patch windy days');
hold on; ylim([0.1, 0.45]);


figure(38)
boxplot(ftle_windy_edge,'Notch','on', 'Labels', 'FTLE observed edge windy days');
hold on; ylim([0.1, 0.45]);

figure(39)
boxplot(ftle_windy_randedge,'Notch','on', 'Labels', 'FTLE random edge windy days');
hold on; ylim([0.1, 0.45]);

%% plotty plot plot plot

figure(40)
boxplot(ftle_calm_patch,'Notch','on', 'Labels', 'FTLE observed patch calm days');
hold on; ylim([0.1, 0.45]);

figure(41)
boxplot(ftle_calm_randpatch,'Notch','on', 'Labels', 'FTLE random patch calm days');
hold on; ylim([0.1, 0.45]);


figure(42)
boxplot(ftle_calm_edge,'Notch','on', 'Labels', 'FTLE observed edge calm days');
hold on; ylim([0.1, 0.45]);

figure(43)
boxplot(ftle_calm_randedge,'Notch','on', 'Labels', 'FTLE random edge calm days');
hold on; ylim([0.1, 0.45]);

%% example figure of patch center vs edges on a survey day

% index for just the large and medium patches
ind_ml = find(patchID_Adelie.len_patch >= 4);
patchID_ml = nan(size(patchID_Adelie.patchID));
for i = 1:length(ind_ml)
    ind_patch = find(patchID_Adelie.patchID == ind_ml(i));
    patchID_ml(ind_patch) = ind_ml(i);
end
patch_binary_ml = nan(size(patchID_Adelie.patch));
patch_binary_ml(~isnan(patchID_ml)) = 1;

% index for one example day
ind_plot = find(patchID_Adelie.timestamp >= datenum(survey_days(9, :)) & patchID_Adelie.timestamp < datenum(survey_days(9,:))+1);
lon = patchID_Adelie.lon(ind_plot);
lat = patchID_Adelie.lat(ind_plot);
patch_binary = patch_binary_ml(ind_plot);
edge_binary = patchID_Adelie.edgeID(ind_plot);
lon_patch = lon(patch_binary == 1);
lat_patch = lat(patch_binary == 1);
lon_edge = lon(~isnan(edge_binary));
lat_edge = lat(~isnan(edge_binary));

addpath(genpath('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting'));
%%
figure(44)

 % plot antarctic bathymetry    
	bathy=load ('antarctic_bathy_2.mat');
	ind2= bathy.depthi==99999;
	bathy.depthi(ind2)=[];
	bathylines1=0:-10:-100;
	bathylines2=0:-200:-1400;
	bathylines=[bathylines2];
	
	[cs, h1] = contour(bathy.loni,bathy.lati, bathy.depthi,bathylines, 'linewidth', .25);
	clabel(cs,h1,'fontsize',6);
	set(h1,'LineColor','black')
	hold on

s(1) = scatter(lon, lat, 'ko');
hold on;
s(2) = scatter(lon_patch, lat_patch, 25, [0.4660 0.6740 0.1880], 'filled');
% s(3) = scatter(lon_edge, lat_edge, 25, [0 0.4470 0.7410],'filled');
s(3) = scatter(lon_edge, lat_edge, 25, [0.8 0.33 0],'filled');

    % the shape file used for this code is split into three different
    % segments, seperated by NaNs. Older versions of MATLAB do not have a
    % problem with this, newer versions (like MATLAB_R2019) do not like
    % this. The code below splits the shape file into three readable
    % pieces.
    tanLand = [240,230,140]./255;
    S1 = shaperead('cst00_polygon_wgs84.shp');
    S2=S1(1:1174);
    ind=[0,find(isnan(S1(1175).X))];
    for x=1:length(ind)-1
        S2(1174+x)=S1(1175);
        S2(1174+x).X=S2(1174+x).X(ind(x)+1:ind(x+1));
        S2(1174+x).Y=S2(1174+x).Y(ind(x)+1:ind(x+1));
    end
    mapshow(S2,'facecolor', tanLand)
    hold on
    
    %marking location of CODAR stations
    s(4) = plot(-64.0554167, -64.7741833, 'g^',...
        'markersize', 12,...
        'markerfacecolor', [0.5 0.5 0.5],...
        'markeredgecolor', 'black');
    s(5) = plot(-64.3604167, -64.7871667, 'gs',...
        'markersize', 12,...
        'markerfacecolor', [0.5 0.5 0.5],...
        'markeredgecolor', 'black');
    s(6) = plot(-64.0446333, -64.9183167, 'gd',...
        'markersize', 12,...
        'markerfacecolor', [0.5 0.5 0.5],...
        'markeredgecolor', 'black');
    
%     scatter(-64.1246, -64.8451, 4000,'MarkerFaceColor', [0.5 0.5 0.5]);  
    s(7) = scatter(-64.1246, -64.8451, 400,'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'black');
    
    ylim([-64.95, -64.75]);
    xlim([-64.45 -63.95]);
    title('Palmer Deep Canyon','FontSize', 16);
%     a=narrow(-63.6928384831366,-64.683125,.3); % place north facing arrow on upper right corner of map
    %l = legend([s,adgr,gegr], 'Palmer Station', 'Wauwermans Islands', 'Joubin Islands', 'Adelie Transect', 'Gentoo Transect','Location', 'SouthEast');
    

    
%     [lon,lat] = meshgrid(CODAR.lon, CODAR.lat);
%     poly_lon = [lon(49,60), lon(48,59), lon(49,58), lon(48,57), lon(49,56), lon(50,55), lon(50,54), lon(51,53),lon(51,52), lon( 51,51), lon(52,50),lon(53,49), lon(54,48), lon(54,47), lon(53,46), lon(53,45), lon(54,45), lon(55,44), lon(59,44), lon(59,45), lon(62,45), lon(64,47), lon(65,47),  lon(65,48), lon(66,49), lon(67,50), lon(68,50), lon(68,52), lon(69, 53), lon(69,55), lon(70,55), lon(70,57), lon(71,57), lon(71,59), lon(72,59), lon(73, 60), lon(73,62), lon(74,62), lon(75, 63), lon(75,66), lon(76,67), lon(76,71), lon(68,75), lon(63,78), lon(63,79), lon(61,79), lon(60,78), lon(56,78), lon(56,76), lon(55,76), lon(55,74), lon(56,73), lon(56,72), lon(54,72), lon(54,71), lon(53,71), lon(53,68), lon(52,67), lon(51,66), lon(51,64), lon(50,63),lon(50,61), lon(49,60)]; 
%     poly_lat = [lat(49,60), lat(48,59), lat(49,58), lat(48,57), lat(49,56), lat(50,55), lat(50,54), lat(51,53),lat(51,52), lat( 51,51), lat(52,50),lat(53,49), lat(54,48), lat(54,47), lat(53,46), lat(53,45), lat(54,45), lat(55,44), lat(59,44), lat(59,45), lat(62,45), lat(64,47), lat(65,47),  lat(65,48), lat(66,49), lat(67,50), lat(68,50), lat(68,52), lat(69, 53), lat(69,55), lat(70,55), lat(70,57), lat(71,57), lat(71,59), lat(72,59), lat(73, 60), lat(73,62), lat(74,62), lat(75, 63), lat(75,66), lat(76,67), lat(76,71), lat(68,75), lat(63,78), lat(63,79), lat(61,79), lat(60,78), lat(56,78), lat(56,76), lat(55,76), lat(55,74), lat(56,73), lat(56,72), lat(54,72), lat(54,71), lat(53,71), lat(53,68), lat(52,67), lat(51,66), lat(51,64), lat(50,63),lat(50,61), lat(49,60)]; 
%     plot(poly_lat, poly_lon);
    
    l = legend([s], 'ACROBAT Profile', 'Phytoplankton Patch', 'Phytoplankton Patch Edge', 'Palmer Station', 'Wauwermans Islands', 'Joubin Islands', 'Stationary Glider', 'Location', 'SouthWest');
%     legend('Palmer Station', 'Wauwermans Islands', 'Joubin Islands', 'ACROBAT transect','Location', 'NorthEast');


    
    project_mercator;
    set(gca, 'Color', [0.6843 0.8157 0.9882]);
    hold on;

%% export data for stats tests

LCS_obs_null.rand_ftle_patch_avg = random_resample_all_percentage.ftle_patch_avg(a_ind_med_large)';
LCS_obs_null.obs_ftle_patch_avg = patchID_Adelie.ftle_patch_avg(a_ind_med_large)';
LCS_obs_null.rand_ftle_edge_avg = random_resample_all_percentage.ftle_edge_avg';
LCS_obs_null.obs_ftle_edge_avg = patchID_Adelie.ftle_edge_avg';
LCS_obs_null.rand_rpd_patch_avg = random_resample_all_percentage.rpd_patch_avg(a_ind_med_large)';
LCS_obs_null.obs_rpd_patch_avg = patchID_Adelie.rpd_patch_avg(a_ind_med_large)';
LCS_obs_null.rand_rpd_edge_avg = random_resample_all_percentage.rpd_edge_avg';
LCS_obs_null.obs_rpd_edge_avg = patchID_Adelie.rpd_edge_avg';
LCS_obs_null.diff_ftle_patch = diff_ftle_med_large';
LCS_obs_null.diff_ftle_edge = diff_edge_ftle';
LCS_obs_null.diff_ftle_strat_patch = ftle_strat_patch_diff';
LCS_obs_null.diff_ftle_strat_edge = ftle_strat_edge_diff';
LCS_obs_null.diff_ftle_mixed_patch = ftle_mixed_patch_diff';
LCS_obs_null.diff_ftle_mixed_edge = ftle_mixed_edge_diff';


% code modified from 'resample_rpd_edges.m' and 'resample_ftle.m' to use the same resampling
% method but for the "fake" randomly rotated/translated ACROBAT data
% created by 'acrobat_background.m' and saved as a struct
% 'random_resample.m'

% modified from 'resample_rpd_ftle_background_acro.m' to perform the same
% analysis for the patchID from a 5m integration

% outfitted for baffin
%% Load in randomly created ACROBAT grid

load '/home/jmv208/ACROBAT/random_resample5m.mat'
random_resample_all5m = random_resample5m;

%%
% gridded RPD from "autocor_timescale_rpd.m" --> values not coords
load '/home/jmv208/SWARM_data/RPD_gridded_season.mat'
load '/home/jmv208/SWARM_data/RPD_coordinates.mat'


load('/home/jmv208/SWARM_data/CODAR_zeroed.mat');


for i = 1:length(CODAR_zeroed.time)
    day_sum(i) = sum(CODAR_zeroed.u(:,:,i), 'all');
end
ind = find(day_sum ~=0);
start_codar = CODAR_zeroed.dnum(ind(1));
end_codar = CODAR_zeroed.dnum(ind(end))-3;
% this time index is the same one I used to create trajectories and
% calculate rpd with.

%% loop through ACRO observations, assign RPD amount to each profile

random_resample_all5m.rpd_match = NaN(size(random_resample_all5m.time));
rpd_2d = reshape(part_dens_gridded, [667,1511]);

counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((random_resample_all5m.time>=datenum(start))&(random_resample_all5m.time<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = random_resample_all5m.lat(ind);
        lon = random_resample_all5m.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(RPD_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = rpd_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            random_resample_all5m.rpd_match(ind_prof) = nanmean(value);
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end





%% loop through ACRO observations, pair to ftle
% load ftle with binary LCS/no LCS variable
load '/home/jmv208/LCS-Tool-master/demo/ocean_dataset/ftle_LR_binary.mat'

start_codar = min(ftle_LR.time);
end_codar = max(ftle_LR.time);

% ECO_all.rpd_match = NaN(size(ECO_all.midtime));
random_resample_all5m.ftle_match = NaN(size(random_resample_all5m.time));
% create grid

ftle_2d = reshape(ftle_LR.ftle, [4400, 1273]);
[X,Y] = meshgrid(ftle_LR.x, ftle_LR.y);
X = reshape(X, [4400,1]);
Y = reshape(Y, [4400,1]);
ftle_coords = [X,Y];

counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((random_resample_all5m.time>=datenum(start))&(random_resample_all5m.time<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = random_resample_all5m.lat(ind);
        lon = random_resample_all5m.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(ftle_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = ftle_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            random_resample_all5m.ftle_match(ind_prof) = nanmean(value);
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end

save('random_resample_all5m' , 'random_resample_all5m');

exit

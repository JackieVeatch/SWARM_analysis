% code modified from 'resample_rpd_edges.m' and 'resample_ftle.m' to use the same resampling
% method but for the "fake" randomly rotated/translated ACROBAT data
% created by 'acrobat_background.m' and saved as a struct
% 'random_resample.m'

% outfitted for baffin
%% Load in randomly created ACROBAT grid

load '/home/jmv208/ACROBAT/random_resample.mat'

%% create a more useable format

random_resample_all.lon = horzcat(random_resample(1).lon, random_resample(2).lon, random_resample(3).lon,... 
    random_resample(4).lon, random_resample(5).lon, random_resample(6).lon, random_resample(7).lon, ...
    random_resample(8).lon, random_resample(9).lon, random_resample(10).lon, random_resample(11).lon, ...
    random_resample(12).lon, random_resample(13).lon, random_resample(14).lon, random_resample(15).lon,...
    random_resample(16).lon);

random_resample_all.lat = horzcat(random_resample(1).lat, random_resample(2).lat, random_resample(3).lat,... 
    random_resample(4).lat, random_resample(5).lat, random_resample(6).lat, random_resample(7).lat, ...
    random_resample(8).lat, random_resample(9).lat, random_resample(10).lat, random_resample(11).lat, ...
    random_resample(12).lat, random_resample(13).lat, random_resample(14).lat, random_resample(15).lat,...
    random_resample(16).lat);

random_resample_all.time = horzcat(random_resample(1).time, random_resample(2).time, random_resample(3).time,... 
    random_resample(4).time, random_resample(5).time, random_resample(6).time, random_resample(7).time, ...
    random_resample(8).time, random_resample(9).time, random_resample(10).time, random_resample(11).time, ...
    random_resample(12).time, random_resample(13).time, random_resample(14).time, random_resample(15).time,...
    random_resample(16).time);

random_resample_all.patch = horzcat(random_resample(1).patch, random_resample(2).patch, random_resample(3).patch,... 
    random_resample(4).patch, random_resample(5).patch, random_resample(6).patch, random_resample(7).patch, ...
    random_resample(8).patch, random_resample(9).patch, random_resample(10).patch, random_resample(11).patch, ...
    random_resample(12).patch, random_resample(13).patch, random_resample(14).patch, random_resample(15).patch,...
    random_resample(16).patch);

random_resample_all.edges = horzcat(random_resample(1).edges, random_resample(2).edges, random_resample(3).edges,... 
    random_resample(4).edges, random_resample(5).edges, random_resample(6).edges, random_resample(7).edges, ...
    random_resample(8).edges, random_resample(9).edges, random_resample(10).edges, random_resample(11).edges, ...
    random_resample(12).edges, random_resample(13).edges, random_resample(14).edges, random_resample(15).edges,...
    random_resample(16).edges);

random_resample_all.patchID = horzcat(random_resample(1).patchID, random_resample(2).patchID, random_resample(3).patchID,... 
    random_resample(4).patchID, random_resample(5).patchID, random_resample(6).patchID, random_resample(7).patchID, ...
    random_resample(8).patchID, random_resample(9).patchID, random_resample(10).patchID, random_resample(11).patchID, ...
    random_resample(12).patchID, random_resample(13).patchID, random_resample(14).patchID, random_resample(15).patchID,...
    random_resample(16).patchID);

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

random_resample_all.rpd_match = NaN(size(random_resample_all.time));
rpd_2d = reshape(part_dens_gridded, [667,1511]);

counter = 1;
for i=start_codar:hours(1):end_codar
    
    start = datenum(i - minutes(30));
    finish = datenum(i + minutes(30));
    ind = find((random_resample_all.time>=datenum(start))&(random_resample_all.time<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = random_resample_all.lat(ind);
        lon = random_resample_all.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(RPD_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = rpd_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            random_resample_all.rpd_match(ind_prof) = nanmean(value);
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
random_resample_all.ftle_match = NaN(size(random_resample_all.time));
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
    ind = find((random_resample_all.time>=datenum(start))&(random_resample_all.time<=datenum(finish)));
    if isempty(ind)
        
    else
        lat = random_resample_all.lat(ind);
        lon = random_resample_all.lon(ind);
        point = NaN([length(ind),2]);
        point(:,1) = lon;
        point(:,2) = lat;
%         k = dsearchn(RPD_coords,point);
        [cIdx,cD] = knnsearch(ftle_coords,point,'K',1,'Distance','chebychev');
        for j = 1:(length(ind))
            value = ftle_2d(cIdx(j,:), counter);
            ind_prof = ind(j);
            random_resample_all.ftle_match(ind_prof) = nanmean(value);
%             gather_index
        end
        
%         value = rpd_2d(ind,counter);
%         ECO_all.rpd_match(ind) = value;
    end
    counter = counter +1;
end

save('random_resample_all' , 'random_resample_all');

exit

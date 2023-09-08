%% Putting newly processed ACROBAT data into workable .mat file
    %%% calculating MLD, Int ML chl, Int ML backscatter, ML chl conc., ML backscatter conc. 
    %%% Int chl, and Int backscatter
% data created using "AcrobatPostProcessingWorksheet_Jan2020.m" and the
% _Feb2020 and _Mar2020 versions found in
% Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/data_reprocessed/ACROBAT-main
% jveatch May2022
% edited 6Jun2023 to calculated int backscatter to 5m

% load in data
Jan = load('/Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/data_reprocessed/raw/ForJack/DATA/Jan2020/PROCESSED/gridded.mat');
Jan = Jan.gridded;
Feb = load('/Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/data_reprocessed/raw/ForJack/DATA/Feb2020/PROCESSED/gridded.mat');
Feb = Feb.gridded;
Mar = load('/Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/data_reprocessed/raw/ForJack/DATA/Mar2020/PROCESSED/gridded.mat');
Mar = Mar.gridded;
%% create one structure
ACRO.z = Jan.z;
ACRO.p = Jan.p;
ACRO.lat = [NaN, Jan.lat, Feb.lat, Mar.lat, NaN]; %lat and lon are two indices shorter than rest of data...
ACRO.lat = -ACRO.lat;
ACRO.lon = [NaN, Jan.lon, Feb.lon, Mar.lon, NaN];
ACRO.dist = [Jan.dist, Feb.dist, Mar.dist];
ACRO.mtime = [Jan.mtime, Feb.mtime, Mar.mtime];
ACRO.updown = [Jan.updown, Feb.updown, Mar.updown];
ACRO.t = [Jan.t, Feb.t, Mar.t];
ACRO.s = [Jan.s, Feb.s, Mar.s];
ACRO.dens = [Jan.dens, Feb.dens, Mar.dens];
ACRO.particle = [Jan.particle, Feb.particle, Mar.particle];
ACRO.CDOM = [Jan.CDOM, Feb.CDOM, Mar.CDOM];
ACRO.chl = [Jan.chl, Feb.chl, Mar.chl];

%% calculate MLD

% path to seawater toolbox
addpath '/Users/jveatch/Documents/MATLAB/seawater_ver3_3'
ACRO.bfrq = nan(size(ACRO.t));

for i = 1:length(ACRO.lat)
    sal = ACRO.s(:,i);
    temp = ACRO.t(:,i);
    pres = ACRO.p;
    lat = ACRO.lat(i);
    ACRO.bfrq(:,i) = vertcat( nan, sw_bfrq(sal, temp, pres, lat)); % no bfrq at surface--> bfrq calculated between two depth bins
    mld_ind = find( ACRO.bfrq(:,i) == max(ACRO.bfrq(:,i)));
    mld = ACRO.z(mld_ind);

    if length(sal(isnan(ACRO.bfrq(:,i)))) > 60 % flag if less than 6 depth bins with data
        ACRO.mld(i) = NaN;
    else
        if mld < 3 % if surface noise is max(bfrq), assign mld to second highest bfrq
            bfrq = sort(ACRO.bfrq(:,i), 'descend');
            bfrq = bfrq(~isnan(bfrq));
            bfrq(1) = NaN;
            mld_ind = find(ACRO.bfrq(:,i) == bfrq(2));
            ACRO.mld(i) = ACRO.z(mld_ind);

            if mld < 3 % if surface noise is still max(bfrq), assign mld to deepest point
                present_ind = ~isnan(ACRO.t(:,i));
                ACRO.mld(i) = max(ACRO.z(present_ind)); % set to deepest point
            else
            end

        elseif max(ACRO.bfrq(:,i)) < 0.00005 % if bfrq doesn't meet this low threshold, MLD is likely below ACRO reach
            max_ind = ~isnan(ACRO.t(:,i));
            ACRO.mld(i) = max(ACRO.z(max_ind)); % set to deepest point
            if max(ACRO.z(max_ind)) < 3
                ACRO.mld(i) = NaN;
            else
            end
        else
            bfrq = sort(ACRO.bfrq(:,i), 'descend');
            ACRO.mld(i) = mld;
        end
        
        bfrq = bfrq(~isnan(bfrq));
        dif = abs(bfrq(1) - bfrq(2));

        if ACRO.mld(i) > 10 && dif < 0.000005 % in a three layered system, we want to take the shallowest transition layer
            ind1 = find(ACRO.bfrq(:,i) == bfrq(1));
            ind2 = find(ACRO.bfrq(:,i) == bfrq(2));
            depth(1) = ACRO.z(ind1);
            depth(2) = ACRO.z(ind2);
            if min(depth) > 3
                ACRO.mld(i) = min(depth);
            else
            end
        end
    end
end

%% calculate integrated ML Chl and ML Backscatter

ACRO.mlchl = nan(size(ACRO.mld));
ACRO.mlparticle = nan(size(ACRO.mld)); % Hank's code calls backscatter "particle", so that's what I'm going to do here
ACRO.mlchl_conc = nan(size(ACRO.mld));
ACRO.mlparticle_conc = nan(size(ACRO.mld));

for i = 1:length(ACRO.lat)
    
    mld = ACRO.mld(i);
    if isnan(mld)
        ACRO.mlchl(i) = NaN;
        ACRO.mlparticle(i) = NaN;
        ACRO.mlchl_conc(i) = NaN;
        ACRO.mlparticle_conc(i) = NaN;
    else
        ind_ml = find(ACRO.z<=mld);
        chl = ACRO.chl(:,i);
        particle = ACRO.particle(:,i);
        ACRO.mlchl(i) = nansum(chl(ind_ml));
        ACRO.mlparticle(i) = nansum(particle(ind_ml));
        ACRO.mlchl_conc(i) = nansum(chl(ind_ml))/length(ind_ml);
        ACRO.mlparticle_conc(i) = nansum(particle(ind_ml))/length(ind_ml);
    end
end

ACRO.intchl = nan(size(ACRO.mld));
ACRO.intparticle = nan(size(ACRO.mld));
ACRO.chl_conc = nan(size(ACRO.mld));
ACRO.particle_conc = nan(size(ACRO.mld));

for i = 1:length(ACRO.lat)

    if length(sal(isnan(ACRO.chl(:,i)))) > 60 % at least 6 datapoints
        ACRO.intchl(i) = NaN;
        ACRO.intparticle(i) = NaN;
    else
        ACRO.intchl(i) = nansum(ACRO.chl(:,i));
        ACRO.intparticle(i) = nansum(ACRO.particle(:,i));
        full_bins = find(~isnan(ACRO.chl(:,i)));
        ACRO.chl_conc(i) = nansum(ACRO.chl(:,i))/length(full_bins);
        ACRO.particle_conc(i) = nansum(ACRO.particle(:,i))/length(full_bins);
    end
end

ACRO.five_particle = nan(size(ACRO.mld));
ind_surf = find(ACRO.z <= 5);

for i = 1:length(ACRO.lat)
    if length(sal(isnan(ACRO.chl(:,i)))) > 60 % at least 6 datapoints
        ACRO.five_particle(i) = NaN;
    elseif length(sal(isnan(ACRO.particle(ind_surf,i)))) > 3 % at least 3 bins in surface
        ACRO.five_particle(i) = NaN;
    else
        ACRO.five_particle(i) = nansum(ACRO.particle(ind_surf,i));
        ACRO.intparticle(i) = nansum(ACRO.particle(:,i));
    end
end

% save('ACRO_reprocessed.mat', 'ACRO');
save('ACRO_reprocessed5m.mat', 'ACRO');

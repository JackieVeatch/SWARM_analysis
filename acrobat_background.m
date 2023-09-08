% Creates random resampling of LCS field using translated and rotated
% ACROBAT data. For each survey day, 100 new surveys are created and saved
% to "random_resample" structure
% Dec2022 Jackie Veatch
% edited to run in baffin
% 17Jan2022 edited to not include edges of domain where LCS are no good


load('/home/jmv208/ACROBAT/patchID_Adelie_matched.mat');
load('/home/jmv208/ACROBAT/LCS_bounds.mat');

bounds = LCS_bounds.Vertices;
lcs_lon = bounds(:,1);
lcs_lat = bounds(:,2);

survey_days = ['15-Jan-2020'; '18-Jan-2020'; '21-Jan-2020'; '24-Jan-2020'; 
    '28-Jan-2020'; '01-Feb-2020';'05-Feb-2020';'07-Feb-2020';
    '12-Feb-2020'; '14-Feb-2020'; '18-Feb-2020'; '21-Feb-2020';'22-Feb-2020'; 
    '25-Feb-2020'; '28-Feb-2020'; '03-Mar-2020' ];

starting_patch = 0;

for i = 1:length(survey_days)

    fprintf('starting new survey day\n');
    day = survey_days(i,:);
    ind = find(patchID_Adelie.timestamp >= datenum(day) & patchID_Adelie.timestamp <= datenum(day)+1);
    ind = ind(2:end);
    lon = patchID_Adelie.lon(ind);
    lat = patchID_Adelie.lat(ind);
    time = patchID_Adelie.timestamp(ind);
    patchID = patchID_Adelie.patch_noNaN(ind); % filled so that 1 NaN 1 becomes 1 1 1, but 1 NaN 0 remains
    edgeID = patchID_Adelie.edges(ind);
    edges_thresh = patchID_Adelie.edges_threshold(ind);
    patchNum = patchID_Adelie.patchID(ind) + starting_patch;
    starting_patch = max(patchNum);
%     len_patch = patchID_Adelie.len_patch(ind); % do this piece when calculating % match
    
    random_resample(i).lon = [];
    random_resample(i).lat = [];
    random_resample(i).time = [];
    random_resample(i).patch = [];
    random_resample(i).edges = [];
    random_resample(i).patchID = [];
    random_resample(i).edges_thresh = [];

    c = 0;
    
    while c < 100

        % define the x- and y-data for the original line we would like to rotate
        x = lon;
        y = lat;

        % x random translate
        a = -0.3;
        b = 0.3;
        N = 1;
        x_trans = a + (b-a).*rand(N,1); % randomize angle between a and b values
        x = x+ x_trans;

        % y random translate
        a = -0.3;
        b = 0.3;
        N = 1;
        y_trans = a + (b-a).*rand(N,1); % randomize angle between a and b values
        y = y + y_trans;

        % create a matrix of these points, which will be useful in future calculations
        v = [x;y];
        
        % choose a point which will be the center of rotation
        x_center = median(x);
        y_center = median(y);
        
        % create a matrix which will be used later in calculations
        center = repmat([x_center; y_center], 1, length(x));
        % define a 60 degree counter-clockwise rotation matrix
        % theta = pi/3;       % pi/3 radians = 60 degrees
        a = 0;
        b = 2*pi;
        N = 1;
        theta = a + (b-a).*rand(N,1); % randomize angle between a and b values
        % theta = 0;

        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        % do the rotation...
        s = v - center;     % shift points in the plane so that the center of rotation is at the origin
        so = R*s;           % apply the rotation about the origin
        vo = so + center;   % shift again so the origin goes back to the desired center of rotation
        % this can be done in one line as:
        % vo = R*(v - center) + center
        % pick out the vectors of rotated x- and y-data
        x_rotated = vo(1,:);
        y_rotated = vo(2,:);
        % make a plot
%         scatter(lon, lat);
%         hold on;
%         scatter(x,y);
%         scatter(x_rotated,y_rotated);
%         plot(median(x), median(y), 'r*');

%         lcs_lon = [-63.99, -63.94, -63.94, -64.04, -64.20, -64.35, -64.40, -64.41, -64.47, -64.49, -64.48, -64.47, -64.32, -64.17, -64.00, -63.99];
%         lcs_lat = [-64.90, -64.85, -64.80, -64.76, -64.74, -64.75, -64.79, -64.83, -64.85, -64.89, -64.91, -65.01, -65.05, -65.02, -64.95, -64.90];
%         plot(lcs_lon, lcs_lat);

%         legend('original', 'translated', 'translated and rotated', 'rotation point', 'LCS bounds');

%         title('Randomly resampled ACROBAT data');

        in = inpolygon(x_rotated, y_rotated, lcs_lon, lcs_lat);

        if sum(in) == length(x_rotated)
            c = c + 1;
            random_resample(i).lon = [random_resample(i).lon, x_rotated];
            random_resample(i).lat = [random_resample(i).lat, y_rotated];
            random_resample(i).time = [random_resample(i).time, time];
            random_resample(i).patch = [random_resample(i).patch, patchID];
            random_resample(i).edges = [random_resample(i).edges, edgeID];
            random_resample(i).patchID = [random_resample(i).patchID, patchNum];
            random_resample(i).edges_thresh = [random_resample(i).edges_thresh, edges_thresh];
%             random_resample(i).len_patch = [random_resample(i).edges, len_patch];

        else
        end

    end

end

save('random_resample' , 'random_resample');

exit


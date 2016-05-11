function [estPaths] = multitrack2D(X,Y,Lr,Lc,numBeams)
%% System parameters
dt          = 1; % Sampling interval
startFrame  = 1; % Starting frame
u           = 0; % Control input
noise       = .1; % process noise intensity
noise_x     = 1;  % measurement noise in the horizontal direction (x axis).
noise_y     = 1;  % measurement noise in the horizontal direction (y axis).
%% Kalman parameters
R   = [noise_x 0; ...
        0 noise_y]; % Covariance of the noise
Q   = [dt^4/4 0 dt^3/2 0; ...
        0 dt^4/4 0 dt^3/2; ...
        dt^3/2 0 dt^2 0; ...
        0 dt^3/2 0 dt^2].*noise^2; % Covariance of the observation noise
P   = Q; % Estimate of initial state
A   = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; % State transition model
B   = [(dt^2/2); (dt^2/2); dt; dt]; % Control-input model
H   = [1 0 0 0; 0 1 0 0];  % Observation model
%% Multi tracking parameters
beamObserv      = [X{startFrame} Y{startFrame} zeros(length(X{startFrame}),1) zeros(length(X{startFrame}),1)]';
beamEstimation  = nan(4,2000);
beamEstimation(:, 1:size(beamObserv, 2)) = beamObserv;  % Initial estimate
beamXestimation = nan(2000); % X estimate
beamYestimation = nan(2000); % Y estimate
P_est           = P; % Initial covariance matrix
trackLost       = zeros(1,2000); % How many times a track was not assigned
numDetect       = size(X{startFrame},1); % Initial number of detections
numBeam         = find(isnan(beamEstimation(1, :)) == 1, 1) - 1; % Initial number of track estimates

%% Start the multi-tracking 
for t = startFrame:length(X)
    beamMeasurement = [X{t} Y{t}]; % Matrix with current measurements
    %% Perform Kalman Filter
    % Time Update (Prediction of state for all the beams)
    numDetect = size(X{t},1); % How many detections in current time
    for beam = 1:numBeam
        beamEstimation(:,beam) = A * beamEstimation(:,beam) + B * u;
    end
    P = A * P* A' + Q;
    
    %% Perform Hungarian Algorithm
    % Create the distance matrix between all the detections
    % For the matrix, it is assigned: rows = tracks & columns = detections
    distMatrix = pdist([beamEstimation(1:2, 1:numBeam)'; beamMeasurement]);
    distMatrix = squareform(distMatrix); % Create the squared distance matrix
    distMatrix = distMatrix(1:numBeam, numBeam+1:end) ; % Do only matching for the tracks detected
    
    [assignment, cost] = assignmentoptimal(distMatrix); % Hungarian Algorithm
    assignment = assignment';
    
    % Check exceptions where matching must be ignored and just estimate
    % In those casses assignment = 0    
    % Detection far from observation
    rejected = [];
    for beam = 1:numBeam
        if assignment(beam) > 0
            rejected(beam) =  distMatrix(beam,assignment(beam)) < 50 ;
        else
            rejected(beam) = 0;
        end
    end
    assignment = assignment .* rejected;
    % Done with matching
    
    % Measurement Update (Correction of state for all the beams)
    K = P * H' / ( H * P * H' + R);
    k = 1;
    for beam = 1:length(assignment)
        if assignment(beam) > 0
            beamEstimation(:, k) = beamEstimation(:, k) + K * (beamMeasurement(assignment(beam), :)' - H * beamEstimation(:, k));
        end
        k = k + 1;
    end
    P = (eye(4)- K * H) * P; % Update error covariance
    
    %% Store data
    beamXestimation(t,1:numBeam) = beamEstimation(1,1:numBeam);
    beamYestimation(t,1:numBeam) = beamEstimation(2,1:numBeam);
    
    % Assigning new detections and lost trackings
    % For new detections: If it wasn't assigned means a new beam
    newTracks = beamMeasurement(~ismember(1: size(beamMeasurement, 1), assignment), :)';
    if ~isempty(newTracks)
        beamEstimation(:, numBeam + 1:numBeam + size(newTracks, 2)) = ...
            [newTracks; zeros(2, size(newTracks, 2))];
        % Number of estimated beams including new ones
        numBeam = numBeam + size(newTracks, 2);
    end
    
    % If a tracking didn't get matched with a detection, 
    % a counter will start
    noTrackInList = find(assignment == 0);
    if ~isempty(noTrackInList)
        trackLost(noTrackInList) = trackLost(noTrackInList) + 1;
    end
    
    % If a track has a counter greater than 6, the tracking will be deleted
    % and reseted to NaN
    bad_trks = find(trackLost > 6);
    beamEstimation(:,bad_trks) = NaN;
    
    %% Visualization
    %{
    clf
    img = ones(500, 500, 3); % Create a blank image for visualization
    imagesc(img);
    hold on;
    plot(Y{t}(:),X{t}(:),'or'); % Plot measurements
    colours = ['r', 'b', 'g', 'c', 'm', 'k'];
    for nB = 1:numBeam
        if ~isnan(beamXestimation(t,nB))
            cIdx = mod(nB,6) + 1; %pick color
            tempX = beamXestimation(1:t,nB);
            tempY = beamYestimation(1:t,nB);
            plot(tempY, tempX, '.-', 'MarkerSize', 3, 'Color', colours(cIdx), 'LineWidth',3)
            axis off
        end
    end
    pause(0.05)
    %}
%     t
end

% Creating the estimated paths for each beam
for i=1:numBeams
    estPaths{i} = [beamXestimation(1:length(X),i) beamYestimation(1:length(Y),i)];
end
    
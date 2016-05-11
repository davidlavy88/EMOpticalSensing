function [estPath] = kalman2D(observ, Lr, Lc)
%% System parameters
dt      = 1; % sampling interval
t       = 1; % starting frame
u       = .005; % control input
x_init  = [observ(t,1); observ(t,2); 0; 0]; % Initial Conditions
x_est   = x_init; % state estimate
noise   = .1; % process noise intensity
noise_x = 1; % noise for x and y are
noise_y = 1; % chosen by user and the same
visualize = 0; % visualize the tracking
numObserv = length(observ);

%% Kalman Filter params
R       = [noise_x 0; ...
            0 noise_y]; %coviarance of the noise
Q       = [dt^4/4 0 dt^3/2 0; ...
            0 dt^4/4 0 dt^3/2; ...
            dt^3/2 0 dt^2 0; ...
            0 dt^3/2 0 dt^2] .* noise^2; % Covariance of the observation noise
P       = Q; % Estimate of initial state
A       = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];  %State transition model
H       = [1 0 0 0; 0 1 0 0];  % Observation model
B       = [(dt^2/2); (dt^2/2); dt; dt];  % Control-input model
z       = []; % The measurements of the node state
x_sta_est = []; %  Initial state estimate
v_est   = []; % Initial velocity estimate
P_est   = P; % Initial covariance matrix

%% Perform Kalman Filter
% figure
for i = t:numObserv
    img     = ones(Lr, Lc, 3); % Create a blank image for visualization
    z(i,:)  = [observ(i,1) observ(i,2)]; % Current measurement coordinates
    
    % Time Update
    x_est   = A * x_est + B * u; % Project the state ahead
    P       = A * P * A' + Q; % Project the error covariance ahead
    % Measurement Update
    K       = P * H' / (H * P * H' + R); % Compute the Kalman Gain
    if ~isnan(z(i,:))
        x_est = x_est + K * (z(i,:)' - H * x_est); % Update estimate with measurement
    end   
    P       =  (eye(4) - K * H) * P; % Update error covariance
    
    x_sta_est   = [x_sta_est; x_est(1:2)'];
    v_est       = [v_est; x_est(3:4)'];
    
    x_estimation(i) = x_est(1); %estimation in horizontal position
    y_estimation(i) = x_est(2); %estimation in vertical position
    
    if visualize==1
        r = 5;
        ang=0:.01:2*pi; %parameters of nodes
        imagesc(img);
        set(gca,'YDir','normal')
    %     axis off
        hold on
        plot(r * cos(ang) + ground(i,1), r * sin(ang) + ground(i,2), '.g'); % Ground truth motion
        plot(r * cos(ang) + z(i,1), r * sin(ang) + z(i,2), '.b'); % The measurement motion
        plot(r * cos(ang) + x_est(1), r * sin(ang) + x_est(2), '.r'); % The kalman filtered motion
        hold off
        legend('Ground truth', 'Measurement', 'Kalman Filter')
        pause(0.05) 
    end
end

x_estimation = x_estimation';
y_estimation = y_estimation';

estPath = [x_estimation y_estimation];

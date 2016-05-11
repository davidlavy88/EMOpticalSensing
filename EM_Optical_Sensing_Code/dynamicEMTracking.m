%% Single beam tracking using static and dynamic EM
clear; close all; clc;
%% Create the path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure you have a path saved. To do this type:
% >> imagesc(ones(500,500,3)), grid on
% >> imfreehand % Draw a random path
% Right click on the drawing and select copy positions, then paste it on 
% the command window and save it as path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For a sample path, uncomment and run the following
load path
%% Static Data
Lr = 500; Lc = 500;
rho = 40;
Lam_s = 50;
Lam_n = 50;
%% Tracking Motion
numFrames = size(path,1);

[sig_pos, matDetect, listDetect,labels] = ...
    fcn_generate_motion_data(Lr,Lc,rho,Lam_s,Lam_n,path);

% Apply staticEM for the first frame
xhats = staticEM(matDetect(:,:,1),listDetect{1},rho,Lam_s,Lam_n,20);
xest = xhats{end};
measurements(1,:) = xest;
% Apply dynamicEM for the second frame
sigma_prior = 3; % This value should be around the number of pixels the beam can move
prev_xhat = xest;
% numFrames = 358;
%%
for ii = 2:numFrames-1
    xhat = dynamicEM(matDetect(:,:,ii),listDetect{ii}, ...
        prev_xhat,rho,Lam_s,Lam_n);
    prev_xhat = xhat;
    measurements(ii,:) = prev_xhat;
end

[estPath] = kalman2D(measurements, Lr, Lc);
%%
%%{
figure; 
for ii = 1:numFrames-1
    imagesc(matDetect(:,:,ii)); axis image; colormap(gray);
    hold on;
    plot(sig_pos(ii,2),sig_pos(ii,1),'g*',...
     'MarkerSize',10,'LineWidth',3)
    plot(measurements(ii,2),measurements(ii,1),'rx',...
     'MarkerSize',10,'LineWidth',3)
    plot(estPath(ii,2),estPath(ii,1),'bo',...
     'MarkerSize',10,'LineWidth',3)
    legend('Ground Truth','EM Estimate', 'Kalman Filter')
    title('Single-beam tracking using dynamic EM and Kalman filter')
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uncomment below for creating a video
    filename = [sprintf('%03d',ii) '.jpg'];
    fullname = fullfile('images',filename);
    saveas(gcf, fullname);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause(0.05)
end
%}

figure
imagesc(matDetect(:,:,35)), plot(sig_pos(:,1),sig_pos(:,2),'g','LineWidth',2), hold on
plot(measurements(:,1),measurements(:,2),'r','LineWidth',2);
plot(estPath(:,1),estPath(:,2),'b','LineWidth',2);
legend('Ground Truth','EM Estimate', 'Kalman Filter')

%% Create video
shuttleVideo = VideoReader('shuttle.avi');
imageNames = dir(fullfile('images','*.jpg'));
imageNames = {imageNames.name}';
outputVideo = VideoWriter(fullfile('images','dynamicEM_Kalman.avi'));
outputVideo.FrameRate = shuttleVideo.FrameRate;
open(outputVideo)
for ii = 1:length(imageNames)
   img = imread(fullfile('images',imageNames{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo)

%% Statistics
% Difference in the X axis between ground truth and estimated position
n=1:numFrames-1;
figure
subplot(211)
plot(n, sig_pos(1:end-1,1), 'b', n, estPath(:,1), 'r');
xlabel('Time'); ylabel('Position in X axis');
title('Position offset in horizontal direction');
legend('Ground truth', 'KF Estimated')

% Difference in the Y axis between ground truth and estimated position
% figure
subplot(212)
plot(n, sig_pos(1:end-1,2), 'b', n, estPath(:,2), 'r');
xlabel('Time'); ylabel('Position in Y axis');
title('Position offset in vertical direction');
legend('Ground truth', 'KF Estimated')

% Distance error between ground truth and estimated position
distX = sig_pos(1:end-1,1) - estPath(:,1);
distY = sig_pos(1:end-1,2) - estPath(:,2);
distTotal = sqrt(distX.^2 + distY.^2);
figure
plot(n,distTotal / (500*sqrt(2)) * 100);
xlabel('Time'); ylabel('Error distance (%)');
title('Error distance between actual and estimation positions');

% Longest offset from estimate and ground truth
maxDistance = max(distTotal);

% MSE and MAE
dynamicMSE = (1/(numFrames-1)) * sum(distX.^2 + distY.^2);
dynamicMAE = (1/(numFrames-1)) * sum(abs(distX) + abs(distY));

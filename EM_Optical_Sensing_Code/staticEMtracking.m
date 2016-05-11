%% Test script for data generation
% Boston University
% EC 503

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

for ii = 1:numFrames
    xhats = staticEM(matDetect(:,:,ii),listDetect{ii},rho,Lam_s,Lam_n,20);
    xest = xhats{end};
    measurements(ii,:) = xest;
end

[estPath] = kalman2D(measurements, Lr, Lc);
%% Visualization
%{
figure; 
for ii = 1:numFrames
    imagesc(matDetect(:,:,ii)); axis image; colormap(gray);
    hold on;
    plot(sig_pos(ii,2),sig_pos(ii,1),'g*',...
     'MarkerSize',10,'LineWidth',3)
    plot(measurements(ii,2),measurements(ii,1),'rx',...
     'MarkerSize',10,'LineWidth',3)
    plot(estPath(ii,2),estPath(ii,1),'bo',...
     'MarkerSize',10,'LineWidth',3)
    legend('Ground Truth','EM Estimate', 'Kalman Filter')
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uncomment below for creating a video
%     filename = [sprintf('%03d',ii) '.jpg'];
%     fullname = fullfile('images2',filename);
%     saveas(gcf, fullname);
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
% shuttleVideo = VideoReader('shuttle.avi');
% imageNames = dir(fullfile('images2','*.jpg'));
% imageNames = {imageNames.name}';
% outputVideo = VideoWriter(fullfile('images2','EM_Kalman.avi'));
% outputVideo.FrameRate = shuttleVideo.FrameRate;
% open(outputVideo)
% for ii = 1:length(imageNames)
%    img = imread(fullfile('images2',imageNames{ii}));
%    writeVideo(outputVideo,img)
% end
% close(outputVideo)

%% Statistics
% Difference in the X axis between ground truth and estimated position
n=1:numFrames;
figure
plot(n, sig_pos(:,1), 'b', n, estPath(:,1), 'r');
xlabel('Time'); ylabel('Position in X axis');
title('Position offset in horizontal direction');
legend('Ground truth', 'KF Estimated')

% Difference in the Y axis between ground truth and estimated position
figure
plot(n, sig_pos(:,2), 'b', n, estPath(:,2), 'r');
xlabel('Time'); ylabel('Position in Y axis');
title('Position offset in vertical direction');
legend('Ground truth', 'KF Estimated')

% Distance error between ground truth and estimated position
distX = sig_pos(:,1) - estPath(:,1);
distY = sig_pos(:,2) - estPath(:,2);
distTotal = sqrt(distX.^2 + distY.^2);
figure
plot(n,distTotal);
xlabel('Time'); ylabel('Error distance');
title('Error distance between actual and estimation positions');

% Longest offset from estimate and ground truth
maxDistance = max(distTotal);

% MSE and MAE
staticMSE = (1/numFrames) * sum(distX.^2 + distY.^2);
staticMAE = (1/numFrames) * sum(abs(distX) + abs(distY));

%% Test script for data generation
% Joshua Rapp
% Boston University
% EC 503

clear; close all; clc;

%% Static Data
Lr = 500; Lc = 500;
rho = 40;
Lam_s = 50;
Lam_n = 50;

[ sig_pos, matDetect, listDetect,labels  ] = fcn_generate_data( Lr,Lc,rho,Lam_s,Lam_n  );
centroid = mean(listDetect);
%% Plot Detections as Image
figure; imagesc(matDetect); axis ij image; colormap(gray);
hold on;
plot(sig_pos(2),sig_pos(1),'g+',...
     'MarkerSize',10,'LineWidth',3)

euclid_dist = sqrt(sum((sig_pos-centroid).^2));

%% Apply EM
xhats = staticEM(matDetect,listDetect,rho,Lam_s,Lam_n,20);
xest = xhats{end};
figure; imagesc(matDetect); axis ij image; colormap(gray);
hold on;
plot(sig_pos(2),sig_pos(1),'g+','MarkerSize',10,'LineWidth',3)
plot(xest(2),xest(1),'rx',...
     'MarkerSize',10,'LineWidth',3)

%% Apply k-means
numClusters = 2;
[idx,C,sumd] = kmeans(listDetect,numClusters);
CCR = mean((2-labels)==idx);

cmap = hsv(numClusters);

figure;
plot(listDetect(idx==1,1),listDetect(idx==1,2),'r.','MarkerSize',12)
hold on

for ii = 2:numClusters
    plot(listDetect(idx==ii,1),listDetect(idx==ii,2),'.','Color',cmap(ii,:),'MarkerSize',12);
end

plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
plot(sig_pos(1),sig_pos(2),'g+','MarkerSize',10,'LineWidth',3)
plot(xest(1),xest(2),'mx','MarkerSize',10,'LineWidth',3)
legend('Cluster 1','Cluster 2','Centroids',...
       'True Position','EM Estimate','Location','NW')
title 'Cluster Assignments and Centroids'
hold off

EMerror = norm(sig_pos-xest');
KMeansError = norm(sig_pos-C(1,:));
%% Tracking Motion
% numFrames = 10;
% speed = 40;
% 
% [sig_pos, matDetect, listDetect,labels] = ...
%     fcn_generate_motion_data(Lr,Lc,rho,Lam_s,Lam_n,numFrames,speed);
% 
% for ii = 1:numFrames
%     figure; imagesc(matDetect(:,:,ii)); axis image; colormap(gray);
%     hold on;
%     plot(sig_pos(ii,2),sig_pos(ii,1),'g*',...
%      'MarkerSize',10,'LineWidth',3)
% end

% implay(matDetect,1);
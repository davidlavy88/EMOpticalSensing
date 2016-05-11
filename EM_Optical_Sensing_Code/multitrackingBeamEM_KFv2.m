%% Multitracking beams 
clc, clear all, close all
rng(1);
%% Load data and create path
% Laser
load laserPoints

paths{1} = laser1;
paths{2} = laser2;
%% Parameters
Lr      = 500; 
Lc      = 500;
rho     = 40;
Lam_s   = 50;
Lam_n   = 50;
numBeam = 2;

[sig_pos_all, matDetect_all, listDetect_all,labels_all] = ...
    fcn_generate_motion_data_MB(Lr,Lc,rho,Lam_s,Lam_n,paths);

numFrames = size(paths{1},1);
% numFrames = 20;

% one row per beam in these results
max_error = zeros(2,100);
MSE = zeros(2,100);
MAE = zeros(2,100);
sigma_factors = 0.1:0.1:10;
% for MAsize = 1:length(sigma_factors);
%     varlist = {'measurements','xhats','prev_xhats',...
%         'all_xhats','dists','sorted_dists','sigma_list','sigma',...
%         'Xmb','Ymb','estPaths','dist_X','dist_Y','distTotal',...
%         'maxDistancePerBeam','xdet','ydet','multiMSE','multiMAE'};
%     clear(varlist{:})
    MAsize = 1;

    xhats = multibeamEM(matDetect_all{1},listDetect_all{1},rho,Lam_s,Lam_n,[],[],2);
    prev_xhats = xhats;
    sigma_list = [];
    for ii = 1:numFrames
        N = min(35,ii-1);
        if ii > 10
            all_xhats = cell2mat(measurements);
            all_xhats = reshape(all_xhats',2,2*size(all_xhats,1))';
            dists = pdist(all_xhats(end-N:end,:));
            sorted_dists = sort(dists);
            sorted_dists = sorted_dists(sorted_dists ~= 0);
            sigma = 10*mean(sorted_dists(1:N));
            sigma_list = [sigma_list; sigma];
            xhats = multibeamEM(matDetect_all{ii},listDetect_all{ii},rho,Lam_s,Lam_n,prev_xhats,sigma,2);
        else
            xhats = multibeamEM(matDetect_all{ii},listDetect_all{ii},rho,Lam_s,Lam_n,prev_xhats,[],2);
        end
        prev_xhats = xhats;
        measurements{1}(ii,:) = xhats{1};
        measurements{2}(ii,:) = xhats{2};
        Xmb{ii} = [xhats{1}(1); xhats{2}(1)];
        Ymb{ii} = [xhats{1}(2); xhats{2}(2)];
    end

    % Preparing data for KF-Hungarian
    estPaths = multitrack2D(Xmb,Ymb,Lr,Lc,numBeam);
    %% Visualization
    %%{
    figure; 
    for ii = 1:numFrames
        imagesc(matDetect_all{ii}); axis image; colormap(gray);
        hold on;
        for beam = 1:numBeam
            plot(sig_pos_all{beam}(ii,2),sig_pos_all{beam}(ii,1),'g*',...
            'MarkerSize',10,'LineWidth',3)
            plot(measurements{beam}(ii,2),measurements{beam}(ii,1),'rx',...
            'MarkerSize',10,'LineWidth',3)
            plot(estPaths{beam}(ii,2),estPaths{beam}(ii,1),'bo',...
            'MarkerSize',10,'LineWidth',3)
        end
        legend('Ground Truth','EM Estimate', 'Kalman Filter')
        hold off
        pause(0.001)
    %     pause
    end
    %}

    % match true path and estimated path indices

    estpath11 = estPaths{1}(1,:);
    estpath21 = estPaths{2}(1,:);

    beam11 = sig_pos_all{1}(1,:);
    % beam12 = sig_pos_all{2}(1,:);

    estpath1beam1dist = norm(estpath11 - beam11);
    estpath2beam1dist = norm(estpath21 - beam11);
    % estpath1beam2dist = norm(estpath11 - beam21);
    % estpath2beam2dist = norm(estpath21 - beam21);

    if estpath1beam1dist > estpath2beam1dist
        estPaths = {estPaths{2} estPaths{1}};
    end


%     %%
%     figure, imagesc(ones(Lr,Lc,3))
%     % colours = ['r' 'b' 'g' 'c' 'm' 'y'];
%     for beam=1:numBeam
%     %     cIdx = mod(nB,6)+1; % Pick color
%         plot(sig_pos_all{beam}(:,1), sig_pos_all{beam}(:,2), 'g', 'LineWidth', 2), hold on
%         % Plotting the EM estimates make the plot look bad so avoid it
%     %     plot(measurements{beam}(:,1),measurements{beam}(:,2), 'r', 'LineWidth',2);
%         plot(estPaths{beam}(:,1),estPaths{beam}(:,2),'b','LineWidth',2);
%     end
%     % legend('Ground Truth','EM Estimate', 'Kalman Filter')
%     legend('Ground Truth', 'Kalman Filter')

    %% Statistics
    % Difference in the X axis between ground truth and estimated position
    n=1:numFrames;
%     figure
%     for beam = 1:numBeam
%         subplot(numBeam,1,beam)
%         plot(n, sig_pos_all{beam}(:,1), 'b', n, estPaths{beam}(:,1), 'r');
%         xlabel('Time'); ylabel('Position in X axis');
%     %     title('Position offset in horizontal direction');
%         str = sprintf('Beam %d: Offset y axis', beam);
%         title(str);
%         legend('Ground truth', 'KF Estimated')
%     end

%     % Difference in the Y axis between ground truth and estimated position
%     figure
%     for beam = 1:numBeam
%         subplot(numBeam,1,beam)
%         plot(n, sig_pos_all{beam}(:,2), 'b', n, estPaths{beam}(:,2), 'r');
%         xlabel('Time'); ylabel('Position in Y axis');
%         str = sprintf('Beam %d: Offset Y axis', beam);
%         title(str);
%         legend('Ground truth', 'KF Estimated')
%     end

    % Distance error between ground truth and estimated position
%     figure
    for beam = 1:numBeam
        subplot(numBeam,1,beam)
        distX(:,beam) = sig_pos_all{beam}(:,1) - estPaths{beam}(:,1);
        distY(:,beam) = sig_pos_all{beam}(:,2) - estPaths{beam}(:,2);
        distTotal(:,beam) = sqrt(distX(:,beam).^2 + distY(:,beam).^2);
%         plot(n,distTotal(:,beam));
%         xlabel('Time'); ylabel('Error distance');
%         str = sprintf('Beam %d: Error distance', beam);
%         title(str);
    end

    % Longest offset from estimate and ground truth
    maxDistancePerBeam = max(distTotal);
    max_error(:,MAsize) = maxDistancePerBeam;

    % MSE and MAE
    for beam=1:numBeam
        xdet = distX(:,beam);
        xdet = xdet(~isnan(xdet));
        ydet = distY(:,beam);
        ydet = ydet(~isnan(ydet));
        multiMSE(beam) = (1/length(xdet)) * (sum(xdet.^2) + sum(ydet.^2));
        MSE(beam,MAsize) = multiMSE(beam);
        multiMAE(beam) = (1/length(xdet)) * (sum(abs(xdet)) + sum(abs(ydet)));
        MAE(beam,MAsize) = multiMAE(beam);
    end
% end

% plot errors as function of moving average filter length

% figure
% subplot(2,1,1)
% plot(sigma_factors,max_error(1,:))
% title('Maximum L^2 Distance Between Estimate and True Position (Beam 1)')
% xlabel('Mean Speed Multiplier')
% ylabel('Distance')
% subplot(2,1,2)
% plot(sigma_factors,max_error(2,:))
% title('Maximum L^2 Distance Between Estimate and True Position (Beam 2)')
% xlabel('Mean Speed Multiplier')
% ylabel('Distance')
% 
% figure
% subplot(2,1,1)
% plot(sigma_factors,MAE(1,:))
% set(gca,'YScale','log')
% title('MAE vs. MA Size(Beam 1)')
% xlabel('Mean Speed Multiplier')
% ylabel('MAE')
% subplot(2,1,2)
% plot(sigma_factors,MAE(2,:))
% set(gca,'YScale','log')
% title('MAE vs. MA Size (Beam 2)')
% xlabel('Mean Speed Multiplier')
% ylabel('MAE')
% 
% figure
% subplot(2,1,1)
% plot(sigma_factors,MSE(1,:))
% set(gca,'YScale','log')
% title('MSE vs. MA Size(Beam 1)')
% xlabel('Mean Speed Multiplier')
% ylabel('MSE')
% subplot(2,1,2)
% plot(sigma_factors,MSE(2,:))
% set(gca,'YScale','log')
% title('MSE vs. MA Size (Beam 2)')
% xlabel('Mean Speed Multiplier')
% ylabel('MSE')
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % figure
% % subplot(2,1,1)
% % plot(1:100,max_error(1,:))
% % title('Maximum L^2 Distance Between Estimate and True Position (Beam 1)')
% % xlabel('MA Size')
% % ylabel('Distance')
% % subplot(2,1,2)
% % plot(1:100,max_error(2,:))
% % title('Maximum L^2 Distance Between Estimate and True Position (Beam 2)')
% % xlabel('MA Size')
% % ylabel('Distance')
% % 
% % figure
% % subplot(2,1,1)
% % plot(1:100,MAE(1,:))
% % title('MAE vs. MA Size(Beam 1)')
% % xlabel('MA Size')
% % ylabel('MAE')
% % subplot(2,1,2)
% % plot(1:100,MAE(2,:))
% % title('MAE vs. MA Size (Beam 2)')
% % xlabel('MA Size')
% % ylabel('MAE')
% % 
% % figure
% % subplot(2,1,1)
% % plot(1:100,MSE(1,:))
% % title('MSE vs. MA Size(Beam 1)')
% % xlabel('MA Size')
% % ylabel('MSE')
% % subplot(2,1,2)
% % plot(1:100,MSE(2,:))
% % title('MSE vs. MA Size (Beam 2)')
% % xlabel('MA Size')
% % ylabel('MSE')
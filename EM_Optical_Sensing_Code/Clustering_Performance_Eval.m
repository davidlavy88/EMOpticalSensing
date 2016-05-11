%% Explore EM Performance
% Joshua Rapp
% April 23, 2016

clear; close all; clc;

%%
numS = 3;
numN = 10;
numTrials = 100;

Lr = 500; Lc = 500;
Lam_s = [50,100,500];
rho = 40;
itmax = 100;

meanCOG = zeros(numS,numN);
meanEM = zeros(numS,numN);
meanKMEANS2 = zeros(numS,numN);
meanKMEANS3 = zeros(numS,numN);
meanSPECT = zeros(numS,numN);

for jj = 1:numS
    lam_s = Lam_s(jj);
    disp(num2str(lam_s));
    Lam_n = round(logspace(log10(lam_s/10),log10(10*lam_s),numN));

    for kk = 1:numN
        lam_n = Lam_n(kk); 
        disp(num2str(lam_n));

        EMdist = zeros(numTrials,1);
        COGdist = zeros(numTrials,1);
        K2MEANSdist = zeros(numTrials,1);
        K3MEANSdist = zeros(numTrials,1);
        SPECTdist = zeros(numTrials,1);

        for t = 1:numTrials 
            [ sig_pos, matDetect, listDetect,labels  ] = fcn_generate_data(Lr,Lc,rho,lam_s,lam_n);
            numDetect = length(labels);
            
%             figure; plot(listDetect(labels==1,1),listDetect(labels==1,2),...
%             'r.',listDetect(labels==0,1),listDetect(labels==0,2),'b.');
            
            % EM Estimate
            xhats = staticEM(matDetect,listDetect,rho,lam_s,lam_n,itmax);
            EM_est = xhats{end};
            
            % COG Estimate
            COG_est = mean(listDetect);
            
            % kmeans estimate, k=2
            KMEANS2_est = kmeans_estimate( listDetect, 2);
            
            % kmeans estimate, k=3
            KMEANS3_est = kmeans_estimate( listDetect, 3 );

            % Spectral Clustering
            %SPECT_est = kmeans_estimate( listDetect, 2, 'spectral', 1);
            
            % Euclidean Distance
            EMdist(t) = sqrt((EM_est(1)-sig_pos(1))^2+(EM_est(2)-sig_pos(2))^2);
            COGdist(t) = sqrt((COG_est(1)-sig_pos(1))^2+(COG_est(2)-sig_pos(2))^2);
            K2MEANSdist(t) = sqrt((KMEANS2_est(1)-sig_pos(1))^2+(KMEANS2_est(2)-sig_pos(2))^2);
            K3MEANSdist(t) = sqrt((KMEANS3_est(1)-sig_pos(1))^2+(KMEANS3_est(2)-sig_pos(2))^2);

        end

        meanEM(jj,kk) = mean(EMdist);
        meanCOG(jj,kk) = mean(COGdist);
        meanKMEANS2(jj,kk) = mean(K2MEANSdist);
        meanKMEANS3(jj,kk) = mean(K3MEANSdist);

    end
end

save('slocumb_validation4.mat');

%%
SNRi = 10*log10(fliplr(logspace(-1,1,numN)));

BW = rho*sqrt(log(4));
SNRo_EM = 10*log10(BW^2./meanEM);
SNRoCOG = 10*log10(BW^2./meanCOG);
SNRoKMEANS2 = 10*log10(BW^2./meanKMEANS2);
SNRoKMEANS3 = 10*log10(BW^2./meanKMEANS3);
SNRoSPECT = 10*log10(BW^2./meanSPECT);

for ii = 1:numS
    figure; plot(SNRi,SNRo_EM(ii,:),SNRi,SNRoCOG(ii,:),SNRi,SNRoKMEANS2(ii,:),SNRi,SNRoKMEANS3(ii,:));
    xlabel('Input SNR (dB)');
    ylabel('Output SNR (dB)');
    legend('EM','COG','K2','K3','Location','northwest');
    title(['Performance of EM, Kmeans, and COG Estimators for \Lambda_s = ' num2str(Lam_s(ii))]);
end
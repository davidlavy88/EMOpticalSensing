%% Explore EM Performance
% Joshua Rapp
% April 23, 2016

clear; close all; clc;

%%
numS = 3;
numN = 10;
numTrials = 100;

Lr = 500; Lc = 500;

noise_min = 1;
noise_max = 100;

Lams = [50,100,500];
Lam_s = [];
Lam_n = [];

for ii = 1:numS
    Lam_s = [Lam_s,Lams(ii)*ones(1,numN)];
    Lam_n = [Lam_n,round(logspace(log10(Lams(ii)*noise_min),...
        log10(noise_max*Lams(ii)),numN))];
end

rho = 40;
itmax = 100;

rows = 1:Lr;
cols = 1:Lc;

%% Theoretical Marginals - Uniform (Noise Only)
unif_rows = ones(Lr,1)/Lr;
marg_unif_rows = cumsum(unif_rows)/sum(unif_rows);
unif_cols = ones(Lc,1)/Lc;
marg_unif_cols = cumsum(unif_cols)/sum(unif_cols);

CCRs = zeros(numS*numN,1);

parfor jj = 1:numS*numN
    lam_s = Lam_s(jj);
    lam_n = Lam_n(jj);
    disp(['Sig: ' num2str(lam_s) ' , Noise: ' num2str(lam_n)]);

    true_labels = zeros(numTrials,1);
    predict_labels = zeros(numTrials,1);

    for t = 1:numTrials 
        [ sig_pos, matDetect,listDetect,label] = fcn_generate_distribution(Lr,Lc,rho,lam_s,lam_n);
        true_labels(t) = label;

        % EM Gaussian Center Estimation
        xhats = staticEM(matDetect,listDetect,rho,lam_s,lam_n,itmax);
        xest = xhats{end};

        % Theoretical Marginal - Gaussian + Uniform (Noise and Signal)
        gauss_rows = (lam_s*normpdf(rows,xest(1),rho)+lam_n*unif_rows')/(lam_s+lam_n);
        marg_gauss_rows = cumsum(gauss_rows)/sum(gauss_rows);
        gauss_cols = (lam_s*normpdf(cols,xest(2),rho)+lam_n*unif_cols')/(lam_s+lam_n);
        marg_gauss_cols = cumsum(gauss_cols)/sum(gauss_cols);

        % Compute Marginals
        data_rows = sum(matDetect,2);
        marg_cdf_rows = cumsum(data_rows)/sum(data_rows);
        data_cols = sum(matDetect,1)';
        marg_cdf_cols = cumsum(data_cols)/sum(data_cols);

        MSE_Noise_cols = mean((marg_cdf_cols-marg_unif_cols).^2);
        MSE_Noise_rows = mean((marg_cdf_rows-marg_unif_rows).^2);
        MSE_Signal_cols = mean((marg_cdf_cols-marg_gauss_cols').^2);
        MSE_Signal_rows = mean((marg_cdf_rows-marg_gauss_rows').^2);

        [~,dist_predict] = min([MSE_Noise_cols+MSE_Noise_rows,MSE_Signal_cols+MSE_Signal_rows]);
        predict_labels(t) = dist_predict-1;

    end

    CCRs(jj) = sum(true_labels==predict_labels)/numTrials;
end

save('distribution_detection3.mat');

%%
SNRi = 10*log10(fliplr(logspace(log10(1/noise_max),log10(noise_min),numN)));
CCRs = reshape(CCRs,numN,numS);

figure; plot(SNRi,CCRs);
xlabel('Input SNR (dB)');
ylabel('Detection CCR');
legend('\Lambda_S = 50','\Lambda_S = 100','\Lambda_S = 500','Location','northwest');
title('Distribution Detection');

%% Performance vs SNR, Unknown spot size
% Joshua Rapp
% April 23, 2016

clear; close all; clc;

%%
numS = 3;
numN = 10;
numTrials = 100;

Lr = 500; Lc = 500;
Lams = [50,100,500];

rho_hat1 = 40;
rho_hat2 = 45;
rho_hat3 = 35;

itmax = 100;

meanCOG = zeros(numS,numN)';
meanStaticEM = zeros(numS,numN)';
meanVar1 = zeros(numS,numN)';
meanVar2 = zeros(numS,numN)';
meanVar3 = zeros(numS,numN)';

noise_min = 0.1;
noise_max = 10;

Lams = [50,100,500];
Lam_s = [];
Lam_n = [];

for ii = 1:numS
    Lam_s = [Lam_s,Lams(ii)*ones(1,numN)];
    Lam_n = [Lam_n,round(logspace(log10(Lams(ii)*noise_min),...
        log10(noise_max*Lams(ii)),numN))];
end
%%
parfor jj = 1:numS*numN
    lam_s = Lam_s(jj);
    lam_n = Lam_n(jj); 
    disp(num2str(lam_n));

    StaticDist = zeros(numTrials,1);
    Var1Dist = zeros(numTrials,1);
    Var2Dist = zeros(numTrials,1);
    Var3Dist = zeros(numTrials,1);
    COGdist = zeros(numTrials,1);

    for t = 1:numTrials 
        [sig_pos,~,matDetect,listDetect,labels] = ...
            fcn_generate_correlated_data( Lr,Lc,rho_hat1,lam_s,lam_n );
        numDetect = length(labels);

        % EM Estimate
        xhats = staticEM(matDetect,listDetect,rho_hat1,lam_s,lam_n,itmax);
        EM_est = xhats{end};

        x_var1 = variableEM(matDetect,listDetect,rho_hat1,lam_s,lam_n,itmax);
        Var1_est = x_var1{end};

        x_var2 = variableEM(matDetect,listDetect,rho_hat2,lam_s,lam_n,itmax);
        Var2_est = x_var2{end};

        x_var3 = variableEM(matDetect,listDetect,rho_hat3,lam_s,lam_n,itmax);
        Var3_est = x_var3{end};

        % COG Estimate
        COG_est = mean(listDetect);

        % Euclidean Distance
        StaticDist(t) = sqrt((EM_est(1)-sig_pos(1))^2+(EM_est(2)-sig_pos(2))^2);
        Var1Dist(t) = sqrt((Var1_est(1)-sig_pos(1))^2+(Var1_est(2)-sig_pos(2))^2);
        Var2Dist(t) = sqrt((Var2_est(1)-sig_pos(1))^2+(Var2_est(2)-sig_pos(2))^2);
        Var3Dist(t) = sqrt((Var3_est(1)-sig_pos(1))^2+(Var3_est(2)-sig_pos(2))^2);
        COGdist(t) = sqrt((COG_est(1)-sig_pos(1))^2+(COG_est(2)-sig_pos(2))^2);

    end

    meanStaticEM(jj) = mean(StaticDist);
    meanVar1(jj) = mean(Var1Dist);
    meanVar2(jj) = mean(Var2Dist);
    meanVar3(jj) = mean(Var3Dist);
    meanCOG(jj) = mean(COGdist);


end

save('slocumb_variable5.mat');

%%
SNRi = 10*log10(fliplr(logspace(-1,1,numN)));

BW = rho_hat1*sqrt(log(4));
SNRo_Static = 10*log10(BW^2./meanStaticEM);
SNRo_Var1 = 10*log10(BW^2./meanVar1);
SNRo_Var2 = 10*log10(BW^2./meanVar2);
SNRo_Var3 = 10*log10(BW^2./meanVar3);
SNRoCOG = 10*log10(BW^2./meanCOG);

for ii = 1:numS
    figure; plot(SNRi,SNRo_Static(:,ii),SNRi,SNRo_Var1(:,ii),...
        SNRi,SNRo_Var2(:,ii),SNRi,SNRo_Var3(:,ii),SNRi,SNRoCOG(:,ii));
    xlabel('Input SNR (dB)');
    ylabel('Output SNR (dB)');
    legend('Static , \rho = 40',['Variable, \rho = ' num2str(rho_hat1)],...
        ['Variable, \rho = ' num2str(rho_hat2)],...
        ['Variable, \rho = ' num2str(rho_hat3)],...
        'COG','Location','southeast');
    title(['Performance of EM Estimators for \Lambda_s = ' num2str(Lams(ii))]);
end
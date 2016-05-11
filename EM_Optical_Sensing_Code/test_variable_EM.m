%% Test script for data generation
% Joshua Rapp
% Boston University
% EC 503

clear; close all; clc;

%% Static Data
Lr = 500; Lc = 500;
rho_hat = 40;

Lam_s = 100;
Lam_n = 1000;

[sig_pos,Sigma_cov,matDetect,listDetect,labels] = ...
    fcn_generate_correlated_data( Lr,Lc,rho_hat,Lam_s,Lam_n );
centroid = mean(listDetect);

% Apply EM
[x_var, Rvar] = variableEM(matDetect,listDetect,0,Lam_s,Lam_n);
x_var_est = x_var{end};
R_var_est = Rvar{end};

x_stat = staticEM(matDetect,listDetect,rho_hat,Lam_s,Lam_n);
x_stat_est = x_stat{end};

Err_var = sqrt(sum((sig_pos-x_var_est').^2));
Err_stat = sqrt(sum((sig_pos-x_stat_est').^2));

disp(['Variable Improvement: ' num2str(Err_stat-Err_var)]);

%% Plot

figure; imagesc(matDetect); axis ij image; colormap(gray);
hold on;
plot(sig_pos(2),sig_pos(1),'g+','MarkerSize',10,'LineWidth',3)
plot(x_var_est(2),x_var_est(1),'rx','MarkerSize',10,'LineWidth',3)
plot(x_stat_est(2),x_stat_est(1),'mp','MarkerSize',10,'LineWidth',3)
plot(centroid(2),centroid(1),'b*','MarkerSize',10,'LineWidth',3)
legend('Truth','Variable EM','Static EM','Centroid');
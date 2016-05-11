%% Test script for data generation
% Joshua Rapp
% Boston University
% EC 503

clear; close all; clc;

%% Static Data
Lr = 500; Lc = 500;
rho = 40;
Lam_s = 50;
Lam_n = 100;

[ sig_pos, matDetect,listDetect,label] = fcn_generate_distribution(Lr,Lc,rho,Lam_s,Lam_n);
centroid = mean(listDetect);

%% Apply EM
xhats = staticEM(matDetect,listDetect,rho,Lam_s,Lam_n,20);
xest = xhats{end};
figure; imagesc(matDetect); axis ij image; colormap(gray);
hold on;
plot(sig_pos(2),sig_pos(1),'g+','MarkerSize',10,'LineWidth',3);
plot(xest(2),xest(1),'rx','MarkerSize',10,'LineWidth',3);
plot(centroid(2),centroid(1),'b*','MarkerSize',10,'LineWidth',3);

legend('Truth','EM est.','Centroid');
%% 
rows = 1:Lr;
cols = 1:Lc;

% Theoretical Marginals
unif_rows = ones(Lr,1)/Lr;
marg_unif_rows = cumsum(unif_rows)/sum(unif_rows);
unif_cols = ones(Lc,1)/Lc;
marg_unif_cols = cumsum(unif_cols)/sum(unif_cols);

gauss_rows = (Lam_s*normpdf(rows,xest(1),rho)+Lam_n*unif_rows')/(Lam_s+Lam_n);
marg_gauss_rows = cumsum(gauss_rows)/sum(gauss_rows);
gauss_cols = (Lam_s*normpdf(cols,xest(2),rho)+Lam_n*unif_cols')/(Lam_s+Lam_n);
marg_gauss_cols = cumsum(gauss_cols)/sum(gauss_cols);

% Compute Marginals
data_rows = sum(matDetect,2);
marg_cdf_rows = cumsum(data_rows)/sum(data_rows);
figure; plot(rows,marg_cdf_rows,rows,marg_unif_rows,rows,marg_gauss_rows);
title('Empirical Marginal CDF (Rows)');
xlabel('Row'); 
legend('Empirical CDF','Noise CDF','Noise+Signal CDF','Location','northwest');

data_cols = sum(matDetect,1)';
marg_cdf_cols = cumsum(data_cols)/sum(data_cols);
figure; plot(cols,marg_cdf_cols,cols,marg_unif_cols,cols,marg_gauss_cols);
title('Empirical Marginal CDF (Columns)');
xlabel('Column'); 
legend('Empirical CDF','Noise CDF','Noise+Signal CDF','Location','northwest');
%%
MSE_Noise_cols = mean((marg_cdf_cols-marg_unif_cols).^2);
MSE_Noise_rows = mean((marg_cdf_rows-marg_unif_rows).^2);

MSE_Signal_cols = mean((marg_cdf_cols-marg_gauss_cols').^2);
MSE_Signal_rows = mean((marg_cdf_rows-marg_gauss_rows').^2);

[~,dist_predict] = min([MSE_Noise_cols+MSE_Noise_rows,MSE_Signal_cols+MSE_Signal_rows]);



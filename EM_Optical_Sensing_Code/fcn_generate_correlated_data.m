function [sig_pos,Sigma_cov,matDetect,listDetect,labels] = ...
    fcn_generate_correlated_data( Lr,Lc,sig_hat,Lam_s,Lam_n )
%FCN_GENERATE_DATA takes in parameters about signal and noise detection
%rates and generates a dataset based on the model of a circular Gaussian
%signal and uniform noise.
%
% The output includes both a 2D-detector view of detections (better for visualization)
% as well as a vector of detection coordinates (easier to process).
%
%*************************************************************************
% Input Parameters
%------------
% [Lr, Lc] = size of detector array represented by matrix
% rho = signal standard deviation
% sigma = prior standard deviation
% Lam_s = beam photo-conversion rate
% Lam_n Noise photoconversion rate
%
% Output Parameters
% ------------------
% sig_pos = coordinates of true signal position
% matDetect = matrix of signal detections
% listDetect = list of detection coordinates
%*************************************************************************

sig_pos = [sig_hat+round((Lr-2*sig_hat)*rand), sig_hat+round((Lc-2*sig_hat)*rand)];

numSig = poissrnd(Lam_s);
numNoise = poissrnd(Lam_n);

sig1 = randi([sig_hat/2,3*sig_hat/2]);
sig2 = randi([sig_hat/2,3*sig_hat/2]);
rho = -1+2*rand(1);
Sigma_cov = [sig1^2, rho*sig1*sig2; rho*sig1*sig2, sig2^2];
sigPreDetect = round(mvnrnd(sig_pos,Sigma_cov,numSig));

%sigPreDetect = [sig_pos(1)+round(sig_hat*randn(numSig,1)),sig_pos(2)+round(sig_hat*randn(numSig,1))];

sigDetect = sigPreDetect(sigPreDetect(:,1)>0,:);
sigDetect = sigDetect(sigDetect(:,2)>0,:);
sigDetect = sigDetect(sigDetect(:,1)<=Lr,:);
sigDetect = sigDetect(sigDetect(:,2)<=Lc,:);

noiseDetect = [randi(Lr,[numNoise, 1]),randi(Lc,[numNoise, 1])];
listDetect = [sigDetect; noiseDetect];
labels = [ones(length(sigDetect),1);zeros(length(noiseDetect),1)];

matDetect = zeros(Lr,Lc);

for ii= 1:length(listDetect)
    matDetect(listDetect(ii,1),listDetect(ii,2)) = matDetect(listDetect(ii,1),listDetect(ii,2))+1;
end


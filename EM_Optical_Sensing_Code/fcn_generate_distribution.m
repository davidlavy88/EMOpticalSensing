function [sig_pos,matDetect,listDetect,label] = fcn_generate_distribution(Lr,Lc,rho,Lam_s,Lam_n)
%FCN_GENERATE_DISTRIBUTION takes in parameters about signal and noise detection
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

numSig = poissrnd(Lam_s);
numNoise = poissrnd(Lam_n);

if randi([0 1])
    sig_pos = [rho+round((Lr-2*rho)*rand), rho+round((Lc-2*rho)*rand)];
    sigPreDetect = [sig_pos(1)+round(rho*randn(numSig,1)),sig_pos(2)+round(rho*randn(numSig,1))];

    sigDetect = sigPreDetect(sigPreDetect(:,1)>0,:);
    sigDetect = sigDetect(sigDetect(:,2)>0,:);
    sigDetect = sigDetect(sigDetect(:,1)<=Lr,:);
    sigDetect = sigDetect(sigDetect(:,2)<=Lc,:);

    noiseDetect = [randi(Lr,[numNoise, 1]),randi(Lc,[numNoise, 1])];
    listDetect = [sigDetect; noiseDetect];
    label = 1;
else
    sig_pos = [0,0];
    numDetect = numSig+numNoise;
    listDetect = [randi(Lr,[numDetect, 1]),randi(Lc,[numDetect, 1])];
    label = 0;
end

matDetect = zeros(Lr,Lc);

for ii= 1:length(listDetect)
    matDetect(listDetect(ii,1),listDetect(ii,2)) = matDetect(listDetect(ii,1),listDetect(ii,2))+1;
end


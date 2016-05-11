function [ sig_pos_all, matDetect, listDetect,labels  ] = fcn_generate_data( Lr,Lc,rho,Lam_s,Lam_n,num_beams  )
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
% Lam_n = Noise photoconversion rate
% num_beams = number of beams
%
% Output Parameters
% ------------------
% sig_pos = coordinates of true signal position (num_sig x 2 matrix)
% matDetect = matrix of signal detections
% listDetect = list of detection coordinates
%*************************************************************************

sig_pos_all = zeros(num_beams,2);
sigDetect_all = [];
for beam = 1:num_beams
    
    sig_pos = [rho+round((Lr-2*rho)*rand), rho+round((Lc-2*rho)*rand)];
    sig_pos_all(beam,:) = sig_pos;
    numSig = poissrnd(Lam_s);
    
    sigPreDetect = [sig_pos(1)+round(rho*randn(numSig,1)),sig_pos(2)+round(rho*randn(numSig,1))];

    sigDetect = sigPreDetect(sigPreDetect(:,1)>0,:);
    sigDetect = sigDetect(sigDetect(:,2)>0,:);
    sigDetect = sigDetect(sigDetect(:,1)<=Lr,:);
    sigDetect = sigDetect(sigDetect(:,2)<=Lc,:);
    
    sigDetect_all = [sigDetect_all ; sigDetect];

end

numNoise = poissrnd(Lam_n);


noiseDetect = [randi(Lr,[numNoise, 1]),randi(Lc,[numNoise, 1])];
listDetect = [unique(sigDetect_all,'rows'); noiseDetect];
labels = [ones(length(unique(sigDetect_all,'rows')),1);zeros(length(noiseDetect),1)];

matDetect = zeros(Lr,Lc);

for ii= 1:length(listDetect)
    matDetect(listDetect(ii,1),listDetect(ii,2)) = 1;
end
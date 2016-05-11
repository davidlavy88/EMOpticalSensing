function [xhats, Rhats] = variableEM(matDetect,listDetect,sigma_hat,Lam_S,Lam_N,itmax)
% use the expectation maximization algorithm to estimate the position on an
% optical beam on a 2D detector array
%
% Inputs:
%
% matDetect - an MxN matrix of photoevent counts at each grid point
%
% nonzero_coords - Lx2 matrix of coordinate pairs [row col; row col; row
% col; ...] where L is the number of locations where at least one photo count occured
%
% LS - mean # signal photoconversions
%
% LN - mean # noise photoevents
%
% itmax - maximum number of iterations allowed (default Inf)
%
% Outputs:
%
% xhat - the position estimate [row position column position]

if nargin == 5
    itmax = Inf;
end

numDetect = size(listDetect,1);
num_row = size(matDetect,1);
num_col = size(matDetect,2);
Area = num_row*num_col;

%% Initialization
% x estimate
cogsum = 0;
for nzpair = 1:numDetect
    row = listDetect(nzpair,1);
    col = listDetect(nzpair,2);
    rc_pair = [row col];
    counts = matDetect(row,col);
    cogsum = cogsum + counts*rc_pair;
end
xhat = cogsum / sum(matDetect(:));
xhat = round(xhat)';

% Covariance estimate
if sigma_hat==0
    Rhat = ((listDetect-ones(numDetect,1)*xhat')'*(listDetect-ones(numDetect,1)*xhat'))/numDetect;
else
    Rhat = diag([sigma_hat^2;sigma_hat^2]);
end
Rinv = inv(Rhat);

% apply EM
lambdaN = Lam_N/Area;
itnum = 1;
xhats{1} = xhat;
Rhats{1} = Rhat;
while true
    sigma_hat = sqrt(det(Rhat));

    num_sum = 0;
    denom_sum = 0;
    R_num_sum = 0;
    
    for nzpair = 1:numDetect
        
        row = listDetect(nzpair,1);
        col = listDetect(nzpair,2);            
        d = [row; col];
        lambdaS = Lam_S/(2*pi*sigma_hat) * exp(-0.5*(d - xhat)'*Rinv*(d - xhat));
        weight = lambdaS/(lambdaS + lambdaN);

        counts = matDetect(row,col);
        num_sum = num_sum + counts * weight * d;
        denom_sum = denom_sum + counts * weight;
        
        R_num_sum = R_num_sum+ weight*(d - xhat)*(d - xhat)';
        
    end
    
    xhat_tminus1 = xhat;
    xhat = round(num_sum/denom_sum);
    Rhat = R_num_sum/denom_sum;
    
    if isequal(xhat,xhat_tminus1) || itnum == itmax
        break;
    end
    itnum = itnum + 1;
    xhats{itnum,1} = xhat;
    Rhats{itnum,1} = Rhat;
end
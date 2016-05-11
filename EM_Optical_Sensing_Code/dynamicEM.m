function xhat = dynamicEM(detector_data,nonzero_coords,prev_xhat,sigma,rho,LS,LN)
% use the expectation maximization algorithm to estimate the position on an
% optical beam on a 2D detector array
%
% Inputs:
%
% detector_data - an MxN matrix of photoevent counts at each grid point
%
% prev_xhat - the most recent position estimate (2x1 or 1x2 vector)
%
% rho - beam spatial variance (scalar double) (probably in units of num
% grid points)
%
% LS - mean # signal photoconversions
%
% LN - mean $ noise photoevents
%
% Outputs:
%
% xhat - the position estimate [row position column position]

if size(prev_xhat) == size(ones(1,2))
    prev_xhat = prev_xhat';
end
L = size(nonzero_coords,1);
num_row = size(detector_data,1);
num_col = size(detector_data,2);
A = num_row*num_col;
Rinv = inv(diag([rho^2 rho^2]));

xhat = prev_xhat;

% determine sigma

static_xhats = staticEM(detector_data,nonzero_coords,rho,LS,LN);
staticX = static_xhats{end};

if isempty(sigma)
    sigma = 0.25*norm(staticX - prev_xhat);
end

% make prev_xhat the origin of the grid, get estimate and adjust using prev_xhat to get value with true
% origin
lambdaN = LN/A;
while true
    
    num_sum = 0;
    denom_sum = 0;
    xhat = xhat - prev_xhat;
    for nzpair = 1:L
            row = nonzero_coords(nzpair,1);
            col = nonzero_coords(nzpair,2);
            
            d = [row; col] - prev_xhat; % set previous xhat as origin
            lambdaS = LS/(2*pi*rho^2) * exp(-0.5*(d - xhat)'*Rinv*(d-xhat));
            w = lambdaS/(lambdaS + lambdaN);
            
            counts = detector_data(row,col);
            num_sum = num_sum + counts * w * d;
            denom_sum = denom_sum + counts * w;
    end
    
    xhat_tminus1 = xhat + prev_xhat;
    xhat = round(num_sum/(denom_sum + (rho/sigma)^2));
    xhat = xhat + prev_xhat; % reset origin
    
    if xhat == xhat_tminus1
%     if norm(xhat - xhat_tminus1) <= 0
        break;
    end
    
end
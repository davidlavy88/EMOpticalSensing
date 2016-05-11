function xhats = staticEM(detector_data,nonzero_coords,rho,LS,LN,itmax)
% use the expectation maximization algorithm to estimate the position on an
% optical beam on a 2D detector array
%
% Inputs:
%
% detector_data - an MxN matrix of photoevent counts at each grid point
%
% nonzero_coords - Lx2 matrix of coordinate pairs [row col; row col; row
% col; ...] where L is the number of locations where at least one photo count occured
%
% rho - beam spatial variance (scalar double) (probably in units of num
% grid points)
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

L = size(nonzero_coords,1);
num_row = size(detector_data,1);
num_col = size(detector_data,2);
A = num_row*num_col;
R = diag([rho^2 rho^2]);
Rinv = inv(R);

% get initial estimate (COG estimate)
cogsum = 0;
for nzpair = 1:L
    row = nonzero_coords(nzpair,1);
    col = nonzero_coords(nzpair,2);
    rc_pair = [row col];
    counts = detector_data(row,col);
    cogsum = cogsum + counts*rc_pair;
end
xhat = cogsum / sum(detector_data(:));
xhat = round(xhat)';

% apply EM
lambdaN = LN/A;
itnum = 1;
xhats{1} = xhat;
while true
    
    num_sum = 0;
    denom_sum = 0;
    
    for nzpair = 1:L
        
        row = nonzero_coords(nzpair,1);
        col = nonzero_coords(nzpair,2);            
        d = [row; col];
        lambdaS = LS/(2*pi*rho^2) * exp(-0.5*(d - xhat)'*Rinv*(d - xhat));
        w = lambdaS/(lambdaS + lambdaN);

        counts = detector_data(row,col);
        num_sum = num_sum + counts * w * d;
        denom_sum = denom_sum + counts * w;
        
    end
    
    xhat_tminus1 = xhat;
    xhat = round(num_sum/denom_sum);
    
    if isequal(xhat,xhat_tminus1) || itnum == itmax
        break;
    end
    itnum = itnum + 1;
    xhats{itnum} = xhat;
end
function xhat = multibeamEM(detector_data,nonzero_coords,rho,LS,LN,...
    prev_xhats,sigma,num_beams)
% use the expectation maximization algorithm to estimate the position on an
% optical beam on a 2D detector array
%
% need to  assign cluster to prev_xhat entry based on centroid and
% prev_xhat position
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
% prev_xhats - most recent position estimates (cell array, one cell per beam) (use [] for static case)
%
% sigma - standard deviation for beam locationn(use[] for static case)
%
% num_beams - number of beams
%
% Outputs:
%
% xhats - the position estimates for each beam cell array of cell arrays

[clusters,Cs] = kmeans(nonzero_coords,num_beams,'Replicates',20,'Distance','cityblock');

% assign clusters to prev_xhats
cluster_prev_xhats = cell(1,length(prev_xhats)); % idx 1 corresponds to cluster 1, idx 2 to cluster 2, etc.
xhat_choices = prev_xhats;
if ~isempty(prev_xhats)
    for cluster_num = 1:size(Cs,1)
        cent = Cs(cluster_num,:)'; % centroid of cluster
        xhatDiffs = cellfun(@(x) norm(x-cent),xhat_choices); % find distance for each xhat to centroid
        [~,idx] = min(xhatDiffs); % index of the unassigned xhat that is closest to centroid
        cluster_prev_xhats{cluster_num} = xhat_choices{idx}; % assign nearest unassigned xhat to cluster
        xhat_choices(idx) = []; % remove assigned xhat
    end
end
    

xhat = cell(1,num_beams);
for beam_num = 1:num_beams
    beam_data = nonzero_coords(clusters == beam_num,:);
    if isempty(prev_xhats)
        xhats = staticEM(detector_data,beam_data,rho,LS,LN);
        xhat{beam_num} = xhats{end};
    elseif ~isempty(prev_xhats)
        prev_xhat = cluster_prev_xhats{beam_num};
        xhat{beam_num} = dynamicEM(detector_data,beam_data,prev_xhat,sigma,rho,LS,LN);
    end
end   
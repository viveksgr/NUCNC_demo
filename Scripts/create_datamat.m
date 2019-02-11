function [X,y] = create_datamat(behav_data,fir_list,edge_,history)
% Create datamatrix to be used for predictions in the GLM
% X = NxP, N = data points, P = spatial predictors
% Discretization of space done at edges - not accurate. 
% y = labels to be used in the GLM, read out from fir_list in the order provided by behav_data
% Toggle history to involve history terms - toggle LNP vs GLM
% -------------------------------------------------------------------------
% Vivek Sagar, VivekSagar2016@u.northwestern.edu
% Feb 10, 2019

if nargin<4
    history = false; % Append spike history?
end
 
X = zeros(length(behav_data),length(edge_{1})*length(edge_{2}));
y = zeros(length(behav_data),1);
for tt = 1:length(behav_data)
    xbin = find(edge_{1}>behav_data(tt,1),1)-1;
    ybin = find(edge_{2}>behav_data(tt,2),1)-1;
    reduced_pos = sub2ind([length(edge_{1}),length(edge_{2})],xbin,ybin);
    X(tt,reduced_pos)=1;
    y(tt) = fir_list(xbin,ybin);
end
X = sparse(X);

if history
    X_t1 = circshift(X,1);
    X = [X X_t1]; % Append history.
    X(1,:) = []; % Because there is no history for the first term.
    y(1,:) = [];
end




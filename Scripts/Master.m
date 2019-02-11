% LNP encoding and MAP decoding of grid cell data.
% Specify the path of grid cell data in <datapath>
% units_ = Nx1 cell array of last 4 IDs of units or {} for all recorded units.
% Make sure t*c*.mat names are non-degenerate.
% -------------------------------------------------------------------------
% Vivek Sagar, VivekSagar2016@u.northwestern.edu
% Feb 10, 2019

datapath = 'C:\NUCNC\Data';
units_ = {'t5c1','t8c1'}; % Cell array of last 4 characters in the filename of unit data. Put {} for all units.
display_ = 'on'; % Display results
size_maze = [-50 50;-50 50]; % Size of the rectangular maze (1,2) is the right edge and (2,1) is the bottom edge.
nbins = 10;

%% Visualize Grid Fields
% Locate data
try
    if ~isempty(units_) % Not all units are used in analysis.
        cell_names = dir(fullfile(datapath,sprintf('*%s.mat',units_{1}))); % First unit
        if length(units_)>1 % Got more units? No problem!
            for cell_ind = 2:length(units_)
                cell_names(cell_ind) = dir(fullfile(datapath,sprintf('*%s.mat',units_{cell_ind}))); % Struct of datafilescell_{1}
            end
        end
    else % Use all units. Go big!
        cell_names = dir(fullfile(datapath,sprintf('*t*c*.mat'))); % Struct of datafiles
    end
    assert(~isempty(cell_names)) % Double check that the data is actually there.
catch
    error('Wrong datapath... way to start things!')
end

res = struct();
% Calculate mean firing rate
for ii = 1:length(cell_names) % Loop over all units
    unit_data = load(fullfile(datapath,cell_names(ii).name)); % Load unit data
    [res(ii).behav_data, res(ii).spike_data, res(ii).time] = extract_griddata(unit_data); % Extract useful variables
    [res(ii).fir_freq, edge_] = count_griddata(res(ii).behav_data, res(ii).spike_data, res(ii).time,size_maze,nbins); % Firing rate is unblurred. Blur it for extra beauty.
    
    % Plot grid fields
    figure('visible',display_)
    imagesc(edge_{1},edge_{2},res(ii).fir_freq)
    axis tight
    str = sprintf('Spatial firing frequency for unit %s (Hz)',units_{ii});
    title(str);
    xlabel('X-position (cm)')
    ylabel('Y-position (cm)')
    colormap(hot);
    colorbar;
    print(fullfile(datapath,units_{ii}), '-dpng');
end
fprintf('Grid fields saved in datapath\n')

%% LNP/GLM model - Encoding
% This code can definitely be made faster - but hey, don't judge a hastily
% written code for the purpose of live demonstration for its efficiency. 

% The model is overparametrized - maybe we don't need so many basis functions
% for defining a spatial RF?

% This is not truly out-of-sample prediction. 
% Sampling of the stimulus was quite good - the benefit of doing all this will not apparant here.

nfolds = 10; % Number of folds for cross validation. 
ii = 1; % Neuron number. Loop over ii for all neurons.
fprintf('Unit %d/%d\n',ii,length(cell_names))

% LNP
[X,y] = create_datamat(res(ii).behav_data,res(ii).fir_freq,edge_);
fold_ind = crossvalind('Kfold',length(y),nfolds); % Just running 1 fold, because GLM is slow af.
test_mask = fold_ind==2;
X_train = X(~test_mask,:);
X_test = X(test_mask,:);
y_train = y(~test_mask);
y_test = y(test_mask);

% Put prior for regularization
% p_prior = double(y_train>0);

p = glmfit(X_train,y_train,'poisson');%,'weights',p_prior); % Link function is log by default.
y_pred = glmval(p,X_test,'log');

spike_pred = poissrnd(y_pred);
pred_acc = corrcoef(spike_pred,res(ii).spike_data(test_mask));
fprintf('Prediction accuracy for the LNP model: %f\n',pred_acc(2));

% Plot Receptive Field
p(1) = []; %Pop out constant term
p_mask = p<-99;
p(p_mask) = min(p(~p_mask));
p_sq = reshape(p,floor(sqrt(length(p))),[]);
figure('visible',display_)
imagesc(edge_{1},edge_{2},p_sq)
axis tight
str = sprintf('LNP firing frequency for unit %s (Hz)',units_{ii});
title(str);
xlabel('X-position (cm)')
ylabel('Y-position (cm)')
colormap(hot);
colorbar;
print(fullfile(datapath,sprintf('LNP%s',units_{ii})), '-dpng');

% % GLM with history
% [X2,~] = create_datamat(res(ii).behav_data,res(ii).fir_freq,edge_,true);
% X2_train = X2(~test_mask,:);
% X2_test = X2(test_mask,:);
% p = glmfit(X2_train,y_train,'poisson'); % Link function is exponential by default.
% y_pred2 = glmval(p,X2_test,'log');
% pred_acc2 = corrcoef(y_pred2,y_test);
% sprintf('Prediction accuracy for the GLM model: %f',pred_acc2(2));

%% MAP estimate - decoding

% Easy way - estimate the entire encoding table 
% Should only be used for an ensemble. Grid cells do not identify spatial
% locations uniquely. 
% Compromised perfomance for ease of understanding the code. 
% Just the sketch of code - add all loops (n_neurons, n_stim, n_test bins)
% 
% % Ideally:
% tuning_curves = zeros(length(cell_names),(nbins+1)^2);
% for ii = 1:length(cell_names)
%     tuning_curves(ii,:)=res(ii).fir_freq(:);
% end

ii = 1;
tuning_curves = res(ii).fir_freq(:);
% Test spike_data = y_test;
dt = median(diff(res(ii).time));
% Loop over spike_data bins
testbin = 10; % Choose the firing rate
stim_bin = 12; % Choose the spatial location
fir_rate = y_test(testbin);
tuning_location = tuning_curves(stim_bin);

% What is the (unnormalized) probability of being at this tuning_location given fir_rate?
% Assume uniform prior.
poster_ = (tuning_location.^fir_rate)*exp(-1*dt*tuning_location); % Check the position of dt
% This gives posterior probability at location given by stim_bin. Loop over
% stim_bin to obtain Xmap for the firing rate given by testbin. 

% Right way - use the GLM to obtain x_map straight away.
% Maximize log p(X|D,k) which is "trivial" and is left as an exercise for
% the reader.
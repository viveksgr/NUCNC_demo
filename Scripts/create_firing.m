function firing_rates = create_firing(behav_data,fir_list)
% Create labels to be used in the GLM
% Reads out firing rates from fir_list in the order provided by behav_data
% -------------------------------------------------------------------------
% Vivek Sagar, VivekSagar2016@u.northwestern.edu
% Feb 10, 2019


firing_rates = zeros(length(behav_data),1);
for ii = 1:length(behav_data);
    
function [fir_freq, edge_]= count_griddata(behav_data, spike_data, time, size_maze,bin)
% Firing frequency maps of grid cells. 
% behav_data = Nx2 positions
% spike_data = Nx1 binary spikes
% time = Nx1 time stamps
% size_maze = 2x2 dimension of the size of maze. Row 1 for x-size and 2 for y-size.  
% This code doesn't measure the time spent in each bin accurately.
% Also hist3 is weird - last rows and columns are kept nans to maintain size. 
% -------------------------------------------------------------------------
% Vivek Sagar, VivekSagar2016@u.northwestern.edu
% Feb 10, 2019

if nargin <5
    bin = 25;
end

scalex = (size_maze(1,2)-size_maze(1,1))/bin;
scaley = (size_maze(2,2)-size_maze(2,1))/bin;
edge_ = {size_maze(1,1):scalex:size_maze(1,2)
         size_maze(2,1):scaley:size_maze(2,2)};

dt = median(diff(time)); % Approximation - ideally should use the whole time series.
firing_pos = behav_data(logical(spike_data),:); % Positions where the unit spikes.
firing_pos_count = hist3(firing_pos,'Edges', edge_); % Number of spikes in each spatial bin
time_pos = hist3(behav_data,'Edges', edge_).*dt; % Time in each spatial bin
fir_freq = firing_pos_count./time_pos;
end
function [behav_data, spike_data, time] = extract_griddata(unit_data,speed_cut)
% Unit_data is any of the '*t*c*.mat' files loaded in struct form.
% behav_data =  Nx2 behav_data (x-pos, y-pos)
% spike_data =  NX1 binary spike data, N are matched time points
% time = NX1 time points
% Toggle speed_cut, remove data points when the animal was stationary. For
% removing SWR effects.
% -------------------------------------------------------------------------
% Vivek Sagar, VivekSagar2016@u.northwestern.edu
% Feb 10, 2019

if nargin<2
    speed_cut = true;
end

pos = [(unit_data.x1+unit_data.x2)/2 (unit_data.y1+unit_data.y2)/2];
time = union(unit_data.t, unit_data.ts);
behav_data = interp1(unit_data.t, pos, time);
spike_data = double(ismember(time,unit_data.ts));

% Remove NaNs
rem_ind = or(isnan(behav_data(:,1)),isnan(behav_data(:,2)));
behav_data(rem_ind,:) = [];
time(rem_ind) = [];
spike_data(rem_ind) = [];

% Remove data points when the animal was stationary or seemed to move with
% weird speed.
if speed_cut
    speed = sqrt(sum((diff(behav_data)./diff(time)).^2,2));
    speed_id = [or((speed<5),(speed>100))];
    behav_data(speed_id,:) = [];
    time(speed_id) = [];
    spike_data(speed_id) = [];
end
end

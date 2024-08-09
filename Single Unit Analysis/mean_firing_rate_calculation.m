function [binning_firing_rate, mean_firing_rate] = mean_firing_rate_calculation(clu_spike_times)
% calculates mean firing rate in bins for the cluster

time_bins = 100; % in seconds

file_length = clu_spike_times(end) - clu_spike_times(1);
mean_firing_rate = length(clu_spike_times)/ file_length; 
[n, edges] = histcounts(clu_spike_times, 0:time_bins:file_length);
binning_firing_rate = n/time_bins; 

end
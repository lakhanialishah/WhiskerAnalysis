function missed_spikes = missed_spikes_based_on_amplitude(clu_amplitudes)
% calculates the percent of missed spikes based on gaussian fit to amplitudes data
% first input is vector of amplitudes for one cluster

f = figure('visible', 'off');
hold on
h = histogram(clu_amplitudes);
fit = histfit(clu_amplitudes);

detection_threshold = 10.6; %what I see in phy 
missed_spikes = find(fit(2).XData > detection_threshold, 1);

end
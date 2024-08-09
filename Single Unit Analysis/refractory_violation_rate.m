function violation_rate = refractory_violation_rate(time_bin, clu_spike_times)
% calculates the refractory violation rate in percentage of the unit 
% first input is the bin (ex: 1 ms, 2 ms) 
% second input is the spiketimes of that cluster

file_length = clu_spike_times(end) - clu_spike_times(1);
isi = clu_spike_times(2:end) - clu_spike_times(1:end-1);
[n, edges] = histcounts(isi, 0:.001:file_length); %n is the number of spikes in each time bin of 1 ms
violation_rate = n(1:time_bin)/length(clu_spike_times) * 100; 

end
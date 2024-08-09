function spont_fr = spont_firing_rate(allstim_matrix, spike_times)
% calculates spontaneous firing rate based on 1 second before all stims for the recording day

prewindow = [-1 0];
verticalconcat = [];
for stimNum = 1:size(allstim_matrix,1)
    spikes = spike_times(spike_times > allstim_matrix(stimNum,3)+prewindow(1) & spike_times < allstim_matrix(stimNum,3)+prewindow(2));
    verticalconcat = vertcat(verticalconcat, length(spikes));
end
spont_fr = sum(verticalconcat)/size(allstim_matrix,1); %in Hz because looking at a one second window

end




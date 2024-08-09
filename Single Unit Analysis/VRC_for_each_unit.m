function step2_VRC = VRC_for_each_unit(step2_matrix, spike_times)
% to calculate VRC for each unit based on Step2 only (this is the recording
% based on what whisker I heard might be the BW and only that whisker is
% stimulated) 

window = [0 .05]; % in seconds so this is 50 ms
prewindow = [-1 0]; % window pre stim

stim_spike_matrix = [];
if isempty(step2_matrix) == 0
    for stimNum = 1: length(step2_matrix)
        spikes = spike_times(spike_times > step2_matrix(stimNum,3)+window(1) & spike_times < step2_matrix(stimNum,3)+window(2));
        prewindowspikes = spike_times(spike_times > step2_matrix(stimNum,3)+prewindow(1) & spike_times < step2_matrix(stimNum,3)+prewindow(2));
        stim_spike_matrix = vertcat(stim_spike_matrix, [step2_matrix(stimNum, 2) (length(spikes)/.05) length(prewindowspikes) (length(spikes)/.05)-(length(prewindowspikes))]); %velocity, postwindow spike rate, prewindow spike rate, evoked spike rate
    end
    
    %collapsing all 50 trials of each velocity
    unique_velocities = unique(stim_spike_matrix(:,1));
    for vel = 1:length(unique_velocities)
        indices = find(stim_spike_matrix(:,1) == unique_velocities(vel));
        mean_for_vel = mean(stim_spike_matrix(indices, 4));
        step2_VRC(vel,1) = unique_velocities(vel);
        step2_VRC(vel,2) = mean_for_vel;
    end
else 
    step2_VRC = [];
end
end
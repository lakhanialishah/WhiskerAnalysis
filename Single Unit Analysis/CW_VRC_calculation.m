function CW_VRC = CW_VRC_calculation(step3_matrix, spike_times, CW_dac)
%getting the VRC for the CW for the unit despite what the BW is
%only used when for when the probe definitely ended up in a barrel so when it did not, the CW_dac is 0

window = [0 .05]; % in seconds so this is 50 ms
prewindow = [-.05 0]; % window pre stim

if isempty(step3_matrix) == 0 
    vels = unique(step3_matrix(:,2));
%     CW_VRC(:,1) = vels;
    for v = 1: length(vels)
        indices = find(step3_matrix(:,1) == CW_dac & step3_matrix(:,2) == vels(v));
        spikefr_in_vels = [];
        for stimnum = 1:length(indices)
            spontspikes = spike_times(spike_times > step3_matrix(indices(stimnum),3)+prewindow(1) & spike_times < step3_matrix(indices(stimnum),3)+prewindow(2));
            evokedspikes = spike_times(spike_times > step3_matrix(indices(stimnum),3)+window(1) & spike_times < step3_matrix(indices(stimnum),3)+window(2));
            true_evoked_spikes_fr = (length(evokedspikes) - length(spontspikes))/ .05;
            spikefr_in_vels = vertcat(true_evoked_spikes_fr, spikefr_in_vels);
        end
        CW_VRC(v,1) = mean(spikefr_in_vels);
    end
else
    CW_VRC = [];
end
end
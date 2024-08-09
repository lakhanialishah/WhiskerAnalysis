function threshold = vel_threshold_for_unit(step3_matrix, spike_times, interested_dac)
% to determine the velocity threshold for each unit, using BW stims

window = [0 .05]; % in seconds so this is 50 ms
prewindow = [-1 0]; % window pre stim
velocities = [0, 65, 195, 326, 456, 587, 797];

if isempty(step3_matrix) == 0
    
    %calculate spont fr for all trials
    all_spont_fr = [];
    for trial = 1:length(step3_matrix)
        prewindowspikes = spike_times(spike_times > step3_matrix(trial,3)+prewindow(1) & spike_times < step3_matrix(trial,3)+prewindow(2));
        spont_fr_trial = length(prewindowspikes); 
        all_spont_fr = vertcat(all_spont_fr, spont_fr_trial);
    end
    
    %calculate evoked fr at each velocity
    vels = unique(step3_matrix(:,2));
    matrix_pvalues = [];
    for v = 1:length(vels)
        holding_matrix = [];
        for dac = interested_dac %only looking at the BW for the threshold
            indices = find(step3_matrix(:, 1) == dac & step3_matrix(:,2) == vels(v));
            for stimNum = 1: length(indices)
                spikes = spike_times(spike_times > step3_matrix(indices(stimNum),3)+window(1) & spike_times < step3_matrix(indices(stimNum),3)+window(2));
                holding_matrix = vertcat(holding_matrix, length(spikes)/.05); %saving evoked fr
            end
            matrix = nan(1575,2);
            matrix(:,1) = all_spont_fr;
            matrix(1:25, 2) = holding_matrix;
            p = ranksum(matrix(:,1), matrix(:,2), 'tail', 'left');
            matrix_pvalues(v,1) = p;
            clear indices h p matrix
        end
    end
    index = find(matrix_pvalues < (.05/7), 1); %7 velocities
    threshold = velocities(index);
else
    threshold = NaN;
end
end
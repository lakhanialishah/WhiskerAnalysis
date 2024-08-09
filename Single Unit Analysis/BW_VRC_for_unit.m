function [BW_VRC_unit, BW_dac_unit] = BW_VRC_for_unit(step3_matrix, spike_times, BW_day_dac)
% finding BW for each unit instead of BW for the probe's penetration
% taking the data from step3 so that I have all the whiskers to compare
% if two whiskers elicit the same response in the unit, then the code uses either

% wilcoxon rank sum test to see the whisker is responsive
% doing .05/9 to take into account all 9 whiskers

% to test one unit:
% spike_times = IDlist{1013,2}.clu_info.spikeTimes;
% step3_matrix = IDlist{1013,2}.clu_info.step3_matrix;

window = [0 .05]; % in seconds so this is 50 ms
prewindow = [-1 0]; % window pre stim

if isempty(step3_matrix) == 0
    
    %calculate spont fr for all trials
    all_spont_fr = [];
    for trial = 1:length(step3_matrix)
        prewindowspikes = spike_times(spike_times > step3_matrix(trial,3)+prewindow(1) & spike_times < step3_matrix(trial,3)+prewindow(2));
        spont_fr_trial = length(prewindowspikes); %if I change prewindow to .05 then make sure to divide 
        all_spont_fr = vertcat(all_spont_fr, spont_fr_trial);
    end
    
    %calculate evoked fr at highest velocity 
    vels = unique(step3_matrix(:,2));
    for v = 7 % highest velocity only 
        dac_evoked_fr_matrix = [];
        for dac = 1:9
            indices = find(step3_matrix(:, 1) == dac & step3_matrix(:,2) == vels(v));
            dac_holding_matrix = [];
            for stimNum = 1: length(indices)
                spikes = spike_times(spike_times > step3_matrix(indices(stimNum),3)+window(1) & spike_times < step3_matrix(indices(stimNum),3)+window(2));
                dac_holding_matrix = vertcat(dac_holding_matrix, length(spikes)/.05); %saving evoked fr
            end
            matrix = nan(1575,2);
            matrix(:,1) = all_spont_fr;
            matrix(1:25, 2) = dac_holding_matrix;
            p = ranksum(matrix(:,1), matrix(:,2), 'tail', 'left'); %kruskalwallis(matrix, [], 'off');
            dac_evoked_fr_matrix(dac,1) = dac;
            dac_evoked_fr_matrix(dac,2) = mean(dac_holding_matrix); %evoked
            dac_evoked_fr_matrix(dac,3) = mean(all_spont_fr); %spont
            dac_evoked_fr_matrix(dac,4) = p; %pvalue 
            clear indices h p
        end
    end

    %to easily tell if whisker reached signficance
    for i = 1: size(dac_evoked_fr_matrix,1)
        if dac_evoked_fr_matrix(i,4) < 0.0056 % for the 9 whiskers 
            dac_evoked_fr_matrix(i,5) = 1; %1 if it reaches significance
        else
            dac_evoked_fr_matrix(i,5) = 0;
        end
    end
    
    %finding the BW for the unit
    if sum(dac_evoked_fr_matrix(:,5)) > 0
        BW_dac = find(max(dac_evoked_fr_matrix(:,2)) == dac_evoked_fr_matrix(:,2));
    else
        BW_dac = NaN;
    end
    
    %to make sure it's only getting one whisker, taking the one that was the BW for the day
    if length(BW_dac) > 1
        if sum(BW_dac == BW_day_dac) == 1 %if one of these dacs was also the BW for that day
            index = find(BW_dac == BW_day_dac);
            BW_dac = BW_dac(index);
        else
            BW_dac = BW_dac(1);
        end
    end
    
    %for the output of the function
    if BW_dac > 0
        BW_dac_unit = BW_dac;
    else
        BW_dac_unit = NaN;
    end
    
    %calculate VRC for BW of that unit based on step3
    BW_VRC_unit(:,1) = vels;
    for v = 1: length(vels)
        indices = find(step3_matrix(:,1) == BW_dac & step3_matrix(:,2) == vels(v));
        spikefr_in_vels = [];
        for stimnum = 1:length(indices)
            spontspikes = spike_times(spike_times > step3_matrix(indices(stimnum),3)+prewindow(1) & spike_times < step3_matrix(indices(stimnum),3)+prewindow(2));
            evokedspikes = spike_times(spike_times > step3_matrix(indices(stimnum),3)+window(1) & spike_times < step3_matrix(indices(stimnum),3)+window(2));
            true_evoked_spikerate = (length(evokedspikes)/.05) - length(spontspikes);
            spikefr_in_vels = vertcat(true_evoked_spikerate, spikefr_in_vels);
        end
        BW_VRC_unit(v,2) = mean(spikefr_in_vels);
    end
else
    BW_VRC_unit = [];
    BW_dac_unit = NaN;
end
end


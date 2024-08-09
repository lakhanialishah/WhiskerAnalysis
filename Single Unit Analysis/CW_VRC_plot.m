function [groups, IDlist] = CW_VRC_plot(groups, IDlist, fields_in_groups)
% plotting CW evoked firing rate in different ways from step3
% only when the probe is in a barrel

% individual units displayed on each subplot, responders and non-responders
f = figure;
title('Responders and nonresponders')
hold on
units_in_groups_all = [];
mice_in_groups = [];
for g = 1:length(fields_in_groups)
    subplot(3,ceil(length(fields_in_groups)/3),g) %change the dimensions based on how many groups I have
    hold on
    subtract_unit = 0;
    total_units = [];
    mice_days = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        if isnan(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC) == 0
            plot(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit(:,1), IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC)
            mice_days = vertcat(mice_days, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.expday);
            total_units = vertcat(total_units, unit);
        end
    end
    plot(groups.(fields_in_groups{g}).BW_VRC_of_unit(:,1), groups.(fields_in_groups{g}).CW_VRC(:,1), 'k', 'LineWidth', 2)
    actual_mice = length(unique(mice_days));
    actual_units = length(total_units);
    text(.02,.95, [num2str(actual_units) 'units' '\newline' num2str(actual_mice) 'mice'],'Units','normalized')
    title(fields_in_groups(g))
    units_in_groups_all = vertcat(units_in_groups_all, actual_units);
    mice_in_groups = vertcat(mice_in_groups, actual_mice);
end
linkaxes()
xlim([0 .6])
han = axes(f,'visible','off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han,'Firing Rates (Hz)');
title(han,'Evoked Firing Rates (All Units)');

% individual units displayed on each subplot, responders only
f1 = figure;
% title('Responders only')
hold on
units_in_groups_resp_only = [];
mice_in_groups = [];
window = [0 .05]; % in seconds so this is 50 ms
prewindow = [-1 0]; % window pre stim
group_VRCs_responders = [];
group_VRCs_responders_std = [];
for g = 1:length(fields_in_groups)
    subplot(3,ceil(length(fields_in_groups)/3),g) %change the dimensions based on how many groups I have
    hold on
    %     subtract_unit = 0;
    mice_days = [];
    compile_VRCs = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        step3_matrix = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.step3_matrix;
        if isempty(step3_matrix) == 1
            continue
        end
        spike_times = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spikeTimes;
        vels = unique(step3_matrix(:,2));
        all_spont_fr = [];
        for trial = 1:length(step3_matrix)
            prewindowspikes = spike_times(spike_times > step3_matrix(trial,3)+prewindow(1) & spike_times < step3_matrix(trial,3)+prewindow(2));
            spont_fr_trial = length(prewindowspikes); %if I change prewindow to .05 then make sure to divide
            all_spont_fr = vertcat(all_spont_fr, spont_fr_trial);
        end
        for v = 7 % which velocity I'm looking at 
            for dac = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_dac %only interested in the CW data
                indices = find(step3_matrix(:, 1) == dac & step3_matrix(:,2) == vels(v));
                dac_holding_matrix = [];
                for stimNum = 1: length(indices)
                    spikes = spike_times(spike_times > step3_matrix(indices(stimNum),3)+window(1) & spike_times < step3_matrix(indices(stimNum),3)+window(2));
                    dac_holding_matrix = vertcat(dac_holding_matrix, length(spikes)/.05); %saving evoked fr
                end
                if isempty(dac_holding_matrix) == 1
                    continue
                end
                matrix = nan(1575,2);
                matrix(:,1) = all_spont_fr;
                matrix(1:25, 2) = dac_holding_matrix;
                p = ranksum(matrix(:,1), matrix(:,2), 'tail', 'left');
                IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.isCWresp = NaN;
                if p < .0056 %for the 9 whiskers
                    plot(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit(:,1), IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC)
                    mice_days = vertcat(mice_days, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.expday);
                    compile_VRCs = horzcat(compile_VRCs, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC);
                    IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.isCWresp = 1; %so I can find if it was CW responsive later
                end
            end
        end
        clear matrix dac_holding_matrix
    end
    groups.(fields_in_groups{g}).CW_VRC_responders_table = compile_VRCs;
    title(fields_in_groups(g))
    if isempty(compile_VRCs) == 1
        actual_units = 0;
        units_in_groups_resp_only = vertcat(units_in_groups_resp_only, actual_units);
        avg_group_VRC = zeros(1,7)' - 5;
        std_group_VRC = zeros(1,7)' - 1;
        group_VRCs_responders = horzcat(group_VRCs_responders, avg_group_VRC);
        group_VRCs_responders_std = horzcat(group_VRCs_responders_std, std_group_VRC);
        continue 
    end
    avg_group_VRC = mean(compile_VRCs, 2);
    std_group_VRC = std(compile_VRCs, 0, 2)/(sqrt(size(compile_VRCs,2)));
    group_VRCs_responders = horzcat(group_VRCs_responders, avg_group_VRC);
    group_VRCs_responders_std = horzcat(group_VRCs_responders_std, std_group_VRC);
    if isempty(avg_group_VRC) == 0
        plot(groups.(fields_in_groups{g}).BW_VRC_of_unit(:,1), mean(compile_VRCs, 2), 'k', 'LineWidth', 2)
    else
        continue
    end
    actual_mice = length(unique(mice_days));
    actual_units = size(compile_VRCs,2);
    text(.02,.95, [num2str(actual_units) 'units' '\newline' num2str(actual_mice) 'mice'],'Units','normalized')
    units_in_groups_resp_only = vertcat(units_in_groups_resp_only, actual_units);
    mice_in_groups = vertcat(mice_in_groups, actual_mice);
end
linkaxes()
xlim([0, 0.6])
ylim([0 27.5])
han2 = axes(f1,'visible','off');
han2.Title.Visible = 'on';
han2.XLabel.Visible = 'on';
han2.YLabel.Visible = 'on';
ylabel(han2,'Firing Rates (Hz)');
title(han2,'Evoked Firing Rates (Responders Only)');

% plotting average of CW VRC responders only
figure;
hold on
colors = turbo(size(group_VRCs_responders,2));
for g = 1:size(group_VRCs_responders,2)
    if isempty(group_VRCs_responders(:,g)) == 1
        continue
    end
    e = errorbar([0, 65, 195, 326, 456, 587, 797], group_VRCs_responders(:,g), group_VRCs_responders_std(:,g), 'LineWidth', 2, 'Color', colors(g,:));
    text(.78,g*.06, [num2str(units_in_groups_resp_only(g)) 'units ' num2str(mice_in_groups(g)) 'mice'],'Units','normalized', 'Color', colors(g,:))
end
title('Average CW VRC (responders only)')
ylabel('CW VRC')
xlabel('Velocity')
xlim([0 800])
ylim([0 27.5])
legend(fields_in_groups, 'Location', 'northwest', 'Box', 'off')


% proportion of CW units that are responders 
figure;
hold on 
% colors = turbo(length(fields_in_groups));
for g = 1:length(fields_in_groups)
    bar(g, units_in_groups_resp_only(g,1)/units_in_groups_all(g,1), 'FaceColor', colors(g,:))
    text(g-.15, .05, (num2str(units_in_groups_all(g,1))))
    groups.(fields_in_groups{g}).CW_proportion_resp = units_in_groups_resp_only(g,1)/units_in_groups_all(g,1); %so I can put it in an excel sheet/prism
end
ylim([0 1]);
xlabel('Experimental Group')
ylabel('Proportion of CW units that respond significantly')
xticks(1:length(fieldnames(groups)))
xticklabels(fields_in_groups)


% all mean VRCs on the same subplot, responders and non-responders
figure
hold on
% colors = turbo(size(group_VRCs_responders,2));
% ko = [0.9290, 0.6940, 0.1250];
% wt = [0, 0.4470, 0.7410];
% wd = [0.9, 0.0580, 0.1];
% colors = [ko ; wd];
for g = 1:size(group_VRCs_responders,2)
%     e = errorbar(groups.(fields_in_groups{g}).BW_VRC_of_unit(:,1), groups.(fields_in_groups{g}).CW_VRC(:,1),groups.(fields_in_groups{g}).CW_VRC(:,2), 'LineWidth', 2); %, 'Color', colors(g,:))
    e = errorbar([0, 65, 195, 326, 456, 587, 797], groups.(fields_in_groups{g}).CW_VRC(:,1), groups.(fields_in_groups{g}).CW_VRC(:,2), 'LineWidth', 2); %, 'Color', colors(g,:))
    e.Color = colors(g,:);
    text(.78,g*.06, [num2str(units_in_groups_all(g)) 'units ' num2str(mice_in_groups(g)) 'mice'],'Units','normalized', 'Color', colors(g,:))
end
legend(fields_in_groups, 'Location', 'northwest', 'Box', 'off')
xlabel('Velocity')
ylabel('Firing Rate (Hz)')
ylim([-.5 27.5])
xlim([0 800])
title('Average CW VRC (responders and nonresponders only)')

end
function evoked_fr_plots_BW_of_unit(groups, IDlist, fields_in_groups)
% plotting evoked firing rate in different ways

% individual units displayed on each subplot
f = figure;
hold on
units_in_groups = [];
mice_in_groups = [];
for g = 1:length(fields_in_groups)
    subplot(3,ceil(length(fields_in_groups)/3),g) %change the dimensions based on how many groups I have
    hold on
    subtract_unit = 0;
    mice_days = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        if isnan(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_unit_dac) == 0
            plot(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit(:,1), IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit(:,2))
            mice_days = vertcat(mice_days, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.expday);
        else
            subtract_unit = subtract_unit + 1;
        end
    end
    plot(groups.(fields_in_groups{g}).BW_VRC_of_unit(:,1), groups.(fields_in_groups{g}).BW_VRC_of_unit(:,2), 'k', 'LineWidth', 2)
    xlim([0, 0.6]);
    actual_mice = length(unique(mice_days));
    actual_units = groups.(fields_in_groups{g}).units_num - subtract_unit;
    text(.02,.95, [num2str(actual_units) 'units' '\newline' num2str(actual_mice) 'mice'],'Units','normalized')
    title(fields_in_groups(g))
    units_in_groups = vertcat(units_in_groups, actual_units);
    mice_in_groups = vertcat(mice_in_groups, actual_mice);
end
linkaxes()
han = axes(f,'visible','off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
ylabel(han,'Firing Rates (Hz)');

% all mean VRCs on the same subplot
figure
hold on
colors = turbo(length(fields_in_groups));
for g = 1:length(fields_in_groups)
    % true velocities used
    e = errorbar([0, 65, 195, 326, 456, 587, 797], groups.(fields_in_groups{g}).BW_VRC_of_unit(:,2), groups.(fields_in_groups{g}).BW_VRC_of_unit(:,3), 'LineWidth', 2);
    e.Color = colors(g,:);
    text(.78,g*.06, [num2str(units_in_groups(g)) 'units ' num2str(mice_in_groups(g)) 'mice'],'Units','normalized', 'Color', colors(g,:))
end
legend(fields_in_groups, 'Location', 'northwest', 'Box', 'off')
xlabel('Velocity')
ylabel('Firing Rate (Hz)')
ylim([-.5 27.5])
xlim([0 800])
end
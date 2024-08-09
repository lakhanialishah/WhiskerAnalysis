function waveform_plots(groups, IDlist, fields_in_groups)
% graphing the waveforms for all units

% individual units displayed on each subplot
f = figure;
hold on 
for g = 1:length(fields_in_groups)
    subplot(3,ceil(length(fields_in_groups)/3),g) %change the dimensions based on how many groups I have
    hold on
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        plot(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.mwf)
    end
    plot(groups.(fields_in_groups{g}).mwf, 'k', 'LineWidth', 2) 
    text(.7,.2, [num2str(groups.(fields_in_groups{g}).units_num) 'units' '\newline' num2str(groups.(fields_in_groups{g}).mice_num) 'mice'],'Units','normalized')
    title(fields_in_groups(g))
end
linkaxes()

% all mean waveforms on the same plot
figure
hold on 
colors = turbo(length(fields_in_groups));
for g = 1:length(fields_in_groups)
    plot(groups.(fields_in_groups{g}).mwf(:,1), 'LineWidth', 2, 'Color', colors(g,:))
end
legend(fields_in_groups)

end
function spont_fr_plots(groups, IDlist, fields_in_groups)

figure
hold on
all_spont_fr = [];
for g = 1:length(fields_in_groups)
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        spont_fr = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spont_fr;
        all_spont_fr = vertcat(all_spont_fr, spont_fr);
    end
    groups.(fields_in_groups{g}).spont_fr_table = all_spont_fr;
    c = cdfplot(all_spont_fr);
    c.LineWidth = 2;
end
legend(fieldnames(groups))
grid off
title('Cumulative Distribution of Spontaneous Firing Rates')
ylabel('Cumulative Fraction')
xlabel('Spontaneous Firing Rate (Hz)')
clear all_spont_fr

figure
hold on
for g = 1:length(fields_in_groups)
    all_spont_fr = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        spont_fr = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spont_fr;
        all_spont_fr = vertcat(all_spont_fr, spont_fr);
    end
    hold on
    s = swarmchart(ones(length(all_spont_fr),1)*g, all_spont_fr, 'filled', 'MarkerFaceAlpha',0.5);
    s.XJitterWidth = .2;
    plot(g,mean(all_spont_fr), '.k', 'MarkerSize', 20)
    clear spontfr all_spont_fr
end
xticks(1:length(fieldnames(groups)))
xticklabels(fields_in_groups)
title('Spontaneous firing rates for all units')
ylabel('Spontaneous Firing Rate (Hz)')
xlabel('Experiment Condition')

end


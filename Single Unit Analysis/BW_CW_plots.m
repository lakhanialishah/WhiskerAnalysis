function BW_CW_plots(groups, IDlist, fields_in_groups)
% graphing data to see if BW equals CW 
% so looking at unit data, NOT BW for the day data 
% also only looking at days when the probe was directly in a barrel, not septa or unclear 

figure;
hold on 

% colors 
ko = [0.9290, 0.6940, 0.1250];
wt = [0, 0.4470, 0.7410];
wd = [0.9, 0.0580, 0.1];
colors = [ko ; wd];

for g = 1:length(fields_in_groups)
    bar(g, groups.(fields_in_groups{g}).BW_is_CW_proportion, 'FaceColor', colors(g,:))
    mice_days = [];
    all_answers = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        answer = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_CW; 
        all_answers = vertcat(all_answers, answer);
        mice_days = vertcat(mice_days, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.expday);
    end
    indices = find(isnan(all_answers) == 1); 
    all_answers(indices) = [];
    mice_days(indices) = [];
    if isempty(all_answers) == 0 
        units_num = length(all_answers); %length(find(all_answers == 0 | all_answers == 1));
    else 
        units_num = 0; 
    end
    text(g-.15, .15, [num2str(units_num)])
    text(g-.15, .05, [num2str(length(unique(mice_days)))])
    clear units_num
end

xticks(1:length(fieldnames(groups)))
xticklabels(fields_in_groups)
ylabel('Proportion of units where CW = BW')
xlabel('Experiment Condition')
ylim([0 1])
end
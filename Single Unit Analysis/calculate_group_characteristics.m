function groups = calculate_group_characteristics(fields_in_groups, g, groups, IDlist)
% calculate all the things I want for the group/experimental condition data

%spontaneous firing rate
all_spont_fr = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    spont_fr = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spont_fr;
    all_spont_fr = vertcat(all_spont_fr, spont_fr);
end
groups.(fields_in_groups{g}).spont_fr = mean(all_spont_fr);
groups.(fields_in_groups{g}).spont_fr_sem = std(all_spont_fr)/ sqrt(length(all_spont_fr));

%mean waveform
all_wf = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    wf = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.mwf;
    all_wf = horzcat(all_wf, wf);
end
groups.(fields_in_groups{g}).mwf = mean(all_wf,2);

%evoked firing rate from step2 only
all_evoked_firing_rates = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    evoked_fr = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.step2_VRC;
    all_evoked_firing_rates = vertcat(all_evoked_firing_rates, evoked_fr);
end
unique_velocities = unique(all_evoked_firing_rates(:,1));
group_step2_VRC(:,1) = unique_velocities;
for vel = 1: length(unique_velocities)
    indices = find(all_evoked_firing_rates(:,1) == unique_velocities(vel));
    group_step2_VRC(vel, 2) = mean(all_evoked_firing_rates(indices,2));
    group_step2_VRC(vel, 3) = (std(all_evoked_firing_rates(indices,2)))/sqrt(length(all_evoked_firing_rates(indices,2)));
end
groups.(fields_in_groups{g}).step2_VRC = group_step2_VRC; %velocity, mean evoked firing rate, standard deviation of evoked firing rate

%total number of units
units_num = length(groups.(fields_in_groups{g}).IDlist_indices);
groups.(fields_in_groups{g}).units_num = units_num;

%total number of mice
all_exp_days = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    day = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.expday;
    all_exp_days = vertcat(all_exp_days, day);
end
unique_days = unique(all_exp_days);
groups.(fields_in_groups{g}).mice_num = length(unique_days);

%mean VRC for group using BW VRC from each unit data
all_VRCs = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    if isempty(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit) == 0
        fr = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit;
        all_VRCs = vertcat(all_VRCs, fr);
    else
        continue
    end
end
%to remove units that don't have BW 
the_answer = double(isnan(all_VRCs(:,2))); %nan values 
indices = find(the_answer == 1); %to remove
all_VRCs(indices,:) = [];
% for some reason this doesnt work anymore
% unique_velocities = unique(all_VRCs(:,1));
velocities = [0, .1, .2, .3, .4, .5, .6];
group_BW_VRC(:,1) = velocities;
for vel = 1: length(velocities)
    indices = vel:7:length(all_VRCs);
    frs = all_VRCs(indices,2);
    group_BW_VRC(vel, 2) = nanmean(frs);
    group_BW_VRC(vel, 3) = nanstd(frs)/ (sqrt(length(frs))); %this is the sem
    group_BW_VRC(vel, 4) = nanstd(frs); %this is the std
end
groups.(fields_in_groups{g}).BW_VRC_of_unit = group_BW_VRC; %velocity, mean evoked firing rate, sem of evoked firing rate

%calculating percentage of when BW = CW in each group 
all_answers = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    answer = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_CW;
    all_answers = vertcat(all_answers, answer);
end
indices = find(isnan(all_answers) == 1);
all_answers(indices) = [];
proportion = (sum(all_answers))/(length(all_answers));
groups.(fields_in_groups{g}).BW_is_CW_proportion = proportion;

%calcuating CW VRC for each group (only when the probe is in a barrel and when the CW is in a peizo) 
all_CW_VRC = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    CW_VRC = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC;
    CW_VRC = CW_VRC';
    all_CW_VRC = vertcat(all_CW_VRC, CW_VRC);
end
the_average = nanmean(all_CW_VRC,1)';
if length(groups.(fields_in_groups{g}).IDlist_indices) <= 1 %if there is just one unit in the group category
    the_sem = [0 0 0 0 0 0 0]';
else
    the_sem = (nanstd(all_CW_VRC,1) / (sqrt(size(all_CW_VRC, 1))))';
end
groups.(fields_in_groups{g}).CW_VRC = [the_average the_sem]; %mean CW VRC, sem of CW VRC

%table with BW_VRC_for_unit specifically for stats 
all_VRCs = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    if isempty(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit) == 0
        fr = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit(:,2);
        all_VRCs = horzcat(all_VRCs, fr);
    else
        continue
    end
end
indices = [];
for i = 1:size(all_VRCs,2)
    if isnan(all_VRCs(1,i)) == 1
        indices = vertcat(indices, i);
    end
end
all_VRCs(:,indices) = [];
groups.(fields_in_groups{g}).BW_VRC_of_unit_table = all_VRCs;

%table with CW_VRC_for_unit specifically for stats 
all_CW_VRCs = [];
for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
    if isempty(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC) == 0
        fr = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC;
        all_CW_VRCs = horzcat(all_CW_VRCs, fr);
    else
        continue
    end
end
indices = [];
for i = 1:size(all_CW_VRCs,2)
    if isnan(all_CW_VRCs(1,i)) == 1
        indices = vertcat(indices, i);
    end
end
all_CW_VRCs(:,indices) = [];
groups.(fields_in_groups{g}).CW_VRC_of_unit_table = all_CW_VRCs;
end
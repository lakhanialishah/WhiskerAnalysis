%% final script
% final script with all functions  
% units saved from the postphy scripts

%% load single units

folder = 'C:\Users\lakha\Documents\Units'; %strcat('C:\Users\lakha\Documents\Units'); is the default
list = dir(folder);
list(1:2) = [];
IDlist = char();

for i = 1:length(list)
    IDlist{i,1} = extractBefore(list(i).name, '.mat');
end

for m = 1:length(list)
    variable = load(fullfile(list(m).folder, list(m).name));
    IDlist{m, 2} = variable;
end

%% KSfile_data_extraction and assigning to clusters

expsummaryfile = 'C:\Users\lakha\Documents\ExperimentsSummary';
exp_files = table2cell(readtable(expsummaryfile, 'Sheet', 1));

% FOR WHEN I HAVE NEW DATA TO EXTRACT
% [stim_matrices_for_each_day] = KSfile_data_extraction(exp_files);
%and then save the variable in the singleunitanalysis scripts folder

% format for stim_matrices_for_each_day:
%     column1 = date; column2 = step1 matrix; column3 = step2 matrix; column4 = step2_2 matrix; column5 = step3 matrix
%     column6 = all stim matrix, column7 = experimental condition, column8 = BW, column9 = BWdac, column10 = probe position, column11 = CW_dac

load('stim_matrices_for_each_day.mat');

for day = 1:size(stim_matrices_for_each_day,1)
    for clu = 1: size(IDlist, 1)
        if stim_matrices_for_each_day{day,1} == IDlist{clu,2}.clu_info.expday
            IDlist{clu,2}.clu_info.step1_matrix = stim_matrices_for_each_day{day,2};
            IDlist{clu,2}.clu_info.step2_matrix = stim_matrices_for_each_day{day,3};
            IDlist{clu,2}.clu_info.step2_2_matrix = stim_matrices_for_each_day{day,4};
            IDlist{clu,2}.clu_info.step3_matrix = stim_matrices_for_each_day{day,5};
            IDlist{clu,2}.clu_info.allstim_matrix = stim_matrices_for_each_day{day,6};
            IDlist{clu,2}.clu_info.exp_cond = stim_matrices_for_each_day{day,7};
            IDlist{clu,2}.clu_info.BW_day = stim_matrices_for_each_day{day,8};
            IDlist{clu,2}.clu_info.BW_day_dac = stim_matrices_for_each_day{day,9};
            IDlist{clu,2}.clu_info.probe_position = stim_matrices_for_each_day{day,10};
            IDlist{clu,2}.clu_info.CW_dac = stim_matrices_for_each_day{day,11};
        end
    end
end

%% determine what layer a unit is in
% started using cambridge probe on 221005

csdfilename = 'C:\Users\lakha\Documents\LayerChansFromCSD';
chanfile = table2array(readtable(csdfilename, 'Sheet', 5));

for i = 1: size(IDlist,1)
    date = IDlist{i,2}.clu_info.expday;
    idx = find(date == chanfile(:,1)); %index for the date in chanfile
    if IDlist{i,2}.clu_info.truechan > chanfile(idx,2) & IDlist{i,2}.clu_info.truechan < chanfile (idx, 3)
        IDlist{i,2}.clu_info.layer = 'L4';
    elseif IDlist{i,2}.clu_info.truechan > chanfile(idx,4) & IDlist{i,2}.clu_info.truechan < chanfile (idx, 5)
        IDlist{i,2}.clu_info.layer = 'L2/3';
    elseif IDlist{i,2}.clu_info.truechan > chanfile(idx,6) & IDlist{i,2}.clu_info.truechan < chanfile (idx, 7)
        IDlist{i,2}.clu_info.layer = 'L5/6';
    else
        IDlist{i,2}.clu_info.layer = 'not assigned';
    end
end

%% spontaneous firing rate for each unit using spont_firing_rate function

for i = 1: size(IDlist,1)
    allstim_matrix = IDlist{i,2}.clu_info.allstim_matrix;
    spike_times = IDlist{i,2}.clu_info.spikeTimes;
    spont_fr = spont_firing_rate(allstim_matrix, spike_times);
    IDlist{i,2}.clu_info.spont_fr = spont_fr;
end

%% VRC for each unit using VRC_for_each_unit function for step2 (and this is using the BW of the day)

for i = 1: size(IDlist,1)
    step2_matrix = IDlist{i,2}.clu_info.step2_matrix;
    spike_times = IDlist{i,2}.clu_info.spikeTimes;
    step2_VRC = VRC_for_each_unit(step2_matrix, spike_times);
    IDlist{i,2}.clu_info.step2_VRC = step2_VRC;
end

%% giving each unit the BW VRC and the BW for the unit (so not the BW for the day, but the BW for the unit)

for i = 1:size(IDlist,1)
    step3_matrix = IDlist{i,2}.clu_info.step3_matrix;
    spike_times = IDlist{i,2}.clu_info.spikeTimes;
    BW_day_dac = IDlist{i,2}.clu_info.BW_day_dac;
    [BW_VRC_unit, BW_dac_unit] = BW_VRC_for_unit(step3_matrix, spike_times, BW_day_dac);
    IDlist{i,2}.clu_info.BW_VRC_for_unit = BW_VRC_unit;
    IDlist{i,2}.clu_info.BW_unit_dac = BW_dac_unit;
end

%% is BW for the unit the CW for the unit

for clu = 1:size(IDlist,1)
    if IDlist{clu,2}.clu_info.CW_dac == IDlist{clu,2}.clu_info.BW_unit_dac
        IDlist{clu,2}.clu_info.BW_CW = 1; %means the BW and the CW are the same dac
    elseif isnan(IDlist{clu,2}.clu_info.BW_unit_dac) == 1
        IDlist{clu,2}.clu_info.BW_CW = NaN; %there is no BW
    elseif IDlist{clu,2}.clu_info.CW_dac == 0
        IDlist{clu,2}.clu_info.BW_CW = NaN; %there is no CW
    elseif IDlist{clu,2}.clu_info.CW_dac ~= IDlist{clu,2}.clu_info.BW_unit_dac
        IDlist{clu,2}.clu_info.BW_CW = 0; %CW does not equal BW
    end
end

%% giving each unit where there is a CW in the first place a CW_VRC from step3
%so this is only when probe ended up in a barrel

for i = 1:size(IDlist,1)
    step3_matrix = IDlist{i,2}.clu_info.step3_matrix;
    spike_times = IDlist{i,2}.clu_info.spikeTimes;
    CW_dac = IDlist{i,2}.clu_info.CW_dac;
    CW_VRC = CW_VRC_calculation(step3_matrix, spike_times, CW_dac);
    IDlist{i,2}.clu_info.CW_VRC = CW_VRC;
end

%% ISI elimination round using refractory_violation_rate function
%function takes refractory time bin I'm interested in and spike times of cluster
for clu = 1: size(IDlist,1)
    clu_spike_times = IDlist{clu,2}.clu_info.spikeTimes;
    IDlist{clu,2}.clu_info.violation_rate_1ms = refractory_violation_rate(1, clu_spike_times);
end

cutoff = 1.5; %this is a % cutoff 

for i = 1:size(IDlist,1)
    if IDlist{i,2}.clu_info.violation_rate_1ms <= cutoff %looking at the 1 ms bin
        keep{i,:} = IDlist{i,:};
    else
        keep{i,:} = NaN;
    end
end

% to delete the units from the list so that we can move on to the next parameter
index = cellfun(@isnan,keep,'uni',false);
index = cellfun(@any,index);
keep(index) = [];

clear cutoff index

%% amplitude: gaussian fit/missed spikes elimination round using missed_spikes_based_on_amplitude function

for clu = 1: size(IDlist,1)
    clu_amplitudes = IDlist{clu,2}.clu_info.amplitudes;
    IDlist{clu,2}.clu_info.missed_spikes = missed_spikes_based_on_amplitude(clu_amplitudes);
end

cutoff = 20; %percent of missed spikes

for clu = 1:size(keep,1)
    %to find the index of that cluster in the IDlist
    findind = strfind(IDlist(:,1), keep{clu,1});
    cluindex = [];
    for k = 1:numel(findind)
        if findind{k} == 1
            cluindex = k ;
        end
    end
    
    if IDlist{cluindex,2}.clu_info.missed_spikes <= cutoff
        keep2{clu,:} = IDlist{cluindex,1};
    else
        keep2{clu,:} = NaN;
    end
end

% to delete the units from the list so that we can move on to the next parameter
index = cellfun(@isnan,keep2,'uni',false);
index = cellfun(@any,index);
keep2(index) = [];

clear cutoff index

%% amplitude: most spikes above the 10.6 uV detection threshold

cutoff = 3; %percent of spikes below the detection threshold
threshold = 11;

for clu = 1:size(keep2,1)
    %code to find the index of that cluster in the IDlist
    findind = strfind(IDlist(:,1), keep2{clu,1});
    cluindex = [];
    for k = 1:numel(findind)
        if findind{k} == 1
            cluindex = k ;
        end
    end
    
    clu_amplitudes = IDlist{cluindex,2}.clu_info.amplitudes;
    percent_below_threshold = length(find(clu_amplitudes < threshold))/length(clu_amplitudes) * 100;
    
    if percent_below_threshold <= cutoff
        keep3{clu,:} = IDlist{cluindex,1};
    else
        keep3{clu,:} = NaN;
    end
end

% to delete the units from the list so that we can move on to the next parameter
index = cellfun(@isnan,keep3,'uni',false);
index = cellfun(@any,index);
keep3(index) = [];

clear cutoff index

%% mean firing rate criteria using mean_firing_rate function

for clu = 1: size(IDlist,1)
    clu_spike_times = IDlist{clu,2}.clu_info.spikeTimes;
    [binning_firing_rate, mean_firing_rate] = mean_firing_rate_calculation(clu_spike_times);
    indices = find(binning_firing_rate < mean_firing_rate * .20); % 20% below mean firing rate
    total = length(indices) / length(binning_firing_rate) * 100;
    IDlist{clu,2}.clu_info.mean_firing_rate = mean_firing_rate;
    IDlist{clu,2}.clu_info.percent_below_mean_firing_rate = total;
end

cutoff = 10; %this is a % cutoff

for clu = 1:size(keep3,1)
    %code to find the index of that cluster in the IDlist
    findind = strfind(IDlist(:,1), keep3{clu,1});
    cluindex = [];
    for k = 1:numel(findind)
        if findind{k} == 1
            cluindex = k ;
        end
    end
    
    if IDlist{cluindex,2}.clu_info.percent_below_mean_firing_rate <= cutoff
        keep4{clu,:} = IDlist{cluindex,1};
    else
        keep4{clu,:} = NaN;
    end
end

% to delete the units from the list so that we can move on to the next parameter
index = cellfun(@isnan,keep4,'uni',false);
index = cellfun(@any,index);
keep4(index) = [];

clear cutoff index

%% assigning FS vs RS and graphing all waveforms I will use 

IDlist = assignFSvRS(keep4, IDlist); %blue iS FS, red is RS

%% separating good units based on experimental groups using the separating_units_into_conditions function

groups = separating_units_into_conditions(keep4, IDlist);

%% calculating group data and saving it in the groups structure

fields_in_groups = fieldnames(groups);
for g = 1:length(fields_in_groups)
    groups = calculate_group_characteristics(fields_in_groups, g, groups, IDlist);
end

%% graphing waveforms in groups

fields_in_groups = fieldnames(groups);
indices = find((contains(fields_in_groups, 'L5') == 1 & contains(fields_in_groups, 'FS') == 1 ) == 1); %indices within groups for a specific subgroup & contains(fields, 'k') == 1
% indices = find(contains(fields_in_groups, 'L5kwd0_21')); %to put just the one on graph
fields_in_groups = fields_in_groups(indices);
waveform_plots(groups, IDlist, fields_in_groups)

%% plotting spontaneous firing rate with spont_fr_plots function

%specific groups on one graph
fields_in_groups = fieldnames(groups);
indices = find((contains(fields_in_groups, 'L5') == 1 & contains(fields_in_groups, 'FS') == 0 ) == 1); %indices within groups for a specific subgroup & contains(fields, 'k') == 1
fields_in_groups = fields_in_groups(indices);
% indices = [find(contains(fields_in_groups, 'L5cwd0_21')) find(contains(fields_in_groups, 'L5kwd0_21'))]; %to put just the 2 on one graph
% fields_in_groups = fields_in_groups(indices);
spont_fr_plots(groups, IDlist, fields_in_groups)

%putting fields_in_groups back to normal
fields_in_groups = fieldnames(groups);

%% plotting evoked firing rate with evoked_fr_plots_BW_of_unit function using BW of the UNIT

%for all experimental conditions on one graph
evoked_fr_plots_BW_of_unit(groups, IDlist, fields_in_groups)

%specific groups on one graph
fields_in_groups = fieldnames(groups);
% indices = [find(contains(fields_in_groups, 'L5cwd0_21')) find(contains(fields_in_groups, 'L5kwd0_21'))]; %to put just the 2 on one graph
indices = find((contains(fields_in_groups, 'L5') == 1 & contains(fields_in_groups, 'FS') == 1 ) == 1); %indices within groups for a specific subgroup & contains(fields, 'k') == 1
fields_in_groups = fields_in_groups(indices);
evoked_fr_plots_BW_of_unit(groups, IDlist, fields_in_groups)

%putting fields_in_groups back to normal
fields_in_groups = fieldnames(groups);

%% plotting proportions of if BW is the CW for the unit

fields_in_groups = fieldnames(groups);
indices = find((contains(fields_in_groups, 'L5') == 1 & contains(fields_in_groups, 'FS') == 0 ) == 1); %indices within groups for a specific subgroup & contains(fields, 'k') == 1
% indices = [find(contains(fields_in_groups, 'L5kwd0_21')) find(contains(fields_in_groups, 'L5kwd7'))]; %to put just the 2 on one graph
fields_in_groups = fields_in_groups(indices);
BW_CW_plots(groups, IDlist, fields_in_groups)

%% plotting CW_VRCs (responders only and all combined graphs)

%for some specific groups on one graph
fields_in_groups = fieldnames(groups);
indices = find((contains(fields_in_groups, 'L5') == 1 & contains(fields_in_groups, 'FS') == 1 ) == 1); %indices within groups for a specific subgroup & contains(fields, 'k') == 1
fields_in_groups = fields_in_groups(indices);
% indices = [find(contains(fields_in_groups, 'L5FSkwd0_21')) find(contains(fields_in_groups, 'L5FSkwd7'))]; %to put just the 2 on one graph
% fields_in_groups = fields_in_groups(indices);
[groups, IDlist] = CW_VRC_plot(groups, IDlist, fields_in_groups);

%putting fields_in_groups back to normal
fields_in_groups = fieldnames(groups);

%% plotting VRCs for one experimental condition, but individual VRCs for each mouse, with all units from that day on one graph

fields_in_groups = fieldnames(groups);

% manually changing the number to get the experimental condition I want to look at throughout this whole section of script
all_dates = [];
for unit = 1:length(groups.(fields_in_groups{9}).IDlist_indices)
    date = IDlist{groups.(fields_in_groups{9}).IDlist_indices(unit,1),2}.clu_info.expday;
    all_dates = vertcat(all_dates, date);
end
unique_dates = unique(all_dates);

figure
for d = 1:length(unique_dates)
    subplot(2,5,d)
    title(unique_dates(d))
    hold on
    for unit = 1:length(groups.(fields_in_groups{9}).IDlist_indices)
        date = IDlist{groups.(fields_in_groups{9}).IDlist_indices(unit,1),2}.clu_info.expday;
        if date == unique_dates(d)
            plot(IDlist{groups.(fields_in_groups{9}).IDlist_indices(unit,1),2}.clu_info.CW_VRC)
        else
            continue
        end
    end
end
linkaxes()

%% power analysis to determine how many units I need to include superficial layers

% using the middle velocity (because there can be shape difference) for the power analysis, manually change the experimental condition
mean1 = groups.L5cwd0_16.BW_VRC_of_unit(4,2);
mean2 = groups.L5cwd2.BW_VRC_of_unit(4,2);
std1 = groups.L5cwd0_16.BW_VRC_of_unit(4,4);
std2 = groups.L5cwd2.BW_VRC_of_unit(4,4);
difference = abs(mean1 - mean2);

%one sample t-test, because only one part of the conditions is changing with all our comparisons
nout = sampsizepwr('t', [mean1 std1], mean2, 0.80); 

%% threshold of each whisker-responsive unit 
%comment out if I'm using BW or CWresp

%going through the units
for clu = 1:size(IDlist,1)
    step3_matrix = IDlist{clu,2}.clu_info.step3_matrix;
    spike_times = IDlist{clu,2}.clu_info.spikeTimes;
    %with using CW responders
%     if isfield(IDlist{clu,2}.clu_info, 'isCWresp') == 1
%         if IDlist{clu,2}.clu_info.isCWresp == 1
%             interested_dac = IDlist{clu,2}.clu_info.CW_dac;
%         end
%     end
%     if isnan(interested_dac) == 1
%         threshold = NaN;
%         IDlist{clu,2}.clu_info.vel_thresh_for_unit_with_CWresp = threshold;
%         continue
%     else
%         threshold = vel_threshold_for_unit(step3_matrix, spike_times, interested_dac);
%         IDlist{clu,2}.clu_info.vel_thresh_for_unit_with_CWresp = threshold;
%     end
    
    %with using BW
    interested_dac = IDlist{clu,2}.clu_info.BW_unit_dac; %could be BW or CW responders based on what I'm interested in
    if isnan(interested_dac) == 1
        threshold = NaN;
        IDlist{clu,2}.clu_info.vel_thresh_for_unit_with_BW = threshold;
        continue
    else
        threshold = vel_threshold_for_unit(step3_matrix, spike_times, interested_dac);
        IDlist{clu,2}.clu_info.vel_thresh_for_unit_with_BW = threshold;
    end
end

%going through the experimental groups to get a mean and std
fields_in_groups = fieldnames(groups);
for g = 1:length(fields_in_groups)
    all_answers = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        threshold = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.vel_thresh_for_unit_with_BW; %either with_BW or with_CWresp
        if isnan(threshold) == 1 %to make sure no nans are included
            continue
        end
        all_answers = vertcat(all_answers, threshold);
    end
    mean_threshold = mean(all_answers); %taking the mean of the thresholds
    groups.(fields_in_groups{g}).velocity_threshold = mean_threshold;
    groups.(fields_in_groups{g}).velocity_threshold_std = std(all_answers);
end

% plotting the threshold of groups, units threshold as dots
fields_in_groups = fieldnames(groups);
indices = find((contains(fields_in_groups, 'L5') == 1 & contains(fields_in_groups, 'FS') == 0 ) == 1); %indices within groups for a specific subgroup & contains(fields, 'k') == 1
fields_in_groups = fields_in_groups(indices);
figure 
hold on 
all_threshold_matrix_boxplot = [];
groups_boxplot = {};
for g = 1:length(fields_in_groups)
    all_thresholds = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        threshold = IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.vel_thresh_for_unit_with_CWresp; %either with_BW or with_CWresp
        if isnan(threshold) == 1 
            continue
        end
        all_thresholds = vertcat(all_thresholds, threshold);
    end
    hold on
    s = swarmchart(ones(length(all_thresholds),1)*g, all_thresholds, 'filled', 'MarkerFaceAlpha',0.5);
    s.XJitterWidth = .8;
    plot(g,mean(all_thresholds), '.k', 'MarkerSize', 20)
    
    all_threshold_matrix_boxplot = vertcat(all_threshold_matrix_boxplot, all_thresholds);
    groups_boxplot = vertcat(groups_boxplot, repmat({fields_in_groups{g}},length(all_thresholds),1));
    clear threshold all_thresholds
end
boxplot(all_threshold_matrix_boxplot, groups_boxplot, 'Whisker',0, 'Symbol', '') %to make the whiskers go away and hide the outlier symbol
xticks(1:length(fieldnames(groups)))
xticklabels(fields_in_groups)
ylabel('Velocity threshold of individual units')

%% to write data onto an excel file

fields_in_groups = fieldnames(groups);
indices = find((contains(fields_in_groups, 'L5') == 1 & contains(fields_in_groups, 'FS') == 0 ) == 1); %indices within groups for a specific subgroup & contains(fields, 'k') == 1
%indices = [find(contains(fields_in_groups, 'L5cwd0_21')) find(contains(fields_in_groups, 'L5cwd7'))]; %to put just the 2 on one graph
fields_in_groups = fields_in_groups(indices);

%BW of the unit
for g = 1:length(fields_in_groups)
    tab = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        if isempty(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit) == 0
            tab = vertcat(tab, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_VRC_for_unit(:,2)');
        end
    end
    writematrix(tab, 'C:\Users\lakha\Documents\ExcelSheetsForPrism\BW_unit.xlsx', 'Sheet', fields_in_groups{g})
    clear tab
end

%CW of the unit
for g = 1:length(fields_in_groups)
    writematrix(groups.(fields_in_groups{g}).CW_VRC_of_unit_table', 'C:\Users\lakha\Documents\ExcelSheetsForPrism\CW_unit.xlsx', 'Sheet', fields_in_groups{g})
end

%proportion of responders with CW stimulation
for g = 1:length(fields_in_groups)
    writematrix(groups.(fields_in_groups{g}).CW_proportion_resp, 'C:\Users\lakha\Documents\ExcelSheetsForPrism\CW_proportion_responsive.xlsx', 'Sheet', fields_in_groups{g})
end

%CW VRCs of responders only
for g = 1:length(fields_in_groups)
    writematrix(groups.(fields_in_groups{g}).CW_VRC_responders_table', 'C:\Users\lakha\Documents\ExcelSheetsForPrism\CW_unit_responders.xlsx', 'Sheet', fields_in_groups{g})
end

%spontaneous activity 
for g = 1:length(fields_in_groups)
    tab = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        if isnan(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spont_fr) == 0
            tab = vertcat(tab, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spont_fr);
        end
    end
    writematrix(tab, 'C:\Users\lakha\Documents\ExcelSheetsForPrism\Spont_fr.xlsx', 'Sheet', fields_in_groups{g})
    clear tab
end

%spontaneous activity of units that have a BW 
for g = 1:length(fields_in_groups)
    tab = [];
        for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
            if isnan(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.BW_unit_dac) == 0
                tab = vertcat(tab, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spont_fr);
            end
        end
    writematrix(tab, 'C:\Users\lakha\Documents\ExcelSheetsForPrism\Spont_fr_BWonly.xlsx', 'Sheet', fields_in_groups{g})
    clear tab
end

%spontaneous activity of units that have a CW (responders and not)
for g = 1:length(fields_in_groups)
    tab = [];
        for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
            if isnan(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.CW_VRC) == 0
                tab = vertcat(tab, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.spont_fr);
            end
        end
    writematrix(tab, 'C:\Users\lakha\Documents\ExcelSheetsForPrism\Spont_fr_CWall.xlsx', 'Sheet', fields_in_groups{g})
    clear tab 
end

%threshold velocity for units 
for g = 1:length(fields_in_groups)
    tab = [];
    for unit = 1:length(groups.(fields_in_groups{g}).IDlist_indices)
        if isnan(IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.vel_thresh_for_unit_with_CWresp) == 0 %either with_BW or with_CWresp
            tab = vertcat(tab, IDlist{groups.(fields_in_groups{g}).IDlist_indices(unit,1),2}.clu_info.vel_thresh_for_unit_with_CWresp); %either with_BW or with_CWresp
        end
    end
    writematrix(tab, 'C:\Users\lakha\Documents\ExcelSheetsForPrism\Velocity_threshold.xlsx', 'Sheet', fields_in_groups{g})
    clear tab
end

%proportion where BW = CW
for g = 1:length(fields_in_groups)
    writematrix(groups.(fields_in_groups{g}).BW_is_CW_proportion, 'C:\Users\lakha\Documents\ExcelSheetsForPrism\BW_CW_proportion.xlsx', 'Sheet', fields_in_groups{g})
end
function [stim_matrices_for_each_day] = KSfile_data_extraction(exp_files)
% to extract stim time, vel, and dac from all files in the KS file
% takes the input argument of the file name
% output is matrices for all steps including in KS file
% if the step was not used, then it spits out an empty vector/variable for it
% format for stim_matrices_for_each_day:
%     column1 = date
%     column2 = step1 matrix
%     column3 = step2 matrix
%     column4 = step2_2 matrix
%     column5 = step3 matrix
%     column6 = all stim matrix
%     column7 = experimental condition
%     column8 = BW
%     column9 = BWdac
%     column10 = probe position
%     column11 = CWdac

% file = 'C:\Users\lakha\Documents\ExperimentsSummary';
% exp_files = table2cell(readtable(file, 'Sheet', 1));
fs = 24414.0625;

stim_matrices_for_each_day = {};
for day = 1: size(exp_files,1)
    stim_matrices_for_each_day{day,1} = exp_files{day,1};
end

for day = 1:size(exp_files,1)
    fprintf(' day %i', day); %printing what number experiment it's on
    if isnan(exp_files{day,1}) == 1
        continue
    else
        if day ~= 16 || 21 || 32 || 48 || 49 || 52|| 53 || 57 || 61
            allstims = [];
            %step1
            if isempty(exp_files{day, 2}) == 0
                step1_file_name = exp_files{day,2};
                step1_data = TDTbin2mat(step1_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
                step1_time = size(step1_data.streams.Raws.data, 2)/fs;
                step1_dacs = step1_data.epocs.Us1_.data;
                step1_vels = ones(size(step1_data.epocs.Us1_.data))*(max(step1_data.streams.sRmp.data));
                step1_stim_times = step1_data.epocs.Us1_.onset;
                step1 = [step1_dacs step1_vels step1_stim_times];
                allstims = vertcat(allstims, step1);
            else
                step1 = [];
            end
            stim_matrices_for_each_day{day,2} = step1;
            %step2
            if isempty(exp_files{day, 3}) == 0
                step2_file_name = exp_files{day,3};
                step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
                step2_time = size(step2_data.streams.Raws.data, 2)/fs;
                step2_dacs = step2_data.epocs.DacN.data;
                step2_vels = step2_data.epocs.RmpV.data;
                if isempty(exp_files{day,2}) == 1 %if there is no step1
                    step2_stim_times = step2_data.epocs.RmpV.onset;
                else %if there is step1
                    step2_stim_times = step2_data.epocs.RmpV.onset + step1_time;
                end
                step2 = [step2_dacs step2_vels step2_stim_times];
                allstims = vertcat(allstims, step2);
            else
                step2 = [];
            end
            stim_matrices_for_each_day{day,3} = step2;
            %step2_2
            if isempty(exp_files{day, 4}) == 0
                step2_2_file_name = exp_files{day,4};
                step2_2_data = TDTbin2mat(step2_2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
                step2_2_time = size(step2_2_data.streams.Raws.data, 2)/fs;
                step2_2_dacs = step2_2_data.epocs.DacN.data;
                step2_2_vels = step2_2_data.epocs.RmpV.data;
                if isempty(exp_files{day,2}) == 1 %if there is no step1
                    step2_2_stim_times = step2_2_data.epocs.RmpV.onset + step2_time;
                else %if there is step1
                    step2_2_stim_times = step2_2_data.epocs.RmpV.onset + step2_time + step1_time;
                end
                step2_2 = [step2_2_dacs step2_2_vels step2_2_stim_times];
                allstims = vertcat(allstims, step2_2);
            else
                step2_2 = [];
            end
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3
            if isempty(exp_files{day, 5}) == 0
                step3_file_name = exp_files{day,5};
                step3_data = TDTbin2mat(step3_file_name, 'TYPE', {'snips', 'epocs'}); %not bringing in the streams because don't need to calulate the file length of step3
                step3_dacs = step3_data.epocs.DacN.data;
                step3_vels = step3_data.epocs.RmpV.data;
                if isempty(exp_files{day,2}) == 1 & isempty(exp_files{day,4}) == 1 %if there is no step1 and no step2_2
                    step3_stim_times = step3_data.epocs.RmpV.onset + step2_time;
                elseif isempty(exp_files{day,2}) == 0 & isempty(exp_files{day,4}) == 1 % if there is step1 and no step2_2
                    step3_stim_times = step3_data.epocs.RmpV.onset + step2_time + step1_time;
                elseif isempty(exp_files{day,2}) == 1 & isempty(exp_files{day,4}) == 0 % if there is no step1 and there is step2_2
                    step3_stim_times = step3_data.epocs.RmpV.onset + step2_time + step2_2_time;
                else isempty(exp_files{day,2}) == 0 & isempty(exp_files{day,4}) == 0 %if there is step1 and step2_2
                    step3_stim_times = step3_data.epocs.RmpV.onset + step1_time + step2_time + step2_2_time;
                end
                step3 = [step3_dacs step3_vels step3_stim_times];
                allstims = vertcat(allstims, step3);
            else
                step3 = [];
            end
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 16 %special case for 230301
            allstims = [];
            %step1_0
            step1_0_file_name = exp_files{day,10};
            step1_0_data = TDTbin2mat(step1_0_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step1_0_time = size(step1_0_data.streams.Raws.data, 2)/fs;
            %step1
            step1_file_name = exp_files{day,2};
            step1_data = TDTbin2mat(step1_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step1_time = size(step1_data.streams.Raws.data, 2)/fs;
            step1_dacs = step1_data.epocs.Us1_.data;
            step1_vels = ones(size(step1_data.epocs.Us1_.data))*(max(step1_data.streams.sRmp.data));
            step1_stim_times = step1_data.epocs.Us1_.onset + step1_0_time;
            step1 = [step1_dacs step1_vels step1_stim_times];
            allstims = vertcat(allstims, step1);
            stim_matrices_for_each_day{day,2} = step1;
            %step2
            step2_file_name = exp_files{day,3};
            step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_time = size(step2_data.streams.Raws.data, 2)/fs;
            step2_dacs = step2_data.epocs.DacN.data;
            step2_vels = step2_data.epocs.RmpV.data;
            step2_stim_times = step2_data.epocs.RmpV.onset + step1_time + step1_0_time;
            step2 = [step2_dacs step2_vels step2_stim_times];
            allstims = vertcat(allstims, step2);
            stim_matrices_for_each_day{day,3} = step2;
            stim_matrices_for_each_day{day,4} = [];
            stim_matrices_for_each_day{day,5} = [];
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 21 %special case for 230621
            allstims = [];
            %step1
            step1_file_name = exp_files{day,2};
            step1_data = TDTbin2mat(step1_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step1_time = size(step1_data.streams.Raws.data, 2)/fs;
            step1_dacs = step1_data.epocs.Us1_.data;
            step1_vels = ones(size(step1_data.epocs.Us1_.data))*(max(step1_data.streams.sRmp.data));
            step1_stim_times = step1_data.epocs.Us1_.onset;
            step1 = [step1_dacs step1_vels step1_stim_times];
            allstims = vertcat(allstims, step1);
            stim_matrices_for_each_day{day,2} = step1;
            %step2
            step2_file_name = exp_files{day,3};
            step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_time = size(step2_data.streams.Raws.data, 2)/fs;
            step2_dacs = step2_data.epocs.DacN.data;
            step2_vels = step2_data.epocs.RmpV.data;
            step2_stim_times = step2_data.epocs.RmpV.onset + step1_time;
            step2 = [step2_dacs step2_vels step2_stim_times];
            allstims = vertcat(allstims, step2);
            stim_matrices_for_each_day{day,3} = step2;
            %step2_2
            step2_2 = [];
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3 (broken into two files for this day)
            part1 = exp_files{day,5};
            part2 = exp_files{day,11};
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 1583.51);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs'});
            part1_dacs = part1_data.epocs.DacN.data;
            part2_dacs = part2_data.epocs.DacN.data;
            dacs = [part1_dacs; part2_dacs];
            part1_vels = part1_data.epocs.RmpV.data;
            part2_vels = part2_data.epocs.RmpV.data;
            vels = [part1_vels; part2_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset;
            part2_stimtimes = part2_data.epocs.RmpV.onset + 1583.51;
            stimtimes = [part1_stimtimes; part2_stimtimes];
            stimtimes = stimtimes + step1_time + step2_time;
            step3 = [dacs vels stimtimes];
            allstims = vertcat(allstims, step3);
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 32 %special case for 230816
            allstims = [];
            %step1
            step1_file_name = exp_files{day,2};
            step1_data = TDTbin2mat(step1_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step1_time = size(step1_data.streams.Raws.data, 2)/fs;
            step1_dacs = step1_data.epocs.Us1_.data;
            step1_vels = ones(size(step1_data.epocs.Us1_.data))*(max(step1_data.streams.sRmp.data));
            step1_stim_times = step1_data.epocs.Us1_.onset;
            step1 = [step1_dacs step1_vels step1_stim_times];
            allstims = vertcat(allstims, step1);
            stim_matrices_for_each_day{day,2} = step1;
            %step2
            step2_file_name = exp_files{day,3};
            step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_time = size(step2_data.streams.Raws.data, 2)/fs;
            step2_dacs = step2_data.epocs.DacN.data;
            step2_vels = step2_data.epocs.RmpV.data;
            step2_stim_times = step2_data.epocs.RmpV.onset + step1_time;
            step2 = [step2_dacs step2_vels step2_stim_times];
            allstims = vertcat(allstims, step2);
            stim_matrices_for_each_day{day,3} = step2;
            %step2_2
            step2_2 = [];
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3 (broken into three files for this day)
            part1 = exp_files{day,5};
            part2 = exp_files{day,11};
            part3 = exp_files{day,12};
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 945.75);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 1307.44);
            part3_data = TDTbin2mat(part3, 'TYPE', {'snips', 'epocs'});
            part1_dacs = part1_data.epocs.DacN.data;
            part2_dacs = part2_data.epocs.DacN.data;
            part3_dacs = part3_data.epocs.DacN.data;
            dacs = [part1_dacs; part2_dacs; part3_dacs];
            part1_vels = part1_data.epocs.RmpV.data;
            part2_vels = part2_data.epocs.RmpV.data;
            part3_vels = part3_data.epocs.RmpV.data;
            vels = [part1_vels; part2_vels; part3_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset;
            part2_stimtimes = part2_data.epocs.RmpV.onset + 945.75;
            part3_stimtimes = part3_data.epocs.RmpV.onset + 945.75 + 1307.44;
            stimtimes = [part1_stimtimes; part2_stimtimes; part3_stimtimes];
            stimtimes = stimtimes + step1_time + step2_time;
            step3 = [dacs vels stimtimes];
            allstims = vertcat(allstims, step3);
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 48 %special case for 240130
            allstims = [];
            %step1, no step1 
            stim_matrices_for_each_day{day,2} = [];
            %step2
            step2_file_name = exp_files{day,3};
            step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_time = size(step2_data.streams.Raws.data, 2)/fs;
            step2_dacs = step2_data.epocs.DacN.data;
            step2_vels = step2_data.epocs.RmpV.data;
            step2_stim_times = step2_data.epocs.RmpV.onset ;
            step2 = [step2_dacs step2_vels step2_stim_times];
            allstims = vertcat(allstims, step2);
            stim_matrices_for_each_day{day,3} = step2;
            %step2_2
            step2_2 = [];
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3 (broken into two files for this day)
            part1 = exp_files{day,5};
            part2 = exp_files{day,11};
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 2260.77);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 911);
            part1_dacs = part1_data.epocs.DacN.data;
            part2_dacs = part2_data.epocs.DacN.data(2:end);
            dacs = [part1_dacs; part2_dacs];
            part1_vels = part1_data.epocs.RmpV.data;
            part2_vels = part2_data.epocs.RmpV.data(2:end);
            vels = [part1_vels; part2_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset;
            part2_stimtimes = part2_data.epocs.RmpV.onset(2:end) + 2260.77;
            stimtimes = [part1_stimtimes; part2_stimtimes];
            stimtimes = stimtimes + step2_time;
            step3 = [dacs vels stimtimes];
            allstims = vertcat(allstims, step3);
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 49 %special case for 240131
            allstims = [];
            %step1, no step1 
            stim_matrices_for_each_day{day,2} = [];
            %step2 (broken into two files for this day)
            part1 = exp_files{day,3};
            part2 = exp_files{day,14}; %the second step2
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 547.07);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs', 'streams'});
            part1_dacs = part1_data.epocs.DacN.data;
            part2_dacs = part2_data.epocs.DacN.data;
            dacs = [part1_dacs; part2_dacs];
            part1_vels = part1_data.epocs.RmpV.data;
            part2_vels = part2_data.epocs.RmpV.data;
            vels = [part1_vels; part2_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset;
            part2_stimtimes = part2_data.epocs.RmpV.onset + 547.07;
            stimtimes = [part1_stimtimes; part2_stimtimes];
            step2 = [dacs vels stimtimes];
            step2time = 547.07 + size(part2_data.streams.Raws.data, 2)/fs;
            allstims = vertcat(allstims, step2); 
            stim_matrices_for_each_day{day,3} = step2;
            clear part1 part2 part1_data clear part2_data part1_dacs part2_dacs dacs part1_vels part2_vels vels part1_stimtimes part2_stimtimes stimtimes
            %step2_2
            step2_2 = [];
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3 (broken into two files for this day)
            part1 = exp_files{day,5};
            part2 = exp_files{day,11};
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 557.89);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs'});
            part1_dacs = part1_data.epocs.DacN.data;
            part2_dacs = part2_data.epocs.DacN.data;
            dacs = [part1_dacs; part2_dacs];
            part1_vels = part1_data.epocs.RmpV.data;
            part2_vels = part2_data.epocs.RmpV.data;
            vels = [part1_vels; part2_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset;
            part2_stimtimes = part2_data.epocs.RmpV.onset + 557.89;
            stimtimes = [part1_stimtimes; part2_stimtimes];
            stimtimes = stimtimes + step2_time;
            step3 = [dacs vels stimtimes];
            allstims = vertcat(allstims, step3);
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 52 %special case for 240216
            allstims = [];
            %step1, no step1 
            stim_matrices_for_each_day{day,2} = [];
            %step2
            step2_file_name = exp_files{day,3};
            step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_time = size(step2_data.streams.Raws.data, 2)/fs;
            step2_dacs = step2_data.epocs.DacN.data;
            step2_vels = step2_data.epocs.RmpV.data;
            step2_stim_times = step2_data.epocs.RmpV.onset ;
            step2 = [step2_dacs step2_vels step2_stim_times];
            allstims = vertcat(allstims, step2);
            stim_matrices_for_each_day{day,3} = step2;
            %step2_2
            step2_2_file_name = exp_files{day,4};
            step2_2_data = TDTbin2mat(step2_2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_2_time = size(step2_2_data.streams.Raws.data, 2)/fs;
            step2_2_dacs = step2_2_data.epocs.DacN.data;
            step2_2_vels = step2_2_data.epocs.RmpV.data;
            step2_2_stim_times = step2_2_data.epocs.RmpV.onset + step2_time;
            step2_2 = [step2_2_dacs step2_2_vels step2_2_stim_times];
            allstims = vertcat(allstims, step2_2);
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3 (broken into two files for this day)
            part1 = exp_files{day,5};
            part2 = exp_files{day,11};
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 2160);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs'});
            part1_dacs = part1_data.epocs.DacN.data(1:1074);
            part2_dacs = part2_data.epocs.DacN.data;
            dacs = [part1_dacs; part2_dacs];
            part1_vels = part1_data.epocs.RmpV.data(1:1074);
            part2_vels = part2_data.epocs.RmpV.data;
            vels = [part1_vels; part2_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset(1:1074);
            part2_stimtimes = part2_data.epocs.RmpV.onset + 2160;
            stimtimes = [part1_stimtimes; part2_stimtimes];
            stimtimes = stimtimes + step2_time + step2_2_time;
            step3 = [dacs vels stimtimes];
            allstims = vertcat(allstims, step3);
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 53 %special case for 240220
            allstims = [];
            %step1, no step1 
            stim_matrices_for_each_day{day,2} = [];
            %step2
            step2_file_name = exp_files{day,3};
            step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_time = size(step2_data.streams.Raws.data, 2)/fs;
            step2_dacs = step2_data.epocs.DacN.data;
            step2_vels = step2_data.epocs.RmpV.data;
            step2_stim_times = step2_data.epocs.RmpV.onset ;
            step2 = [step2_dacs step2_vels step2_stim_times];
            allstims = vertcat(allstims, step2);
            stim_matrices_for_each_day{day,3} = step2;
            %step2_2
            stim_matrices_for_each_day{day,4} = [];
            %step3 (broken into two files for this day)
            part1 = exp_files{day,5};
            part2 = exp_files{day,11};
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 2062.13);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs'});
            part1_dacs = part1_data.epocs.DacN.data(1:1025);
            part2_dacs = part2_data.epocs.DacN.data;
            dacs = [part1_dacs; part2_dacs];
            part1_vels = part1_data.epocs.RmpV.data(1:1025);
            part2_vels = part2_data.epocs.RmpV.data;
            vels = [part1_vels; part2_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset(1:1025);
            part2_stimtimes = part2_data.epocs.RmpV.onset + 2062.13;
            stimtimes = [part1_stimtimes; part2_stimtimes];
            stimtimes = stimtimes + step2_time;
            step3 = [dacs vels stimtimes];
            allstims = vertcat(allstims, step3);
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
       
        if day == 57 %special case for 240320
            allstims = [];
            %step1
            step1_file_name = exp_files{day,2};
            step1_data = TDTbin2mat(step1_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step1_time = size(step1_data.streams.Raws.data, 2)/fs;
            step1_dacs = step1_data.epocs.Us1_.data;
            step1_vels = ones(size(step1_data.epocs.Us1_.data))*(max(step1_data.streams.sRmp.data));
            step1_stim_times = step1_data.epocs.Us1_.onset;
            step1 = [step1_dacs step1_vels step1_stim_times];
            allstims = vertcat(allstims, step1);
            stim_matrices_for_each_day{day,2} = step1;
            %step2
            stim_matrices_for_each_day{day,3} = [];
            %step2_2
            step2_2 = [];
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3 (broken into three files for this day)
            part1 = exp_files{day,5};
            part2 = exp_files{day,11};
            part3 = exp_files{day,12};
            part1_data = TDTbin2mat(part1, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 717.78);
            part2_data = TDTbin2mat(part2, 'TYPE', {'snips', 'epocs'}, 'T1', 0, 'T2', 2305.83);
            part3_data = TDTbin2mat(part3, 'TYPE', {'snips', 'epocs'});
            part1_dacs = part1_data.epocs.DacN.data;
            part2_dacs = part2_data.epocs.DacN.data;
            part3_dacs = part3_data.epocs.DacN.data;
            dacs = [part1_dacs; part2_dacs; part3_dacs];
            part1_vels = part1_data.epocs.RmpV.data;
            part2_vels = part2_data.epocs.RmpV.data;
            part3_vels = part3_data.epocs.RmpV.data;
            vels = [part1_vels; part2_vels; part3_vels];
            part1_stimtimes = part1_data.epocs.RmpV.onset;
            part2_stimtimes = part2_data.epocs.RmpV.onset + 717.78;
            part3_stimtimes = part3_data.epocs.RmpV.onset + 717.78 + 2305.83;
            stimtimes = [part1_stimtimes; part2_stimtimes; part3_stimtimes];
            stimtimes = stimtimes + step1_time ;
            step3 = [dacs vels stimtimes];
            allstims = vertcat(allstims, step3);
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
        if day == 61 %special case for 240629
            allstims = [];
            %step1 (ends early)
            step1_file_name = exp_files{day,2};
            step1_data = TDTbin2mat(step1_file_name, 'TYPE', {'streams', 'snips', 'epocs'}, 'T1', 0, 'T2', 913);
            step1_time = size(step1_data.streams.Raws.data, 2)/fs;
            step1_dacs = step1_data.epocs.Us1_.data;
            step1_vels = ones(size(step1_data.epocs.Us1_.data))*(max(step1_data.streams.sRmp.data));
            step1_stim_times = step1_data.epocs.Us1_.onset;
            step1 = [step1_dacs step1_vels step1_stim_times];
            allstims = vertcat(allstims, step1);
            stim_matrices_for_each_day{day,2} = step1;
            %step2
            step2_file_name = exp_files{day,3};
            step2_data = TDTbin2mat(step2_file_name, 'TYPE', {'streams', 'snips', 'epocs'});
            step2_time = size(step2_data.streams.Raws.data, 2)/fs;
            step2_dacs = step2_data.epocs.DacN.data;
            step2_vels = step2_data.epocs.RmpV.data;
            step2_stim_times = step2_data.epocs.RmpV.onset + step1_time;
            step2 = [step2_dacs step2_vels step2_stim_times];
            allstims = vertcat(allstims, step2);
            stim_matrices_for_each_day{day,3} = step2;            
            %step2_2
            step2_2 = [];
            stim_matrices_for_each_day{day,4} = step2_2;
            %step3 (one file)
            step3_file_name = exp_files{day,5};
            step3_data = TDTbin2mat(step3_file_name, 'TYPE', {'snips', 'epocs'});
            step3_dacs = step3_data.epocs.DacN.data;
            step3_vels = step3_data.epocs.RmpV.data;
            step3_stim_times = step3_data.epocs.RmpV.onset + step2_time + step1_time;
            step3 = [step3_dacs step3_vels step3_stim_times];
            allstims = vertcat(allstims, step3);            
            stim_matrices_for_each_day{day,5} = step3;
            stim_matrices_for_each_day{day,6} = allstims;
            stim_matrices_for_each_day{day,7} = exp_files{day,6};
            stim_matrices_for_each_day{day,8} = exp_files{day,7};
            stim_matrices_for_each_day{day,9} = exp_files{day,8};
            stim_matrices_for_each_day{day,10} = exp_files{day,9};
            stim_matrices_for_each_day{day,11} = exp_files{day,13};
        end
        
    end
end
end
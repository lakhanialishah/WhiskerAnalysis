%% step3 (all whiskers all velocities) script

%% VRC 

BLOCKPATH = 'F:\LAKHANI\240710\WT-240710-143303';
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'snips', 'epocs'});

velocities = unique(data.epocs.RmpV.data);
results = {};

for dac = 1:9
    dacindices = find(data.epocs.DacN.data == dac);
    dacvelocity = data.epocs.RmpV.data(dacindices);
    dactime = data.epocs.DacN.onset(dacindices);
    for i = 1:length(velocities)
        velindices = find(dacvelocity == velocities(i));
        stimtimes = dactime(velindices);
        verticalconcat = [];
        for k = 1:length(stimtimes)
            endtime = stimtimes(k) + .05;
            indices = find(data.snips.eNe2.ts > stimtimes(k) & data.snips.eNe2.ts < endtime);
            realtimestamps = data.snips.eNe2.ts(indices);
            channels = data.snips.eNe2.chan(indices);
            forplotting = [];
            forplotting = realtimestamps - stimtimes(k);
            forplotting(:,2) = channels;
            forplotting(:,3) = k; 
            verticalconcat = vertcat(verticalconcat, forplotting);
        end
        results{dac, i} = verticalconcat;  
    end    
end

totalspikes = [];
for m = 1:size(results,1)
    for j = 1:size(results,2)
        totalspikes(m,j) = length(results{m,j});
    end
end

%% total spikes figure 

figure
colors = turbo(9);
for p = 1:9
   hold on  
   plot(velocities, totalspikes(p, :), 'LineWidth', 1.5, 'Color', colors(p,:))
end
xlabel('Velocity (in V)')
ylabel('Total Spikes in 50 ms')
legend('SW1', 'SW2', 'SW3', 'SW4', 'BW', 'SW5', 'SW6', 'SW7', 'SW8'); %need to manually change where BW is

%% raster plots of one dac at different velocities 

dacofinterest = 5;
figure
for i = 1:7
    subplot(3,3, i)
    plot(results{dacofinterest,i}(:,1), results{dacofinterest,i}(:,3), 'k.')
    xlabel('Time')
    ylabel('Trial Number')
    title(velocities(i))
end

%% raster plots of all dacs at one velocity 

velofinterest = 0.3;
velindex = find(velocities >= velofinterest & velocities < velofinterest + 0.01);

figure 
for i = 1:9
    subplot(3,3,i)
    plot(results{i,velindex}(:,1), results{i,velindex}(:,3), 'k.')
    xlabel('Time')
    ylabel('Trial Number')
    title(i)
end

%% example for when this recording is broken up into multiple files 
% 220912 example 

BLOCKPATH1 = 'F:\LAKHANI\220912\WT-220912-120151';
BLOCKPATH2 = 'F:\LAKHANI\220912\WT-220912-125416';
data1 = TDTbin2mat(BLOCKPATH1, 'TYPE', {'epocs','snips', 'streams'}); %needs all the original recording files
data2 = TDTbin2mat(BLOCKPATH2, 'TYPE', {'epocs','snips', 'streams'}); 

dacs1 = data1.epocs.DacN.data;
dacs2 = data2.epocs.DacN.data(162:end);
alldacs = [dacs1; dacs2];
stimtimes1 = data1.epocs.DacN.onset;
stimtimes2 = data2.epocs.DacN.onset(162:end);
stimtimes2 = stimtimes2 + (size(data1.streams.Raws.data, 2)/24414.0625);
allstimtimes = [stimtimes1; stimtimes2];
vel1 = data1.epocs.RmpV.data;
vel2 = data2.epocs.RmpV.data(162:end);
vel = [vel1; vel2];
velocities = unique(vel);
timestamps1 = data1.snips.eNe2.ts; 
timestamps2 = data2.snips.eNe2.ts; 
timestamps2 = timestamps2 + (size(data1.streams.Raws.data, 2)/24414.0625);
allts = [timestamps1; timestamps2];
chanstamps1 = data1.snips.eNe2.chan; 
chanstamps2 = data2.snips.eNe2.chan; 
allchans = [chanstamps1; chanstamps2];

results = {};

for dac = 1:9
    dacindices = find(alldacs == dac);
    dacvelocity = vel(dacindices);
    dactime = allstimtimes(dacindices);
    for i = 1:length(velocities)
        velindices = find(dacvelocity == velocities(i));
        stimtimes = dactime(velindices);
        verticalconcat = [];
        for k = 1:length(stimtimes)
            endtime = stimtimes(k) + .05;
            indices = find(allts > stimtimes(k) & allts < endtime);
            realtimestamps = allts(indices);
            channels = allchans(indices);
            forplotting = [];
            forplotting = realtimestamps - stimtimes(k);
            forplotting(:,2) = channels;
            forplotting(:,3) = k; 
            verticalconcat = vertcat(verticalconcat, forplotting);
        end
        results{dac, i} = verticalconcat;  
    end    
end

totalspikes = [];
for m = 1:size(results,1)
    for j = 1:size(results,2)
        totalspikes(m,j) = length(results{m,j});
    end
end

figure
colors = turbo(9);
for p = 1:9
   hold on  
   plot(velocities, totalspikes(p, :), 'LineWidth', 1.5, 'Color', colors(p,:))
end
xlabel('Velocity (in V)')
ylabel('Total Spikes in 50 ms')
legend('SW1', 'SW2', 'SW3', 'SW4', 'PW', 'SW5', 'SW6', 'SW7', 'SW8'); %need to manually change where PW is

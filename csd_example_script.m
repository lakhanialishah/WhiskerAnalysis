%% csd with MUA BW one velocity only 

PWdac = 5;
excludedchans = [37, 38, 39, 40, 41, 42, 56, 57];
BLOCKPATH = 'F:\LAKHANI\240705\WT-240705-112721'; %step1 
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'snips', 'epocs', 'streams'});

%find a stim and calculate what lfp points I need to take 
indices = find(data.epocs.Us1_.data == PWdac);
stimtimes = data.epocs.RAMP.onset(indices);
totaltime = size(data.streams.LFP1.data,2)/data.streams.LFP1.fs;
eachstep = totaltime/size(data.streams.LFP1.data,2);
timevector = 0:eachstep:totaltime;

%to interpolate between the channels after excluding non-working chans 
newlfp = double(data.streams.LFP1.data);
% newlfp(excludedchans,:) = NaN;
% for chan = 1:size(newlfp,1)
%     if chan == 40
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.75 (newlfp(42,i))*.25];
%             newlfp(chan,i) = mean(toavg);
%         end
%     elseif chan == 41
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.25 (newlfp(42,i))*.75];
%             newlfp(chan,i) = mean(toavg);
%         end
%     else
%         for i = 1:size(newlfp,2)
%                 if isnan(newlfp(chan,i)) == 1
%                     toavg = [newlfp(chan+1,i) newlfp(chan-1,i)];
%                     newlfp(chan,i) = mean(toavg);
%                 end
%         end
%     end
% end

%to get stim associated lfp 
allchans = [];
for i = 1:length(stimtimes)
    firstidx = find(timevector>stimtimes(i),1);
    lastidx = firstidx + round(0.05/eachstep); %to get the .1 second after the stim
    interval = firstidx:lastidx;
    for ch = 1:size(newlfp,1)
        y = newlfp(ch,interval); 
        allchans(i,ch,:) = y;
    end
end
allstimavg = squeeze(mean(allchans,1));

%csd
param = {};
param.iCSD_method = 'delta_source';
param.h = 20; %um for Cambridge probe 
param.sigma_c = .3; %what the paper has 
param.R = 21.65; %um
param.dist_off_center = 0;
iCSD_method = param.iCSD_method;
r_LFP = [];

[iCSD, F_mat_inv, r_LFP] = CSD(allstimavg, param);
figure
subplot(2,2,1)
imagesc(iCSD)
title('BW one velocity')

%% csd with MUA BW all velocities averaged 

excludedchans = [37, 38, 39, 40, 41, 42, 56, 57];
BLOCKPATH = 'F:\LAKHANI\240705\WT-240705-115002'; %step2
data = TDTbin2mat(BLOCKPATH, 'TYPE', {'snips', 'epocs', 'streams'});

N = data.streams.Raws.fs/data.streams.LFP3.fs; 
for chan = 1:64 
    rawdecdata(chan,:) = data.streams.Raws.data(chan, 1:N:end); %decimate raw channel by N
end
data.streams.RAWdec.data = rawdecdata; 
data.streams.RAWdec.fs = data.streams.LFP3.fs; 
data.streams.RAWdec.name = 'RawDec'; 

stimtimes = data.epocs.DacN.onset;
totaltime = size(data.streams.RAWdec.data,2)/data.streams.RAWdec.fs;
eachstep = totaltime/size(data.streams.RAWdec.data,2);
timevector = 0:eachstep:totaltime;

%to interpolate between the channels after excluding non-working chans 
newlfp = double(data.streams.RAWdec.data);
% newlfp(excludedchans,:) = NaN;
% for chan = 1:size(newlfp,1)
%     if chan == 40
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.75 (newlfp(42,i))*.25];
%             newlfp(chan,i) = mean(toavg);
%         end
%     elseif chan == 41
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.25 (newlfp(42,i))*.75];
%             newlfp(chan,i) = mean(toavg);
%         end
%     else
%         for i = 1:size(newlfp,2)
%                 if isnan(newlfp(chan,i)) == 1
%                     toavg = [newlfp(chan+1,i) newlfp(chan-1,i)];
%                     newlfp(chan,i) = mean(toavg);
%                 end
%         end
%     end
% end

%to get stim associated lfp 
allchans = [];
for i = 1:length(stimtimes)
    firstidx = find(timevector>stimtimes(i),1);
    lastidx = firstidx + round(0.05/eachstep); %to get the .1 second after the stim
    interval = firstidx:lastidx;
    for ch = 1:size(newlfp,1)
        y = newlfp(ch,interval); 
        allchans(i,ch,:) = y;
    end
end
allstimavg = squeeze(mean(allchans,1));
% plot(allstimavg')

%csd
param = {};
param.iCSD_method = 'delta_source';
param.h = 20; %um for Cambridge probe 
param.sigma_c = .3; %what the paper has 
param.R = 21.65; %um
param.dist_off_center = 0;
iCSD_method = param.iCSD_method;
r_LFP = [];

[iCSD, F_mat_inv, r_LFP] = CSD(allstimavg, param);
subplot(2,2,2)
imagesc(iCSD)
title('BW all velocities')

%% csd with BW fastest 3 velocities averaged 

fastestvels = find(data.epocs.RmpV.data > 0.3);
stimtimes = data.epocs.DacN.onset(fastestvels);
totaltime = size(data.streams.RAWdec.data,2)/data.streams.RAWdec.fs;
eachstep = totaltime/size(data.streams.RAWdec.data,2);
timevector = 0:eachstep:totaltime;

%to interpolate between the channels after excluding non-working chans 
newlfp = double(data.streams.RAWdec.data);
% newlfp(excludedchans,:) = NaN;
% for chan = 1:size(newlfp,1)
%     if chan == 40
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.75 (newlfp(42,i))*.25];
%             newlfp(chan,i) = mean(toavg);
%         end
%     elseif chan == 41
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.25 (newlfp(42,i))*.75];
%             newlfp(chan,i) = mean(toavg);
%         end
%     else
%         for i = 1:size(newlfp,2)
%                 if isnan(newlfp(chan,i)) == 1
%                     toavg = [newlfp(chan+1,i) newlfp(chan-1,i)];
%                     newlfp(chan,i) = mean(toavg);
%                 end
%         end
%     end
% end

%to get stim associated lfp 
allchans = [];
for i = 1:length(stimtimes)
    firstidx = find(timevector>stimtimes(i),1);
    lastidx = firstidx + round(0.05/eachstep); %to get the __ second after the stim
    interval = firstidx:lastidx;
    for ch = 1:size(newlfp,1)
        y = newlfp(ch,interval); 
        allchans(i,ch,:) = y;
    end
end
allstimavg = squeeze(mean(allchans,1));

%csd
param = {};
param.iCSD_method = 'delta_source';
param.h = 20; %um for Cambridge probe 
param.sigma_c = .3; %what the paper has 
param.R = 21.65; %um
param.dist_off_center = 0;
iCSD_method = param.iCSD_method;
r_LFP = [];

[iCSD, F_mat_inv, r_LFP] = CSD(allstimavg, param);
subplot(2,2,3)
imagesc(iCSD)
title('BW fastest 3 velocities')

%% csd with BW weakest 3 velocities averaged 

slowestvels = find(data.epocs.RmpV.data > 0 & data.epocs.RmpV.data < 0.3);
stimtimes = data.epocs.DacN.onset(slowestvels);
totaltime = size(data.streams.RAWdec.data,2)/data.streams.RAWdec.fs;
eachstep = totaltime/size(data.streams.RAWdec.data,2);
timevector = 0:eachstep:totaltime;

%to interpolate between the channels after excluding non-working chans 
newlfp = double(data.streams.RAWdec.data);
% newlfp(excludedchans,:) = NaN;
% for chan = 1:size(newlfp,1)
%     if chan == 40
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.75 (newlfp(42,i))*.25];
%             newlfp(chan,i) = mean(toavg);
%         end
%     elseif chan == 41
%         for i = 1:size(newlfp,2)
%             toavg = [(newlfp(39,i))*.25 (newlfp(42,i))*.75];
%             newlfp(chan,i) = mean(toavg);
%         end
%     else
%         for i = 1:size(newlfp,2)
%                 if isnan(newlfp(chan,i)) == 1
%                     toavg = [newlfp(chan+1,i) newlfp(chan-1,i)];
%                     newlfp(chan,i) = mean(toavg);
%                 end
%         end
%     end
% end

%to get stim associated lfp 
allchans = [];
for i = 1:length(stimtimes)
    firstidx = find(timevector>stimtimes(i),1);
    lastidx = firstidx + round(0.05/eachstep); %to get the __ second after the stim
    interval = firstidx:lastidx;
    for ch = 1:size(newlfp,1)
        y = newlfp(ch,interval); 
        allchans(i,ch,:) = y;
    end
end
allstimavg = squeeze(mean(allchans,1));
% plot(allstimavg')

%csd
param = {};
param.iCSD_method = 'delta_source';
param.h = 20; %um for Cambridge probe 
param.sigma_c = .3; %what the paper has 
param.R = 21.65; %um
param.dist_off_center = 0;
iCSD_method = param.iCSD_method;
r_LFP = [];

[iCSD, F_mat_inv, r_LFP] = CSD(allstimavg, param);
subplot(2,2,4)
imagesc(iCSD)
title('BW slowest 3 velocities')
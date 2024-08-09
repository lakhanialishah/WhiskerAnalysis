%% 240705 post phy analysis to extract single units
% needs https://github.com/cortex-lab/spikes/ 

expday = 240705;
genotype = 'KO';
WDdays = 0; %how many days of whisker deprivation
age = 22; 
BW = 'D2';
excludedchans = [37, 38, 39, 40, 41, 42, 56, 57];

%% set up 

expt_name = '240705';
dat_file_name = 'file_240705.dat';
nCh = 64; 
fs = 24414.0625;
wfWin = [-24 48]; %window of spike waveforms to extract,in samples; this is ~ -1 ms to +2 ms around each spike
dataType = 'int16';
KS2_dir = strcat('D:\Data\', expt_name);
phy_output_dir = strcat('D:\Data\', ...
    expt_name);
processed_phy_output_dir = strcat('D:\Data\', expt_name);

spikeTimes = readNPY(fullfile(phy_output_dir, 'spike_times.npy'));
spikeClusters = readNPY(fullfile(phy_output_dir, 'spike_clusters.npy'));
gwfparams.expt_name = expt_name;
gwfparams.KS2_dir = KS2_dir;
gwfparams.dat_file_name = dat_file_name;
gwfparams.processed_phy_output_dir = processed_phy_output_dir ;
gwfparams.dataType = dataType;
gwfparams.nCh = nCh;
gwfparams.wfWin = wfWin;
gwfparams.spikeTimes = spikeTimes;
gwfparams.spikeClusters = spikeClusters;
gwfparams.fs = fs;
gwfparams.nWf = 2000;     
gwfparams.dataDir = strcat('D:\Data\', expt_name);
gwfparams.fileName = 'file_240705.dat';

data_phy = loadKSdir(KS2_dir);
gwfparams.dataDir = KS2_dir;

%% calculating and saving what I want for each unit

refractoryPeriod = 1; %in ms
goodunits = data_phy.cids(data_phy.cgs == 2); %only getting the good units from phy, 2 means good, 1 means multiunit 
totalClus = numel(goodunits); 
stats = {};
processedresultspath = ''; 

% to figure out what the true channel is based on phy  
allchans = (1:64)';
probechans = allchans;
if ~isempty(excludedchans)
    allchans(excludedchans) = NaN;
    todelete = find(isnan(allchans(:,1)) == 1);
    probechans(todelete,:) = [];
end
phychans = (1:length(probechans))'; %1 indexed because the getwaveforms function is making me 1 index 
chanmatrix = [phychans, probechans];

all_spike_times = double(readNPY(fullfile('D:\Data\', expt_name, 'spike_times.npy')));
all_spike_times_seconds = double(readNPY(fullfile('D:\Data\', expt_name, 'spike_times.npy')))/24414.0625; %divide by sampling rate to convert to seconds
all_spike_clusters = double(readNPY(fullfile('D:\Data\', expt_name, 'spike_clusters.npy')));
all_amplitudes = double(readNPY(fullfile('D:\Data\', expt_name, 'amplitudes.npy')));

for clu_Num = 1:totalClus
    clu_ID = goodunits(clu_Num); %only going through what I sorted as good in phy 
    clu_spike_times = all_spike_times_seconds(all_spike_clusters == clu_ID);
    clu_amplitudes = all_amplitudes(all_spike_clusters == clu_ID);
    clu_info.expday = expday; 
    clu_info.cluID = clu_ID; 
    clu_info.spikeTimes = clu_spike_times;
    clu_info.amplitudes = clu_amplitudes; 
    
    %using the getwaveform function 
    gwfparams.spikeTimes = all_spike_times(all_spike_clusters == clu_ID);
    indices = find(all_spike_clusters == clu_ID); 
    gwfparams.spikeClusters = all_spike_clusters(indices);
    [waveforms, unitchan, mwf, sdwf] = get_waveforms_and_AP_times_from_phy_output(gwfparams); %need to be in the folder NewScriptsCaleb
%     plot(mwf(:,unitchan), 'k', 'LineWidth', 2)
    minimum = min(mwf(:, unitchan));
    minindex = find(mwf(:,unitchan) == minimum);
    wfpostmin = mwf(minindex+1:end, unitchan);
    maximum = max(wfpostmin);
    maxindex = find(mwf(:, unitchan) == maximum);
    samples = maxindex - minindex;
    tr2pk = (samples/gwfparams.fs)*1000;
    amp = maximum - minimum; 
    idx = find(chanmatrix(:,1) == unitchan); 
    clu_info.mwf = mwf(:,unitchan);
    clu_info.sdwf = sdwf(:, unitchan);
    clu_info.phychan = unitchan; % channel that the unit is largest on in phy 
    clu_info.truechan = chanmatrix(idx,2); %channel that the unit is truly largest on in terms of the probe
    clu_info.tr2pk = tr2pk; %trough to peak time for the unit
    clu_info.mean_amp = mean(clu_amplitudes); %amplitude of the waveform by taking the mean of all the amplitudes
    
    processed_results_filename = [processedresultspath, '\', 'cluster' num2str(clu_ID), '_', num2str(expday) '.mat'];
    save(processed_results_filename, 'clu_info');   
    
end

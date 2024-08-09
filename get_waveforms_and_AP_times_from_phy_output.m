function [waveForms, unitchan, mwf, sdwf] = get_waveforms_and_AP_times_from_phy_output(gwfparams, data_phy)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% This is a slight modification of the function getWaveForms,
% originally contributed by C. Schoonover and A. Fink.
% See https://github.com/cortex-lab/spikes/blob/master/analysis/getWaveForms.m
% for original code.
% Written by Caleb Wright
%
% % EXAMPLE INPUT
% gwfparams.KS2_dir = '/path/to/data/';    % KiloSort output folder
% gwfparams.processed_phy_output_dir = '/path/to/processed/results'
% gwfparams.fileName = 'AllCombinedFiles.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 64;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-24 48];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%
% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.KS2_dir,gwfparams.dat_file_name);          
filenamestruct = dir(fileName);
nCh = gwfparams.nCh;
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = gwfparams.wfWin(2) - gwfparams.wfWin(1);
mmf = memmapfile(fileName, 'Format', {'int16', [double(nCh) double(nSamp)], 'x'});
num_tsteps = size(mmf.Data.x, 2);
chMap = readNPY(fullfile(gwfparams.KS2_dir, 'channel_map.npy'))+1;
%chMap is a (one-indexed) 1 x num_good_channels column vector of channels 
%determined to be good by KiloSort2.  Each entry is a good channel number.  
%Bad channels are skipped over.  For example, if channel 2 is bad, then
%chMap = [1 3 4 ...].  Note that channel numbers is a subset of the (
%one-indexed) chanMap used by KiloSort2, which in turn reflects the 
%order in which data was streamed
%to the .dat file.  NCW always streamed data in order of channel number.
%For the TDT system, this corresponded to probe geometry.  For the Cerebus
%system, this did not.  For instance, for the TDT system, chanMap = [1 2 3
%4 ...].  For the Cerebus system, chanMap = [1 32 2 31 ...].

nChInMap = numel(chMap); %number of good channels only

%design high-pass filter
fs = gwfparams.fs;
dt = 1/fs;
tspc_filt = 25;
hpFilt = designfilt('highpassiir', 'StopbandFrequency', 100, ...
    'PassbandFrequency', 500, 'StopbandAttenuation', 60, ...
    'PassbandRipple', 1, 'SampleRate', fs, 'DesignMethod', 'butter');

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters); %unique(data_phy.cids);
numUnits = size(unitIDs,1);

for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    processed_results_filename = [gwfparams.processed_phy_output_dir, ...
        'cluster' num2str(curUnitID) '.mat'];
    %figure
    %waveFormsMean = nan(nChInMap,wfNSamples);

    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    fprintf('\nprocessing %s', gwfparams.expt_name);
    fprintf(' cluster %i', curUnitID);
    fprintf(' (%i spikes)\n', curUnitnSpikes)

    if curUnitnSpikes > 1
%         if curUnitnSpikes > 2000
%             %Use random subset of all spikes to calculate waveforms (to
%             %speed things up), but
%             %get spike times for ALL spikes
%             AP_indices_to_use = randperm(curUnitnSpikes, 2000);
%             curUnitnSpikes = 2000;
%             fprintf('***Using 2000 randomly-selected spike times to calculate waveforms (to speed things up)***\n')
%             curSpikeTimes_ = curSpikeTimes(AP_indices_to_use); %work w/a smaller array
%         else
%             AP_indices_to_use = 1:curUnitnSpikes;
%             curSpikeTimes_ = curSpikeTimes;
%         end
        curSpikeTimes_ = curSpikeTimes;
        waveForms = nan(curUnitnSpikes,nChInMap,wfNSamples);

        for r = 1:curUnitnSpikes
            curSpikeTime = int64(curSpikeTimes_(r));

            %make sure window surrounding each spike time falls w/in
            %length of recording
            if curSpikeTime + gwfparams.wfWin(2)+tspc_filt - 1 < ...
                    num_tsteps && curSpikeTime...
                    +gwfparams.wfWin(1)-tspc_filt > 0

                tmpWf = mmf.Data.x(1:nCh,curSpikeTime...
                    +gwfparams.wfWin(1)-tspc_filt:curSpikeTime...
                    +gwfparams.wfWin(2)+tspc_filt - 1);

                %numerical order of tmpWf entries corresponds to order
                %in which data was written to .dat file.  For Cerebus
                %system, this does NOT correspond to probe geometry.

                tmpfilt = nan(nCh,wfNSamples); %num channels in dat file x wfNSamples

                for j = 1:size(tmpWf, 1)
                    %cycle over entries in tmpWf
                    %So tmpfilt will be in same order as tmpWf, which
                    %is order in which data was written to .dat file.
                    filt_trace = filtfilt(hpFilt,...
                    double(tmpWf(j, :))); 
                    tmpfilt(j, :) = filt_trace(...
                        tspc_filt+1:tspc_filt+wfNSamples);
                end

                %keep only good channels
                waveForms(r,:,:) = tmpfilt(chMap,:);

                %so waveForms is now in GEOMETRIC order, since chMap is
                %in geometric order, and we're using it to index into
                %tmpfilt.

            end
        end

        fprintf('Completed waveform calculation.\n')

        mwf = squeeze(nanmean(waveForms(:,:,:),1))'; 
        %num_tsteps x num_good_chans.  Note that mwf is in geometric
        %order (bottom-to-top), since waveForms is in geometric order.

        sdwf = squeeze(nanstd(waveForms(:,:,:),1))';            
        wfstruct.mwf = mwf;
        wfstruct.sdwf = sdwf;

        %CONVERT SPIKE TIMES TO SECONDS! (To FOUR decimal places.)
        spts = double(curSpikeTimes).*dt; 

        stats.Vppm = max(mwf,[],1)-min(mwf,[],1); %amplitude of mean waveform, on each good channel
        stats.SNRwf = stats.Vppm./mean(sdwf,1);
        unitchan = find(stats.Vppm==max(stats.Vppm)); %(one-indexed)
        %index of 'biggest' channel, according to probe geometry
        %(bottom-to-top)
        stats.unitchan = unitchan; %channel index (of good chans only) on which mean waveform is largest (one-indexed)
        stats.unitchanspread = sum(stats.SNRwf>1);  %SNR=1 chosen arbitrarily. May need updating!
        isispts = spts(2:end)-spts(1:end-1);
        [n,x] = histc(isispts,0:0.001:1);
        stats.ISI1msbin = n(1)/sum(n);
        stats.ISI2msbin = n(2)/sum(n);   
        
       
%         clearvars wfstruct stats waveForms mwf sdwf
        %waveFormsMean(curUnitInd,:,:) = meanwvfrm;
        disp(['Completed ' int2str(curUnitInd) ' units of ' ...
            int2str(numUnits) ' (' int2str(size(curSpikeTimes,1)) ...
            ' spikes)' '.']);
        %end
    end
end





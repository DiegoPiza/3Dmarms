%{
Max Interval Method to Burst Detection (Fixed threshold-based method)
Bursts are defined using five fixed threshold parameters (max ISI at start
of burst, max ISI in burst, min burst duration, min IBI, and min number of
spikes in burst). The values of these parameters are chosen a priori and
any series of spikes that satisfy these thresholds is classified as a
burst.

KArtik Pradeepan June 2020
Edits - BC
%}
function [allBurstData] = maxIntervalHCstrict(spikeTrain, SAMPLING_RATE,minISI)
% Bursts Parameters
max_begin_ISI = minISI; % seconds
max_end_ISI = minISI; % seconds
min_IBI = minISI+ 0.001; % seconds
min_burst_duration = 0.0015; % seconds
min_spikes_in_burst = 3; % spikes
% Spike train parameters
spikeBin = find(spikeTrain)/SAMPLING_RATE;
nSpikes = length(spikeBin);
% Create output storage
allBurstData = {};
inBurst = false;
burstNum = 0;
currentBurst = [];
%% Phase 1 - Burst Detection
% Here a burst is defined as starting when two consecutive spikes have an
% ISI less than max_begin_ISI apart. The end of the burst is given when two
% spikes have an ISI greater than max_end_ISI.
% Find ISIs closer than max_begin_ISI and end with max_end_ISI.
% The last spike of the previous burst will be used to calculate the IBI.
% For the first burst, there is no previous IBI.
n = 2;
while (n <= nSpikes)
    ISI = spikeBin(n) - spikeBin(n-1); % Calculate ISI
    if inBurst == 1 % currently in burst
        if ISI > max_end_ISI % end the burst
            currentBurst = [currentBurst, spikeBin(n-1)]; %store spike in burst
            burstNum = burstNum + 1;
            allBurstData{burstNum} = currentBurst; % store burst data
            currentBurst = []; % reset for new burst
            inBurst = 0; % no longer in burst
        else
            currentBurst = [currentBurst, spikeBin(n-1)];
        end
    else % currently not in burst
        if ISI < max_begin_ISI % possibly found start of new burst
            currentBurst = [currentBurst, spikeBin(n-1)];
            inBurst = 1;
        end
    end
    n = n + 1; % move to next spike
end
if isempty(allBurstData)
    return
end
% Calculate IBIs
IBI = [];
for b=2:burstNum
    prevBurstEnd = allBurstData{b-1}(end);
    currBurstBeg = allBurstData{b}(1);
    IBI = [IBI, currBurstBeg - prevBurstEnd];
end
%% Phase 2 - Merging of Bursts
% Here we see if any pair of bursts have an IBI less than min_IBI; if so,
% we then merge the bursts. We specifically need to check when say three
% bursts are merged into one.
% tic
allBurstTemp = cell(1,sum(IBI>min_IBI)+1);
allBurstTemp{1} = allBurstData{1};
newInd=1;
oldInd =1;
while newInd < size(allBurstTemp,2)
    nextBurst = allBurstData{oldInd+1};
    currBurst = allBurstTemp{newInd};
    if IBI(oldInd) < min_IBI %IBI is too short to be separate bursts
        %combine this burst with the next burst
        allBurstTemp{newInd} = [currBurst,nextBurst];
%         disp ("Merged burstNum: " + string(i));
    else
        newInd = newInd+1;

        allBurstTemp{newInd} = nextBurst;

    end
     oldInd = oldInd+1;
end
allBurstData = allBurstTemp;
% toc
%% Phase 3 - Quality Control
% Remove small bursts less than min_bursts_duration or having too few
% spikes less than min_spikes_in_bursts. In this phase we have the
% possibility of deleting all spikes.
tooShort = 0;
burstNum = length(allBurstData);
if burstNum > 1
    for b=1:burstNum
        currentBurst = allBurstData{b};
        if length(currentBurst) < min_spikes_in_burst
            % disp ("Deleted burstNum: " + string(b) + " because not enough spikes.");
            currentBurst = [];
        elseif (currentBurst(end)-currentBurst(1)) < min_burst_duration
            % disp ("Deleted burstNum: " + string(b) + " because duration too short.");
%             (currentBurst(end)-currentBurst(1));
            currentBurst = [];
            tooShort = tooShort+1;
        end
        allBurstData{b} = currentBurst;
    end
end

allBurstData = allBurstData(~cellfun('isempty',allBurstData));

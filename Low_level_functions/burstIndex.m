function [outUnit] = burstIndex(unit,raster)
%get the input data, and then run all the periods
% outUnit is slope and then intercept of the bursting period, and then the same
% for the late period. OutFRs is firing rate, and then burst fraction
%https://pubmed.ncbi.nlm.nih.gov/9246432/
%tasks = fieldnames(Data.(unit));
countLimit = 40; %ms
outUnit = table({unit},'VariableNames',"unitIDs");  
%outFRs = table({unit},'VariableNames',{'unitIDs'});
%     currTask = Data.(unit).(tasks{task});
%     periods = fieldnames(currTask);
%     for period = 1:length(periods)
%         uData = currTask.(periods{period});
%         if isempty(uData)
%             outUnit.([tasks{task} periods{period } 'Meas']) = zeros(1,1);
%             outFRs.([tasks{task} periods{period } 'Meas']) = [nan(1,2)];
%         elseif length(vertcat(uData.RastInds{:})) < 50
%             outUnit.([tasks{task} periods{period } 'Meas']) = zeros(1,1);
%             outFRs.([tasks{task} periods{period } 'Meas']) = [nan(1,2)];
%         else
            FRs = sum(raster)/(length(raster)/1000)/1000; %in ms
            %BF = nansum(vertcat(uData.ISI{:})<.021)./length(vertcat(uData.ISI{:}));
            ISIs = diff(find(raster));%ceil(vertcat(uData.ISI{:})*1000);
            shortISIs = ISIs(ISIs<countLimit);
            allCounts = histcounts(shortISIs,'BinEdges',0:1:countLimit);
            allProportions = allCounts/sum(allCounts);
            k = zeros(1,countLimit-1);
            for t = 2:countLimit
                k(t-1) = FRs*exp(-FRs * (t*2));
            end
            proportions = k / sum(k);
            outUnit.BurstInd = getBurstIndex(allProportions,proportions);
            outUnit.FRs=FRs*1000;
            %outFRs.([tasks{task} periods{period } 'Meas']) = [FRs, BF]; 
end


function [difference]=getBurstIndex(proportions,compProportions)
startSearch = 3;
endSearch = 20;
[BurstPeak,burstInd ]= max(proportions(startSearch:endSearch));
compPeak = compProportions(startSearch+burstInd);
burstData = sum(proportions(1:burstInd+3));
compData = sum(compProportions(1:3+burstInd));
difference(:,1) = (BurstPeak - compPeak)/(BurstPeak+compPeak);
difference(:,2) = (burstData - compData)/(burstData + compData);
% sub20data = sum(proportions(startSearch:endSearch));
% sub20comp = sum(compProportions(startSearch:endSearch));
% difference(:,3) = (sub20data - sub20comp)/(sub20data + sub20comp);
burstIndTrue = burstInd+startSearch-1;
sub20data = sum(proportions(max(1,burstIndTrue-3):burstIndTrue+3));
sub20comp = sum(compProportions(max(1,burstIndTrue-3):burstIndTrue+3));
difference(:,3) = (sub20data - sub20comp)/(sub20data + sub20comp);
difference(:,4) = (sum(proportions(2:7)) - sum(compProportions(2:7)))...
    /(sum(proportions(2:7)) + sum(compProportions(2:7)));
end

function SIs = Sp_Inf_Cont_2022_update(OccupancyMap, frs, count)
    % Computes spatial information scores (SIs) based on occupancy and count data
    
    % Input:
    % OccupancyMap - Occupancy map data
    % frs - Firing rates data
    % count - Count data
    
    % Output:
    % SIs - Spatial information scores
    
    % Set unoccupied bins to NaN
    OccupancyMap(OccupancyMap == 0) = nan;
    
    % Calculate occupancy sum and occupancy proportion
    occupancy_sum = nansum(OccupancyMap);
    occupancy_proportion = OccupancyMap / occupancy_sum;
    
    % Replace NaN values in count with NaN in the same positions as unoccupied bins
    count(isnan(OccupancyMap)) = nan;
    
    % Calculate lambda_mean (mean firing rate)
    lambda_mean = nansum(count) / nansum(OccupancyMap);
    
    % Initialize SIs array and counters
    SIs = nan(numel(OccupancyMap), 1);
    skipped = 0;
    computed = 0;

    % Calculate SIs for each bin
    for i = 1:numel(OccupancyMap)
        if isnan(occupancy_proportion(i))
            skipped = skipped + 1;
            continue
        else
            computed = computed + 1;
            % Calculate spatial information score (SI)
            SIs(i) = occupancy_proportion(i) * (frs(i) / lambda_mean) * log2(frs(i) / lambda_mean + eps);
            
            % Check for NaN values and errors
            if isnan(SIs(i))
                error([num2str(i) ' computed a nan SI']);
            end
        end
    end
    
    % Check for negative SIs sum
    if nansum(SIs) < 0
        error([num2str(i) ' computed a negative SIs']);
    end
    
    % Calculate mean, median, and total SIs
    si_mean = nanmean(nanmean(SIs));
    si_median = nanmedian(nanmedian(SIs));
    si_total = nansum(nansum(SIs));
    
    % Plot histogram if needed
    % histogram(SIs)
end

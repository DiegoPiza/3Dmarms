function [fix, start_end_sac, peaks_sac, start_end_fix, peaks_fix, locss, amplitudesac] = fixsacepoch(TrackingData)
%FIXSACEPOCH Process tracking data to determine saccade and fixation epochs
%   This function calculates saccade and fixation epochs based on 
%   angular velocity of the head. It returns the periods and peaks for 
%   both saccades and fixations.

    % Calculate angular velocity
    ang_vel = HeadSpeed(TrackingData, 1/60, 0);
    sr = 60; % Sampling rate

    % Find saccade epochs
    epoch_ang = find_epochs(ang_vel, 200, true, 4);

    % Initialize variables
    [onset, term, pks, locs, amplitude, lengtha] = deal(nan([length(epoch_ang) 1]));

    % Extract parameters for saccade epochs
    for i = 1:length(epoch_ang)
        onset(i) = epoch_ang(i).start_time;

        if ~isempty(epoch_ang(i).end_time)
            term(i) = epoch_ang(i).end_time;
        else
            term(i) = length(ang_vel);
        end

        [pks(i), locs(i)] = max(ang_vel(onset(i): term(i)));
        locs(i) = onset(i) + locs(i) - 1;
        amplitude(i) = sum(ang_vel(onset(i):term(i)) * 0.0167);

        if ~isempty(epoch_ang(i).end_time)
            lengtha(i) = epoch_ang(i).epoch_length;
        else
            lengtha(i) = term(i) - onset(i);
        end
    end

    % Filter out movements with less than 10 deg amplitude
    valid = amplitude >= 10;
    [onset, term, pks, locs, amplitude, lengtha] = arrayfun(@(x) x(valid), {onset, term, pks, locs, amplitude, lengtha}, 'UniformOutput', false);

    % Construct saccade-related outputs
    sac = horzcat(onset{:});
    start_end_sac = [onset{:}; term{:}];
    peaks_sac = pks{:};
    locss = locs{:};
    amplitudesac = amplitude{:};

    % Find fixation epochs
    epoch_ang = find_epochs(ang_vel, 200, false, 1);

    % Re-initialize variables
    [onset, term, pks, locs, amplitude, lengtha] = deal(nan([length(epoch_ang) 1]));

    % Extract parameters for fixation epochs
    for i = 1:length(epoch_ang)
        onset(i) = epoch_ang(i).start_time;

        if ~isempty(epoch_ang(i).end_time)
            term(i) = epoch_ang(i).end_time;
        else
            term(i) = length(ang_vel);
        end

        [pks(i), locs(i)] = max(ang_vel(onset(i): term(i)));
        locs(i) = onset(i) + locs(i) - 1;
        amplitude(i) = sum(ang_vel(onset(i):term(i)) * 0.0167);

        if ~isempty(epoch_ang(i).end_time)
            lengtha(i) = epoch_ang(i).epoch_length;
        else
            lengtha(i) = term(i) - onset(i);
        end
    end

    % Filter out fixations with less than 200ms in duration
    valid = lengtha >= 14;
    [onset, term, pks, locs, amplitude, lengtha] = arrayfun(@(x) x(valid), {onset, term, pks, locs, amplitude, lengtha}, 'UniformOutput', false);

    % Construct fixation-related outputs
    fix = horzcat(onset{:});
    start_end_fix = [onset{:}; term{:}];
    peaks_fix = pks{:};

end

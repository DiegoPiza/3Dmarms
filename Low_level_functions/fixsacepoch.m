function [fix, start_end_sac, peaks_sac, start_end_fix, peaks_fix, locss, amplitudesac] = fixsacepoch(TrackingData)
%FIXSACEPOCH Process tracking data to determine rapid head movements and fixation epochs
%  fix= fixation timestamps
% start_end_sac= onset and offset timestamps of head movements
% peaks_sac= peak head movement velocity 
% start_end_fix= onset and offset timestamps of head fixations
% locss= timestamps for the peak velocity of head movements
% amplitudesac= amplitude in degrees of head movement

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
    valid = lengtha >= floor(sr*0.2);
    [onset, term, pks, locs, amplitude, lengtha] = arrayfun(@(x) x(valid), {onset, term, pks, locs, amplitude, lengtha}, 'UniformOutput', false);

    % Construct fixation-related outputs
    fix = horzcat(onset{:});
    start_end_fix = [onset{:}; term{:}];
    peaks_fix = pks{:};

end

function [h_vel, v_vel, r_vel] = HeadVelocity(TrackingData, sr, gaus_kernel)
    % Computes horizontal, vertical, and roll head movement rotation velocity in deg/sec
    
    % Input:
    % TrackingData - Data containing tracking information
    % sr - Sampling rate
    % gaus_kernel - Flag indicating whether to apply Gaussian kernel smoothing
    
    % Output:
    % h_vel - Horizontal head movement rotation velocity (deg/sec)
    % v_vel - Vertical head movement rotation velocity (deg/sec)
    % r_vel - Roll head movement rotation velocity (deg/sec)
    
    q1 = TrackingData.WQuat;
    q2 = TrackingData.iXQuat;
    q3 = TrackingData.jYQuat;
    q4 = TrackingData.kZQuat;

    [x, y, z] = quat2angle([q1, q2, q3, q4]); % x is horizontal, y is roll, z is vertical head rotations

    h_vel = nan(size(z));
    v_vel = nan(size(z));
    r_vel = nan(size(z));

    h_vel(2:end) = rad2deg(angdiff(x(1:end-1), x(2:end))) / (1/sr);
    v_vel(2:end) = rad2deg(angdiff(z(1:end-1), z(2:end))) / (1/sr);
    r_vel(2:end) = rad2deg(angdiff(y(1:end-1), y(2:end))) / (1/sr);

    % Remove extreme values
    r_vel(h_vel > 2000, :) = nan;
    r_vel(v_vel > 2000, :) = nan;
    r_vel(h_vel < -2000, :) = nan;
    r_vel(v_vel < -2000, :) = nan;
    h_vel(h_vel > 2000, :) = nan;
    v_vel(v_vel > 2000, :) = nan;
    h_vel(h_vel < -2000, :) = nan;
    v_vel(v_vel < -2000, :) = nan;

    % Apply Gaussian smoothing if specified
    if gaus_kernel
        inv = isnan(h_vel);
        h_vel = fillmissing(h_vel, 'linear', 'EndValues', 'none');
        % Gaussian kernel
        sigma = floor(0.25 * sr); % 50 ms window
        L = sigma * 12; % Length of 6 sigma each side
        alpha = ((L - 1) / sigma) / 2;
        h = 1 / (sqrt(2 * pi) * sigma); % Height of Gaussian
        w = gausswin(L, alpha) * h;
        w(w < h * 1.e-3) = []; % Deleting Gaussian edges
        h_vel = conv(h_vel, w, 'same');
        h_vel(inv) = nan;

        inv = isnan(v_vel);
        v_vel = fillmissing(v_vel, 'linear', 'EndValues', 'none');
        sigma = floor(0.25 * sr); % 50 ms window
        L = sigma * 12; % Length of 6 sigma each side
        alpha = ((L - 1) / sigma) / 2;
        h = 1 / (sqrt(2 * pi) * sigma); % Height of Gaussian
        w = gausswin(L, alpha) * h;
        w(w < h * 1.e-3) = []; % Deleting Gaussian edges
        v_vel = conv(v_vel, w, 'same');
        v_vel(inv) = nan;

        inv = isnan(r_vel);
        r_vel = fillmissing(r_vel, 'linear', 'EndValues', 'none');
        sigma = floor(0.25 * sr); % 50 ms window
        L = sigma * 12; % Length of 6 sigma each side
        alpha = ((L - 1) / sigma) / 2;
        h = 1 / (sqrt(2 * pi) * sigma); % Height of Gaussian
        w = gausswin(L, alpha) * h;
        w(w < h * 1.e-3) = []; % Deleting Gaussian edges
        r_vel = conv(r_vel, w, 'same');
        r_vel(inv) = nan;
    end
end

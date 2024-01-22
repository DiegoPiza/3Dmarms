function [ang_vel] = HeadSpeed(TrackingData, f, gaus_kernel)
    % Calculate angular velocity of head movements
    
    % Input:
    % TrackingData - Struct containing tracking data including quaternion values
    % f - Sampling rate in seconds
    % gaus_kernel - true if smoothing with Gaussian kernel
    
    % Output:
    % ang_vel - Angular velocity in degrees/second
    
    q = ([TrackingData.WQuat TrackingData.iXQuat TrackingData.jYQuat TrackingData.kZQuat]);
    q1 = quaternion(q);
    ang_vel = nan([length(TrackingData.WQuat) 1]);
    ang_vel(2:end) = rad2deg(abs(dist(q1(1:end-1), q1(2:end))))./f;
    ang_vel(ang_vel > 2000, :) = nan;

    if gaus_kernel
        inv = isnan(ang_vel);
        ang_vel = fillmissing(ang_vel, 'linear', 'EndValues', 'none');
        sr = 1 / f; % Convert sr to hertz
        
        % Apply position-based filters
        ang_vel(TrackingData.XPosition > 0.6, :) = nan;
        ang_vel(TrackingData.XPosition < 0, :) = nan;
        ang_vel(TrackingData.YPosition < 0, :) = nan; 
        ang_vel(TrackingData.YPosition > 1.14, :) = nan;
        ang_vel(TrackingData.ZPosition < 0.2, :) = nan;
        
        % Gaussian kernel smoothing
        sigma = floor(0.25 * sr); % 250 ms window
        L = sigma * 12; % Length of 6 sigma each side
        alpha = ((L - 1) / sigma) / 2;
        h = 1 / (sqrt(2 * pi) * sigma);   % Height of Gaussian
        w = gausswin(L, alpha) * h; % The exact correspondence with the standard deviation of a Gaussian probability density function is σ = (L – 1)/(2α).
        w(w < h * 1.e-3) = []; % Delete Gaussian edges
        ang_vel = conv(ang_vel, w, 'same');
        ang_vel(inv) = nan;
    end
end

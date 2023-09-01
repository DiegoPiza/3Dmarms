function [eucl, dx, dy, dz] = TracVelocity(TrackingData, f, gaus_kernel)
    % Computes velocity of translation component of tracking data in cm/s
    
    % Input:
    % TrackingData - Data containing tracking information
    % f - Sampling rate (e.g., 1 ms = 0.001)
    % gaus_kernel - Flag indicating whether to smooth with Gaussian kernel
    
    % Output:
    % eucl - Euclidean velocity (cm/s)
    % dx - Velocity along X-axis (cm/s)
    % dy - Velocity along Y-axis (cm/s)
    % dz - Velocity along Z-axis (cm/s)
    
    dx = nan([length(TrackingData.XPosition) 1]);
    dy = nan([length(TrackingData.XPosition) 1]);
    dz = nan([length(TrackingData.XPosition) 1]);
    
    dx(2:end) = diff(TrackingData.XPosition * 100) / f; % Conversion to cm/s
    dy(2:end) = diff(TrackingData.YPosition * 100) / f;
    dz(2:end) = diff(TrackingData.ZPosition * 100) / f;
    
    eucl = sqrt((dx.^2) + (dy.^2) + (dz.^2)); % Pythagoras
    
    % Remove outliers and very high values
    eucl(TrackingData.XPosition > 0.6, :) = nan;
    eucl(TrackingData.XPosition < 0, :) = nan;
    eucl(TrackingData.YPosition < 0, :) = nan;
    eucl(TrackingData.YPosition > 1.14, :) = nan;
    eucl(TrackingData.ZPosition < 0.2, :) = nan;
    eucl(eucl > 300) = nan;
    
    % Apply Gaussian smoothing if specified
    if gaus_kernel == 1
        inv = isnan(eucl);
        eucl = fillmissing(eucl, 'linear', 'EndValues', 'none');
        sr = 1 / f; % Converting sampling rate to hertz
        sigma = floor(0.2 * sr); % 200 ms window
        L = sigma * 12; % Length of 6 sigma each side
        alpha = ((L - 1) / sigma) / 2;
        h = 1 / (sqrt(2 * pi) * sigma); % Height of Gaussian
        w = gausswin(L, alpha) * h; % The exact correspondence with the standard deviation of a Gaussian probability density function is σ = (L – 1) / (2α).
        w(w < h * 1.e-3) = []; % Deleting Gaussian edges
        eucl = conv(eucl, w, 'same');
        eucl(inv) = nan;
    elseif gaus_kernel == 2
        inv = isnan(eucl);
        dx = fillmissing(dx, 'linear', 'EndValues', 'none');
        dy = fillmissing(dy, 'linear', 'EndValues', 'none');
        dz = fillmissing(dz, 'linear', 'EndValues', 'none');
        sr = 1 / f; % Converting sampling rate to hertz
        sigma = floor(0.05 * sr); % 250 ms window
        L = sigma * 12; % Length of 6 sigma each side
        alpha = ((L - 1) / sigma) / 2;
        h = 1 / (sqrt(2 * pi) * sigma); % Height of Gaussian
        w = gausswin(L, alpha) * h; % The exact correspondence with the standard deviation of a Gaussian probability density function is σ = (L – 1) / (2α).
        w(w < h * 1.e-3) = []; % Deleting Gaussian edges
        dx = conv(dx, w, 'same');
        dy = conv(dy, w, 'same');
        dz = conv(dz, w, 'same');
        eucl = sqrt((dx.^2) + (dy.^2) + (dz.^2)); % Pythagoras
        eucl(inv) = nan;
    end
    
    % Additional processing steps could be added here
end

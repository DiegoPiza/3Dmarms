%% Speed Time series analysis  https://www.sciencedirect.com/science/article/pii/S2211124718316437#fig3
%Similar to Zé Henrique T.D. Góis - Characterizing Speed Cells in the Rat Hippocampus
rng shuffle                 % Seed the random number generator
clear all                   % Clear workspace
close all                   % Close all figure windows
%%
load sample_data.mat        % Load sample data
%%
speed_score = table;        % Initialize table to store speed scores

numofshuffle = 1000;        % Number of shuffles for permutation test
sr = 60;                    % Sampling rate of time series

% Gaussian kernel parameters
sigma = floor(0.25 * sr);   % Standard deviation (250 ms window)
L = sigma * 12;             % Length of the window (6 sigma on each side)
alpha = ((L - 1) / sigma) / 2;  % Alpha parameter for Gaussian window
h = 1 / (sqrt(2 * pi) * sigma); % Height of the Gaussian
w = gausswin(L, alpha) * h; % Generate Gaussian window
% The exact correspondence with the standard deviation of a Gaussian probability density function is σ = (L – 1)/(2α).
w(w < h * 1.e-3) = [];      % Delete negligible Gaussian edges
plot_figures = false;       % Flag to control plotting

fields = fieldnames(units.singleunits); % Get field names of single units
score = table;                         % Initialize table to store scores
for j = 1:length(fields)
    unitID = string(fields(j));                % Convert field name to string
    neuron = char(fields(j));                  % Convert field name to character array

    unit = units.singleunits.(neuron);         % Access unit data
    if length(unit) < 250
        continue                               % Skip units with less than 250 data points
    end
    score.unitID{j} = unitID;                  % Store unit ID in score table

    %%
    puls = TrackingData.neural_time;           % Get neural time pulses
    close all

    gapt = find(gap) ./ 30000;                 % Find gaps in data (in seconds)

    inv_frm = ismember(puls, gapt);            % Identify invalid frames

    spikes = histc(unit, puls);                % Histogram of spikes aligned to neural time pulses
    heatm = HeadSpeed(TrackingData, 1/60, 1);  % Calculate head speed
    heatm(TrackingData.XPosition > 0.6, :) = nan;    % Set out-of-bounds X positions to NaN
    heatm(TrackingData.XPosition < 0, :) = nan;
    heatm(TrackingData.YPosition < 0, :) = nan;      % Set out-of-bounds Y positions to NaN
    heatm(TrackingData.YPosition > 1.14, :) = nan;
    heatm(TrackingData.ZPosition < 0, :) = nan;      % Set out-of-bounds Z positions to NaN

    heatm(heatm > 2000, :) = nan;              % Remove unrealistic head speed values

    if length(spikes) > length(heatm)          % Check if neural data is longer than tracking data
        spikes = spikes(1:length(heatm));      % Truncate spikes to match heatm length
        inv_frm = inv_frm(1:length(heatm));    % Truncate invalid frames indicator
    else
        heatm = heatm(1:length(spikes), :);    % Truncate heatm to match spikes length
        inv_frm = inv_frm(1:length(spikes), :);
    end

    heatm(inv_frm, :) = nan;                   % Set invalid frames to NaN

    % Gaussian kernel smoothing
    spikes = conv(spikes, w, 'same');          % Convolve spikes with Gaussian window
    % Remove NaNs
    spikes(isnan(heatm(:, 1))) = [];           % Remove spikes where heatm is NaN
    heatm(isnan(heatm(:, 1)), :) = [];         % Remove corresponding heatm entries

    R = corr(spikes, heatm);                   % Calculate correlation between spikes and head speed
    score.Head_real{j} = R;                    % Store real correlation score

    Rsh = zeros([numofshuffle, 1]);            % Initialize array for shuffled correlations
    shift = randi([-90 * sr, 90 * sr], [numofshuffle, 1]); % Random shifts between -90s and +90s

    parfor i = 1:numofshuffle
        ash = circshift(spikes, shift(i));     % Circularly shift spikes
        Rsh(i) = corr(ash, heatm);             % Calculate correlation with shifted spikes
    end
    score.Head_sh{j} = Rsh;                    % Store shuffled correlation scores
    % Statistical analysis
    mu = mean(Rsh);                            % Mean of shuffled correlations
    sigma_stats = std(Rsh);                    % Standard deviation of shuffled correlations
    p_val = 1 - normcdf(R, mu, sigma_stats);   % Calculate p-value
    reject_h0 = p_val < 0.05 / length(spikes); % Hypothesis testing with Bonferroni correction
    score.head_p_val{j} = p_val;               % Store p-value
    score.head_reject_h0{j} = reject_h0;       % Store hypothesis test result
    %% Plotting

    if plot_figures == true
        h_fig2 = figure('units', 'normalized', 'position', [0 0 1 1]); % Create full-screen figure
        subplot(2, 1, 1)
        histogram(Rsh, 'normalization', 'pdf');                        % Plot histogram of shuffled correlations
        hold all
        shuff_xs = linspace(min(Rsh), max(Rsh), 100);                  % X-axis values for normal curve
        shuff_curve = normpdf(shuff_xs, mu, sigma_stats);              % Normal distribution curve
        h_line1 = plot(shuff_xs, shuff_curve, '--b', 'linewidth', 2);  % Plot normal distribution
        h_line2 = line([R R], get(gca, 'ylim'), 'color', 'r', 'linewidth', 3); % Plot line for real correlation
        legend([h_line1, h_line2], {'Shuffled'; 'Real'}, 'location', 'northeastoutside')
        title('Head Speed Distribution of shuffled total Speed Score')
        xlabel({'Speed Score (Corr)'})
        ylabel('Probability density')
        text(0.2, 10, char(['Head Speed Score: ' num2str(R)]))         % Display real correlation
        text(0.2, 15, char(['Permutation test p-value: ' num2str(p_val)])) % Display p-value
    end

    spikes = histc(unit, puls);                % Recalculate spikes histogram

    heatm = TracVelocity(TrackingData, 1/60, 1); % Calculate body velocity
    heatm(TrackingData.XPosition > 0.6, :) = nan; % Set out-of-bounds X positions to NaN
    heatm(TrackingData.XPosition < 0, :) = nan;
    heatm(TrackingData.YPosition < 0, :) = nan;   % Set out-of-bounds Y positions to NaN
    heatm(TrackingData.YPosition > 1.14, :) = nan;
    heatm(TrackingData.ZPosition < 0, :) = nan;
    heatm(heatm > 300, :) = nan;               % Remove unrealistic body speed values

    if length(spikes) > length(heatm)          % Check if neural data is longer than tracking data
        spikes = spikes(1:length(heatm));      % Truncate spikes to match heatm length
        inv_frm = inv_frm(1:length(heatm));    % Truncate invalid frames indicator
    else
        heatm = heatm(1:length(spikes), :);    % Truncate heatm to match spikes length
        inv_frm = inv_frm(1:length(spikes), :);
    end

    heatm(inv_frm, :) = nan;                   % Set invalid frames to NaN

    % Gaussian kernel smoothing
    spikes = conv(spikes, w, 'same');          % Convolve spikes with Gaussian window
    % Remove NaNs
    spikes(isnan(heatm(:, 1))) = [];           % Remove spikes where heatm is NaN
    heatm(isnan(heatm(:, 1)), :) = [];         % Remove corresponding heatm entries

    R = corr(spikes, heatm);                   % Calculate correlation between spikes and body speed
    score.Body_real{j} = R;                    % Store real correlation score

    Rsh = zeros([numofshuffle, 1]);            % Initialize array for shuffled correlations
    shift = randi([-90 * sr, 90 * sr], [numofshuffle, 1]); % Random shifts between -90s and +90s

    parfor i = 1:numofshuffle
        ash = circshift(spikes, shift(i));     % Circularly shift spikes
        Rsh(i) = corr(ash, heatm);             % Calculate correlation with shifted spikes
    end

    score.Body_sh{j} = Rsh;                    % Store shuffled correlation scores
    % Statistical analysis
    mu = mean(Rsh);                            % Mean of shuffled correlations
    sigma_stats = std(Rsh);                    % Standard deviation of shuffled correlations
    p_val = 1 - normcdf(R, mu, sigma_stats);   % Calculate p-value
    reject_h0 = p_val < 0.05 / length(spikes); % Hypothesis testing with Bonferroni correction
    score.body_p_val{j} = p_val;               % Store p-value
    score.body_reject_h0{j} = reject_h0;       % Store hypothesis test result

    %% Plotting

    if plot_figures == true
        subplot(2, 1, 2)

        histogram(Rsh, 'normalization', 'pdf');                        % Plot histogram of shuffled correlations
        hold all
        shuff_xs = linspace(min(Rsh), max(Rsh), 100);                  % X-axis values for normal curve
        shuff_curve = normpdf(shuff_xs, mu, sigma_stats);              % Normal distribution curve
        h_line1 = plot(shuff_xs, shuff_curve, '--b', 'linewidth', 2);  % Plot normal distribution
        h_line2 = line([R R], get(gca, 'ylim'), 'color', 'r', 'linewidth', 3); % Plot line for real correlation
        legend([h_line1, h_line2], {'Shuffled'; 'Real'}, 'location', 'northeastoutside')
        title('Body Speed Distribution of shuffled total Speed Score')
        xlabel({'Speed Score (Corr)'})
        ylabel('Probability density')
        text(0.2, 10, char(['Body Vel Speed Score: ' num2str(R)]))     % Display real correlation
        text(0.2, 15, char(['Permutation test p-value: ' num2str(p_val)])) % Display p-value
    end

end

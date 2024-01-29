%% Figure S8a
 load('sampleLFP.mat')
% Define parameters
fs = 1000;
filter_order = 4;
ct = [4 10];
ct = ct / (fs/2);

% Design a bandpass Butterworth filter
[b, a] = butter(filter_order, ct, 'bandpass');

% Apply the filter to the data
bpasstheta = filtfilt(b, a, double(head_movLFP)')';

% Compute the analytic signal
y = hilbert(bpasstheta')';

% Extract phase information
sigphase = abs(unwrap(angle(y), [], 2));

% Select a subset of data
a = sigphase(:, 300:700);

% Initialize variables
pval = nan([1 size(a, 2)]);
z = nan([1 size(a, 2)]);
h = [];

% Loop over columns
for i = 1:size(a, 2)
    % Circular statistics test
    [pval(i), z(i)] = circ_rtest(a(:, i));

    % Histogram computation
    [h(:, i), binedges] = histcounts(wrapTo360(rad2deg(a(:, i))), [0:5:360], 'Normalization', 'probability');
end

% Smoothed histogram plot
hsmooth = imgaussfilt(h, 0.85);
imagesc(hsmooth)
set(gca, 'TickDir', 'out');
set(gca, 'YDir', 'normal')
xticks([0:50:400])
xticklabels([-200:50:200])
yticks([0:18:length(h(:, 1))])
yticklabels([0:90:360])
colorbar
set(gca, 'TickDir', 'out');
line([80 80], [0 360], 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--')
line([320 320], [0 360], 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--')


%% Figure S8b

ang=wrapToPi(a(:,320));
bin_n=10;
r=circ_r(ang);
mu=circ_mean(ang);
[h edge]=histcounts(ang,bin_n);
edge_center=edge(1:8);
mr = max(h);
r = r * mr; %from Circular stats toolbox
phi = mu;
zm = r*exp(1i*phi);
% Coordinates of the two points
x1 = 0;
y1 = 0;
x2 = real(zm);
y2 = imag(zm);
% Calculate the differences in x and y coordinates
deltaX = x2 - x1;
deltaY = y2 - y1;
% Calculate the length of the vector using Euclidean distance formula
lengthr = sqrt(deltaX.^2 + deltaY.^2);
% Calculate the angle of the vector in radians
angleRad = atan2(deltaY, deltaX);
% Convert the angle to degrees
angleDeg = rad2deg(angleRad);
his=polarhistogram(ang,bin_n,'Normalization','probability');

hold on

ang=wrapToPi(a(:,80));
bin_n=10;
r=circ_r(ang);
mu=circ_mean(ang);
[h edge]=histcounts(ang,bin_n);
edge_center=edge(1:8);
mr = max(h);
r = r * mr; %from Circular stats toolbox
phi = mu;
zm = r*exp(1i*phi);
% Coordinates of the two points
x1 = 0;
y1 = 0;
x2 = real(zm);
y2 = imag(zm);
% Calculate the differences in x and y coordinates
deltaX = x2 - x1;
deltaY = y2 - y1;
% Calculate the length of the vector using Euclidean distance formula
lengthr = sqrt(deltaX.^2 + deltaY.^2);
% Calculate the angle of the vector in radians
angleRad = atan2(deltaY, deltaX);
% Convert the angle to degrees
angleDeg = rad2deg(angleRad);
his=polarhistogram(ang,bin_n,'Normalization','probability');



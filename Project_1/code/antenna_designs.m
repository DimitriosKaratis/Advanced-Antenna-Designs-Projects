%% Advanced Antenna Designs
% Assignment 1 - April 2025
% KARATIS DIMITRIOS 10775

clc;
clear;
close all;

%% ========================= COMMON SETTINGS =========================
freq = linspace(1.5e9, 6.5e9, 100);  % Frequency range
Z0 = 50;  % Characteristic impedance

% ========== H Antenna Parameters ==========
Wf = 1.5e-3;
Lf = 7e-3;
Wh = 7e-3;
Lh = 10e-3;
Wh1 = 6e-3;
Lh1 = 9e-3;
Wh2 = 2e-3;
Lh2 = 4e-3;
Wh3 = 1e-3;
Lh3 = 4e-3;
Wh4 = 2.25e-3;

h_sub = 0.8e-3;
Wsub = 12e-3;
Lsub = 18e-3;
Lgnd = 4e-3;

% ========== Dielectric Parameters ==========
epsilor_r = 4.4;
loss_tangent = 0.018;
thickness = 0.8e-3;

%% ========================= ANTENNA 1 =========================
% Antenna 1 geometry (Rectangular patch)

L_patch = Lh;
W_patch = Wh1 + 2 * Wh2;

L_feed_line = Lf - Lh3 + Lgnd;
W_feed_line = Wf;

x_middle_patch = Wsub/2;
y_middle_patch = Lgnd + Lf - Lh3 + Lh/2;


x_middle_feed_line = x_middle_patch;
y_middle_feed_line = (Lgnd + Lf - Lh3)/2;

Wgnd = Wsub;

x_middle_ground = x_middle_patch;
y_middle_ground = Lgnd/2;

x_middle_feed = x_middle_patch;
y_middle_feed = Lgnd;

x_middle_pcb = x_middle_patch;
y_middle_pcb = Lsub/2;

% ========== Substrate ==========
substrate = dielectric('Name', 'FR4', ...
                       'EpsilonR', epsilor_r, ...
                       'LossTangent', loss_tangent, ...
                       'Thickness', thickness);  

% ========== Building the squared-shaped PEC patch ==========
square_patch = antenna.Rectangle('Length', W_patch, 'Width', L_patch, ...
                                'Center', [x_middle_patch, y_middle_patch]);

% ========== Feed Line ==========
feedLine = antenna.Rectangle('Length', W_feed_line, 'Width', L_feed_line, ...
                              'Center', [x_middle_feed_line, y_middle_feed_line]);

% Combine the feed line, PEC line, and square-shaped patch
patchWithFeedAndPEC = square_patch + feedLine;

% ========== Ground Plane ==========
groundPlane = antenna.Rectangle('Length', Wgnd, 'Width', Lgnd, ...
                                'Center', [x_middle_ground, y_middle_ground]);

% ========== PCB Shape ==========
pcbShape = antenna.Rectangle('Length', Wsub, 'Width', Lsub, ...
                                'Center', [x_middle_pcb, y_middle_pcb]);

% ========== PCB Stack ==========
ant1 = pcbStack;
ant1.BoardThickness = h_sub;
ant1.BoardShape = pcbShape;
ant1.Layers = {patchWithFeedAndPEC, substrate, groundPlane};
ant1.FeedLocations = [x_middle_feed, y_middle_feed, 1, 3]; 
ant1.FeedDiameter = 0.08e-3; 

% ========== Visualization ==========
figure;
show(ant1);
title('Square Monopole Antenna');

% S-parameters for Antenna 1
mesh(ant1, 'MaxEdgeLength', 0.01, 'MinEdgeLength', 0.001, 'GrowthRate', 0.7);

s1 = sparameters(ant1, freq, Z0);
s11_1 = 20*log10(abs(s1.Parameters(1,1,:)));
s11_1 = squeeze(s11_1);

%% ========================= ANTENNA 2 =========================
% Antenna 2 geometry (H-shaped patch)

L_feed_line = Lf + Lgnd;
W_feed_line = Wf;

x_middle_patch = Wsub/2;
y_middle_patch = L_feed_line - Lh3 + Lh/2;


x_middle_feed_line = x_middle_patch;
y_middle_feed_line = L_feed_line/2;

Wgnd = Wsub;

x_middle_ground = x_middle_patch;
y_middle_ground = Lgnd/2;

x_middle_feed = x_middle_patch;
y_middle_feed = Lgnd;

x_middle_pcb = x_middle_patch;
y_middle_pcb = Lsub/2;

y_middle_left = y_middle_patch;
x_middle_left = x_middle_patch - Wh/2 - Wh3/2;
y_middle_right = y_middle_patch;
x_middle_right = x_middle_patch + Wh/2 + Wh3/2;
x_middle_center = x_middle_patch;
y_middle_center = y_middle_patch;

L_left = Lh;
W_left = Wh2;
L_center = Lh - 2*Lh3;
W_center = Wh1;
L_right = L_left;
W_right = W_left;

% ========== Substrate ==========
substrate = dielectric('Name', 'FR4', ...
                       'EpsilonR', epsilor_r, ...
                       'LossTangent', loss_tangent, ...
                       'Thickness', thickness);  

% ========== Building the H-shaped PEC patch ==========
% Vertical connector of the H
left_part = antenna.Rectangle('Length', W_left, 'Width', L_left, ...
                                'Center', [x_middle_left, y_middle_left]);
% Horizontal connector of the H
center_part = antenna.Rectangle('Length', W_center, 'Width', L_center, ...
                                'Center', [x_middle_center, y_middle_center]);
% Vertical connector of the H
right_part = antenna.Rectangle('Length', W_right, 'Width', L_right, ...
                                'Center', [x_middle_right, y_middle_right]);

% ========== Feed Line ==========
feedLine = antenna.Rectangle('Length', W_feed_line, 'Width', L_feed_line, ...
                              'Center', [x_middle_feed_line, y_middle_feed_line]);

% Combine the feed line, PEC line, and H-shaped patch
patchWithFeedAndPEC = left_part + center_part + right_part + feedLine;

% ========== Ground Plane ==========
groundPlane = antenna.Rectangle('Length', Wgnd, 'Width', Lgnd, ...
                                'Center', [x_middle_ground, y_middle_ground]);

% ========== PCB Shape ==========
pcbShape = antenna.Rectangle('Length', Wsub, 'Width', Lsub, ...
                                'Center', [x_middle_pcb, y_middle_pcb]);

% ========== PCB Stack ==========
ant2 = pcbStack;
ant2.BoardThickness = h_sub;
ant2.BoardShape = pcbShape;
ant2.Layers = {patchWithFeedAndPEC, substrate, groundPlane};
ant2.FeedLocations = [x_middle_feed, y_middle_feed, 1, 3]; 
ant2.FeedDiameter = 0.08e-3;

% ========== Visualization ==========
figure;
show(ant2);
title('H-Ring Monopole Antenna');

% S-parameters for Antenna 2
mesh(ant2, 'MaxEdgeLength', 0.01, 'MinEdgeLength', 0.001, 'GrowthRate', 0.7);

s2 = sparameters(ant2, freq, Z0);
s11_2 = 20*log10(abs(s2.Parameters(1,1,:)));
s11_2 = squeeze(s11_2);

%% ========================= ANTENNA 3 =========================
% Antenna 3 geometry (H proposed antenna)

x = (Wh2 - Wh3)/2;

L_feed_line = Lf + Lgnd;
W_feed_line = Wf;

x_middle_patch = Wsub/2;
y_middle_patch = L_feed_line - Lh3 + Lh/2;

x_middle_feed_line = x_middle_patch;
y_middle_feed_line = L_feed_line/2;

Wgnd = Wsub;

x_middle_ground = x_middle_patch;
y_middle_ground = Lgnd/2;

x_middle_feed = x_middle_patch;
y_middle_feed = Lgnd;

x_middle_pcb = x_middle_patch;
y_middle_pcb = Lsub/2;


x_middle_left_vertical_big = x_middle_patch - Wh/2 - Wh3 - x/2;
y_middle_left_vertical_big = y_middle_patch;
L_left_vertical_big = Lh;
W_left_vertical_big = x;

x_middle_right_vertical_big = x_middle_patch + Wh/2 + Wh3 + x/2;
y_middle_right_vertical_big = y_middle_patch;
L_right_vertical_big = Lh;
W_right_vertical_big = x;


x_middle_left_vertical_small_top = x_middle_patch - Wh1/2 - x/2;
y_middle_left_vertical_small_top = y_middle_patch + Lh/2 - x - Lh2/2;
L_left_vertical_small_top = Lh2;
W_left_vertical_small_top = x;

x_middle_right_vertical_small_top = x_middle_patch + Wh1/2 + x/2;
y_middle_right_vertical_small_top = y_middle_patch + Lh/2 - x - Lh2/2;
L_right_vertical_small_top = Lh2;
W_right_vertical_small_top = x;


x_middle_left_vertical_small_bottom = x_middle_patch - Wh1/2 - x/2;
y_middle_left_vertical_small_bottom = y_middle_patch - Lh/2 + x + Lh2/2;
L_left_vertical_small_bottom = Lh2;
W_left_vertical_small_bottom = x;

x_middle_right_vertical_small_bottom = x_middle_patch + Wh1/2 + x/2;
y_middle_right_vertical_small_bottom = y_middle_patch - Lh/2 + x + Lh2/2;
L_right_vertical_small_bottom = Lh2;
W_right_vertical_small_bottom = x;


x_middle_left_horizontal_small_top = x_middle_patch - Wh/2 - Wh3/2;
y_middle_left_horizontal_small_top = y_middle_patch + Lh/2 - x/2;
L_left_horizontal_small_top = x;
W_left_horizontal_small_top = Wh2; 

x_middle_right_horizontal_small_top = x_middle_patch + Wh/2 + Wh3/2;
y_middle_right_horizontal_small_top = y_middle_patch + Lh/2 - x/2;
L_right_horizontal_small_top = x;
W_right_horizontal_small_top = Wh2; 


x_middle_left_horizontal_small_bottom = x_middle_patch - Wh/2 - Wh3/2;
y_middle_left_horizontal_small_bottom = y_middle_patch - Lh/2 + x/2;
L_left_horizontal_small_bottom = x;
W_left_horizontal_small_bottom = Wh2; 

x_middle_right_horizontal_small_bottom = x_middle_patch + Wh/2 + Wh3/2;
y_middle_right_horizontal_small_bottom = y_middle_patch - Lh/2 + x/2;
L_right_horizontal_small_bottom = x;
W_right_horizontal_small_bottom = Wh2; 


x_middle_horizontal_big_top = x_middle_patch;
y_middle_horizontal_big_top = y_middle_patch + Lh/2 - Lh3 - x/2;
L_horizontal_big_top = x;
W_horizontal_big_top = Wh1; 


x_middle_horizontal_big_bottom = x_middle_patch;
y_middle_horizontal_big_bottom = y_middle_patch - Lh/2 + Lh3 + x/2;
L_horizontal_big_bottom = x;
W_horizontal_big_bottom = Wh1; 


% ========== Substrate ==========
substrate = dielectric('Name', 'FR4', ...
                       'EpsilonR', epsilor_r, ...
                       'LossTangent', loss_tangent, ...
                       'Thickness', thickness);  

% ========== Building the H-shaped PEC patch ==========
left_vertical_big = antenna.Rectangle('Length', W_left_vertical_big, 'Width', L_left_vertical_big, ...
                                'Center', [x_middle_left_vertical_big, y_middle_left_vertical_big]);

right_vertical_big = antenna.Rectangle('Length', W_right_vertical_big, 'Width', L_right_vertical_big, ...
                                'Center', [x_middle_right_vertical_big, y_middle_right_vertical_big]);

left_vertical_small_top = antenna.Rectangle('Length', W_left_vertical_small_top, 'Width', L_left_vertical_small_top, ...
                                'Center', [x_middle_left_vertical_small_top, y_middle_left_vertical_small_top]);


right_vertical_small_top = antenna.Rectangle('Length', W_right_vertical_small_top, 'Width', L_right_vertical_small_top, ...
                                'Center', [x_middle_right_vertical_small_top, y_middle_right_vertical_small_top]);


left_vertical_small_bottom = antenna.Rectangle('Length', W_left_vertical_small_bottom, 'Width', L_left_vertical_small_bottom, ...
                                'Center', [x_middle_left_vertical_small_bottom, y_middle_left_vertical_small_bottom]);

right_vertical_small_bottom = antenna.Rectangle('Length', W_right_vertical_small_bottom, 'Width', L_right_vertical_small_bottom, ...
                                'Center', [x_middle_right_vertical_small_bottom, y_middle_right_vertical_small_bottom]);

left_horizontal_small_top = antenna.Rectangle('Length', W_left_horizontal_small_top, 'Width', L_left_horizontal_small_top, ...
                                'Center', [x_middle_left_horizontal_small_top, y_middle_left_horizontal_small_top]);

right_horizontal_small_top = antenna.Rectangle('Length', W_right_horizontal_small_top, 'Width', L_right_horizontal_small_top, ...
                                'Center', [x_middle_right_horizontal_small_top, y_middle_right_horizontal_small_top]);

left_horizontal_small_bottom = antenna.Rectangle('Length', W_left_horizontal_small_bottom, 'Width', L_left_horizontal_small_bottom, ...
                                'Center', [x_middle_left_horizontal_small_bottom, y_middle_left_horizontal_small_bottom]);

right_horizontal_small_bottom = antenna.Rectangle('Length', W_right_horizontal_small_bottom, 'Width', L_right_horizontal_small_bottom, ...
                                'Center', [x_middle_right_horizontal_small_bottom, y_middle_right_horizontal_small_bottom]);

horizontal_big_top = antenna.Rectangle('Length', W_horizontal_big_top, 'Width', L_horizontal_big_top, ...
                                'Center', [x_middle_horizontal_big_top, y_middle_horizontal_big_top]);

horizontal_big_bottom = antenna.Rectangle('Length', W_horizontal_big_bottom, 'Width', L_horizontal_big_bottom, ...
                                'Center', [x_middle_horizontal_big_bottom, y_middle_horizontal_big_bottom]);

% ========== Feed Line ==========
feedLine = antenna.Rectangle('Length', W_feed_line, 'Width', L_feed_line, ...
                              'Center', [x_middle_feed_line, y_middle_feed_line]);

% Combine the feed line, PEC line, and H-shaped proposed patch
patchWithFeedAndPEC = left_vertical_big + right_vertical_big + left_vertical_small_top ...
+ right_vertical_small_top + left_vertical_small_bottom + right_vertical_small_bottom ...
+ left_horizontal_small_top + right_horizontal_small_top + left_horizontal_small_bottom ...
+ right_horizontal_small_bottom + horizontal_big_top + horizontal_big_bottom + feedLine;

% ========== Ground Plane ==========
groundPlane = antenna.Rectangle('Length', Wgnd, 'Width', Lgnd, ...
                                'Center', [x_middle_ground, y_middle_ground]);

% ========== PCB Shape ==========
pcbShape = antenna.Rectangle('Length', Wsub, 'Width', Lsub, ...
                                'Center', [x_middle_pcb, y_middle_pcb]);

% ========== PCB Stack ==========
ant3 = pcbStack;
ant3.BoardThickness = h_sub;
ant3.BoardShape = pcbShape;
ant3.Layers = {patchWithFeedAndPEC, substrate, groundPlane};
ant3.FeedLocations = [x_middle_feed, y_middle_feed, 1, 3]; 
ant3.FeedDiameter = 0.08e-3;

% ========== Visualization ==========
figure;
show(ant3);
title('Proposed H-Ring Monopole Antenna');


% S-parameters for Antenna 3
mesh(ant3, 'MaxEdgeLength', 0.01, 'MinEdgeLength', 0.001, 'GrowthRate', 0.7);

s3 = sparameters(ant3, freq, Z0);
s11_3 = 20*log10(abs(s3.Parameters(1,1,:)));
s11_3 = squeeze(s11_3);

%% Radiation Pattern at specific frequencies
freqPlot = 3.4e9; 
pattern_data_3_1 = pattern(ant3, freqPlot);
[co_pol_yz_3_1, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 90, 'CoordinateSystem', 'polar', 'Polarization', 'V');
[co_pol_xz_3_1, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 0, 'CoordinateSystem', 'polar', 'Polarization', 'V');
[cross_pol_yz_3_1, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 90, 'CoordinateSystem', 'polar', 'Polarization', 'H');
[cross_pol_xz_3_1, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 0, 'CoordinateSystem', 'polar', 'Polarization', 'H');

freqPlot = 5.6e9; 
pattern_data_3_2 = pattern(ant3, freqPlot);
[co_pol_yz_3_2, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 90, 'CoordinateSystem', 'polar', 'Polarization', 'V');
[co_pol_xz_3_2, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 0, 'CoordinateSystem', 'polar', 'Polarization', 'V');
[cross_pol_yz_3_2, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 90, 'CoordinateSystem', 'polar', 'Polarization', 'H');
[cross_pol_xz_3_2, ~] = pattern(ant3, freqPlot, linspace(0, 360, 361), 0, 'CoordinateSystem', 'polar', 'Polarization', 'H');


%% ========================= Combined S11 Plot =========================
% S11 plot for all antennas in one figure for comparison
figure;
hold on;
plot(freq, s11_1, 'b-', 'LineWidth', 1.5); 
plot(freq, s11_2, 'g-', 'LineWidth', 1.5);
plot(freq, s11_3, 'r-', 'LineWidth', 1.5);

% Add a horizontal dashed black line at -10 dB
yline(-10, 'k--', 'LineWidth', 1.2); 

title('S11 Comparison for Antennas 1, 2, and 3');
xlabel('Frequency (Hz)');
ylabel('S11 (dB)');
legend('Basic Structure (Ordinary Square Monopole Antenna)', 'H-Shaped Monopole Antenna', 'The Proposed Monopole Antenna');
grid on;
hold off;

%% Smith Chart for the proposed antenna
smithFig = figure;
h = smithplot(s3);
title('Simulated Input Impedance - Smith Chart for the proposed antenna');
grid on;
hold on;

% --- Add SWR = 2 circle ---
theta = linspace(0, 2*pi, 500);
radius = 0.333;
x = radius * cos(theta);
y = radius * sin(theta);
swrCircleHandle = plot(x, y, 'r--', 'LineWidth', 1.5); 

% Threshold for SWR = 2
swr_threshold = 0.333;

% Frequencies to test (Hz)
highlight_freqs = [3.4e9, 5.6e9];
colors = ['r', 'g', 'b', 'm'];

color_idx = 1;
pointHandles = []; 
legendLabels = {'SWR Circle'}; 

for k = 1:length(highlight_freqs)
    [~, idx] = min(abs(freq - highlight_freqs(k)));
    gamma = s3.Parameters(1,1,idx);
    gamma_mag = abs(gamma);

    if gamma_mag <= swr_threshold
        pt = plot(real(gamma), imag(gamma), 'o', ...
                  'MarkerSize', 8, ...
                  'MarkerFaceColor', colors(color_idx), ...
                  'MarkerEdgeColor', colors(color_idx));

        text(real(gamma) + 0.05, imag(gamma), ...
             sprintf('%.2f GHz', highlight_freqs(k)/1e9), ...
             'Color', colors(color_idx), ...
             'FontSize', 10);

        pointHandles(end+1) = pt;
        legendLabels{end+1} = sprintf('%.2f GHz', highlight_freqs(k)/1e9);

        color_idx = color_idx + 1;
    end
end

% Add legend 
legend([swrCircleHandle, pointHandles], legendLabels);

hold off;

%% ============ Radiation Pattern Plots for Antenna 3 (Proposed) ===============
figure('Position', [100, 100, 800, 600]);

% Plot y-z plane (E-plane) - Subplot (a) (3.4 GHz)
subplot(2,2,1);
polarplot(deg2rad(linspace(0, 360, 361)), co_pol_yz_3_1, 'b-', 'LineWidth', 1.5); 
hold on;
polarplot(deg2rad(linspace(0, 360, 361)), cross_pol_yz_3_1, 'r--', 'LineWidth', 1.5);
title('(a) y-z plane (E-plane) (3.4 GHz)');
rlim([-40 0]); 
thetalim([0 360]);
legend('Co-pol', 'Cross-pol', 'Location', 'southoutside');
grid on;

% Plot x-z plane (H-plane) - Subplot (a) (3.4 GHz)
subplot(2,2,2);
polarplot(deg2rad(linspace(0, 360, 361)), co_pol_xz_3_1, 'b-', 'LineWidth', 1.5); 
hold on;
polarplot(deg2rad(linspace(0, 360, 361)), cross_pol_xz_3_1, 'r--', 'LineWidth', 1.5);
title('(a) x-z plane (H-plane) (3.4 GHz)');
rlim([-40 0]);
thetalim([0 360]);
legend('Co-pol', 'Cross-pol', 'Location', 'southoutside');
grid on;

% Plot y-z plane (E-plane) - Subplot (b) (5.6 GHz)
subplot(2,2,3);
polarplot(deg2rad(linspace(0, 360, 361)), co_pol_yz_3_2, 'b-', 'LineWidth', 1.5); 
hold on;
polarplot(deg2rad(linspace(0, 360, 361)), cross_pol_yz_3_2, 'r--', 'LineWidth', 1.5);
title('(b) y-z plane (E-plane) (5.6 GHz)');
rlim([-40 0]); 
thetalim([0 360]);
legend('Co-pol', 'Cross-pol', 'Location', 'southoutside');
grid on;

% Plot x-z plane (H-plane) - Subplot (b) (5.6 GHz)
subplot(2,2,4);
polarplot(deg2rad(linspace(0, 360, 361)), co_pol_xz_3_2, 'b-', 'LineWidth', 1.5); 
hold on;
polarplot(deg2rad(linspace(0, 360, 361)), cross_pol_xz_3_2, 'r--', 'LineWidth', 1.5);
title('(b) x-z plane (H-plane) (5.6 GHz)');
rlim([-40 0]);
thetalim([0 360]);
legend('Co-pol', 'Cross-pol', 'Location', 'southoutside');
grid on;

% Add a main title for the figure
sgtitle('Figure 7: Radiation Patterns of the Proposed Antenna (3.4 and 5.6 GHz)', 'FontWeight', 'bold');


%% ============== Figure 8: Gain vs Frequency ==============
% Frequency ranges
freq_low = linspace(2.5e9, 4.3e9, 20);  
freq_high = linspace(4.5e9, 6.3e9, 20);    

% Initialize gain storage
gain_low = zeros(size(freq_low));
gain_high = zeros(size(freq_high));

% Calculate gains for lower band
for i = 1:length(freq_low)
    [gain, ~] = pattern(ant3, freq_low(i), 0, 0, 'Type', 'gain');
    gain_low(i) = gain;
end

% Calculate gains for higher band
for i = 1:length(freq_high)
    [gain, ~] = pattern(ant3, freq_high(i), 0, 0, 'Type', 'gain');
    gain_high(i) = gain;
end

% Create figure
figure('Position', [100, 100, 800, 600]);

% Subplot (a): Lower band
subplot(2,1,1);
plot(freq_low/1e9, gain_low, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Freq [GHz]');
ylabel('Gain [dBi]');
title('(a)');
grid on;

% Subplot (b): Higher band
subplot(2,1,2);
plot(freq_high/1e9, gain_high, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Freq [GHz]');
ylabel('Gain [dBi]');
title('(b)');
grid on;

% Main title
sgtitle('Figure 8: Measured Maximum Gain for the Proposed Antenna', 'FontWeight', 'bold');

%% ========================= Surface Current Distribution =========================
% Create analysis at the desired frequencies
freq1 = 3.4e9;
freq2 = 4.3e9;
freq3 = 5.6e9;

figure;
current(ant2,freq2,Slicer="on");

figure;
current(ant3,freq1,Slicer="on");

figure;
current(ant3,freq3,Slicer="on");

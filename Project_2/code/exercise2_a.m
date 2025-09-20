%% Advanced Antenna Designs
% Assignment 2 - Exercise 2a - June 2025
% KARATIS DIMITRIOS 10775

clear all; clc; close all;  

%% Simulation parameters
lambda = 1550e-9;               % Wavelength in meters (1550 nm)
c = 3e8;                        % Speed of light [m/s]
freq = c / lambda;              % Operating frequency [Hz]
omega = 2 * pi * freq;          % Angular frequency
k = 2 * pi / lambda;            % Wavenumber
d = lambda / 2;                 % Distance between elements [m]
positions = (-3.5:1:3.5) * d;   % Equivalent manual definition of positions

T = 1/freq;               % Period of the wave [s]
az = 8;                   % Number of antenna elements in the array

Ns = 30;                  % Number of spatial samples per wavelength
ds = lambda/Ns;           % Spatial discretization step
Nt = 35;                  % Number of temporal samples per period
dt = T/Nt;                % Temporal discretization step
t = 0:dt:(1*T);           % Time vector (1 period)

% Observation space (radial and angular discretization)
R = (0*lambda):ds:(8*lambda);  % Radius from origin [m]
Ntheta = 240;                  % Angular resolution
dtheta = 2*pi/Ntheta;          % Angular step
theta = 0:dtheta:(2*pi);       % Angular vector [rad]

%% Given phase profiles from Table 1 (in radians)
phases_rad = [
    0.0  0.0     0.0     0.0     0.0     0.0     0.0     0.0;     % θ_0 = 0 deg
    0.0  2.657   2.1724  1.6878  1.2032  0.7186  0.234   5.5319;  % θ_0 = 30 deg
    0.0  2.9921  5.9842  2.6931  5.6852  2.3941  5.3862  2.0951;  % θ_0 = 60 deg
    0.0  1.4077  2.8153  4.2230  5.6307  0.7551  2.1628  3.5705   % θ_0 = 90 deg
];

angle_deg = [0, 30, 60, 90];  % Corresponding steering angles [degrees]

%% Visualization
figure('Color','w'); Fs = 12;  

for idx = 1:4
    phi_n = phases_rad(idx, :);  
    AF = zeros(size(theta));     

    for n = 1:az
        % Add contribution of each antenna element to the Array Factor
        AF = AF + exp(1j * (k * positions(n) * sin(theta) - phi_n(n)));
    end

    AF = abs(AF);          
    AF = AF / max(AF);    

    % Polar plot of the normalized AF
    subplot(1,4,idx);
    polarplot(theta, -AF, 'LineWidth', 2);  
    title(['AF for θ₀ = ' num2str(angle_deg(idx)) '°'], 'FontSize', Fs);
end

sgtitle('Exercise 2a – Normalized Array Factor (AF) for 0°, 30°, 60°, 90°', ...
        'FontSize', Fs+2);

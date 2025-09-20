%% Advanced Antenna Designs
% Assignment 2 - Exercise 1b - June 2025
% KARATIS DIMITRIOS 10775

clear all; clc; close all;

%% Simulation parameters
freq = 1e9;               % Operating frequency [Hz]
c = 3e8;                  % Speed of light [m/s]
lambda = c/freq;          % Wavelength [m]
T = 1/freq;               % Period [s]
omega = 2*pi*freq;        % Angular frequency [rad/s]
k = 2*pi/lambda;          % Wavenumber [rad/m]

Ns = 30;                  % Number of spatial samples per wavelength
ds = lambda/Ns;           % Spatial discretization step
Nt = 35;                  % Number of temporal samples per period
dt = T/Nt;                % Temporal discretization step
t = 0:dt:(1*T);           % Time vector (1 period)

% Observation domain
R = (0*lambda):ds:(8*lambda);
Ntheta = 240;
dtheta = 2*pi/Ntheta;
theta = 0:dtheta:(2*pi);
x = R.'*cos(theta);
y = R.'*sin(theta);


% Antenna sizes and angles
az_all = [5, 9];                % Array sizes
angle_deg = [30, 90];           % Beam steering angles
d = lambda/2;                   % Element spacing

% Loop over array sizes (5 and 9 elements)
for az_idx = 1:length(az_all)
    az = az_all(az_idx);        % Current number of elements
    half = floor(az/2);
    r = (-half:1:(az-1-half)) * d;  % Element x-positions

    figure('Color','w');        % One figure per array size

    for ang = 1:length(angle_deg)
        theta0 = angle_deg(ang) * pi/180;
        delta = -k * d * sin(theta0);
        it = 1;

        %% Compute field E(x,y)
        E = zeros(size(x));
        for n = 1:az
            rx = r(n);
            for ix = 1:length(R)
                for iy = 1:length(theta)
                    Rn = sqrt((x(ix,iy)-rx)^2 + y(ix,iy)^2);
                    E(ix,iy) = E(ix,iy) + cos(omega*t(it) - k*Rn - delta*(n - ceil(az/2)));
                end
            end
        end

        %% Compute array factor
        A = ones(1, az);
        Fa = zeros(1, length(theta));
        for i = 0:(length(A)-1)
            Fa = Fa + (A(i+1) * exp(1i*i*delta + 1i*k*(i*d - ((length(A)-1)/2)*d)*cos(theta)));
        end
        Fa = abs(Fa);
        Fa = Fa / max(Fa);

        %% Visualization
        subplot(2,2,(ang-1)*2 + 1)
        polar(theta, -Fa, 'k'); hold on;
        title(['Array Factor - ' num2str(az) ' elements, \theta_0 = ' num2str(angle_deg(ang)) '°']);

        subplot(2,2,(ang-1)*2 + 2)
        pcolor(x/max(max(x)), y/max(max(y)), E); shading interp;
        pbaspect([1 1 1]); colormap jet;
        title(['E(x,y) - ' num2str(az) ' elements, \theta_0 = ' num2str(angle_deg(ang)) '°']);
        xlabel('x'); ylabel('y');
    end

    sgtitle(['Linear Array with ' num2str(az) ' Elements – Comparison of Steering Angles']);
end

%% Advanced Antenna Designs
% Assignment 2 - Exercise 1c - June 2025
% KARATIS DIMITRIOS 10775

clear all; clc; close all;

%% Simulation parameters
freq = 1e9;                 % Frequency [Hz]
c = 3e8;                    % Speed of light [m/s]
lambda = c/freq;            % Wavelength [m]
T = 1/freq;                 % Period [s]
omega = 2*pi*freq;          % Angular frequency [rad/s]
k = 2*pi/lambda;            % Wave number [rad/m]
az = 7;                     % Number of elements in the array

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
d = lambda/2;

% Color for annotation
teal = [ 0 0.5 0.5];

% Beam steering angle
theta0_deg = 30;
theta0 = theta0_deg * pi / 180;
delta = -k * d * sin(theta0);

%% Binomial amplitude distribution (normalized)
binomial_weights = zeros(1, az);
for i = 0:az-1
    binomial_weights(i+1) = nchoosek(az-1, i);
end
A = binomial_weights / max(binomial_weights);

%% Element positions (symmetric linear array)
r1x = -3*d; r2x = -2*d; r3x = -d; r4x = 0;
r5x = d;    r6x = 2*d;  r7x = 3*d;
r1y = 0; r2y = 0; r3y = 0; r4y = 0;
r5y = 0; r6y = 0; r7y = 0;

%% Compute total field E(x,y)
for ix=1:length(R)
    for iy=1:length(theta)
        R1 = sqrt((x(ix,iy)-r1x)^2 + (y(ix,iy)-r1y)^2);
        E1(ix,iy) = A(1) * cos(omega * t(1) - k*R1 - delta*-3);

        R2 = sqrt((x(ix,iy)-r2x)^2 + (y(ix,iy)-r2y)^2);
        E2(ix,iy) = A(2) * cos(omega * t(1) - k*R2 - delta*-2);

        R3 = sqrt((x(ix,iy)-r3x)^2 + (y(ix,iy)-r3y)^2);
        E3(ix,iy) = A(3) * cos(omega * t(1) - k*R3 - delta*-1);

        R4 = sqrt((x(ix,iy)-r4x)^2 + (y(ix,iy)-r4y)^2);
        E4(ix,iy) = A(4) * cos(omega * t(1) - k*R4 - delta*0);

        R5 = sqrt((x(ix,iy)-r5x)^2 + (y(ix,iy)-r5y)^2);
        E5(ix,iy) = A(5) * cos(omega * t(1) - k*R5 - delta*1);

        R6 = sqrt((x(ix,iy)-r6x)^2 + (y(ix,iy)-r6y)^2);
        E6(ix,iy) = A(6) * cos(omega * t(1) - k*R6 - delta*2);

        R7 = sqrt((x(ix,iy)-r7x)^2 + (y(ix,iy)-r7y)^2);
        E7(ix,iy) = A(7) * cos(omega * t(1) - k*R7 - delta*3);
    end
end

% Total field
E = E1 + E2 + E3 + E4 + E5 + E6 + E7;

%% Compute array factor with binomial weights
Fa = zeros(1, length(theta));
for i = 0:(az-1)
    Fa = Fa + (A(i+1) * exp(1i*i*delta + 1i*k*(i*d - ((az-1)/2)*d)*cos(theta)));
end
Fa = abs(Fa);
Fa = Fa / max(Fa);

%% Visualization
f1 = figure(10); clf; set(gcf,'Color',[1 1 1]); Fs=10;
sp1 = subplot(1,2,1); set(gca,'FontSize',Fs);
polar(theta, -Fa, 'k'); hold on;
title(['Array Factor (Binomial) – \theta_0 = ' num2str(theta0_deg) '°']);

subplot(1,2,2)
pcolor(x/max(max(x)), y/max(max(y)), E); shading interp;
pbaspect([1 1 1]);
colormap jet;
xlabel('x'); ylabel('y');
title('Spatial Field E(x,y) with Binomial Amplitude');

sgtitle('7-Element Linear Array with Binomial Weighting');

%% Advanced Antenna Designs
% Assignment 2 - Exercise 1a - June 2025
% KARATIS DIMITRIOS 10775

clear all; clc; close all;

%% Simulation parameters
freq = 1e9;               % Operating frequency [Hz]
c = 3e8;                  % Speed of light in free space [m/s]
lambda = c/freq;          % Wavelength [m]
T = 1/freq;               % Period of the wave [s]
omega = 2*pi*freq;        % Angular frequency [rad/s]
k = 2*pi/lambda;          % Wave number [rad/m]
az = 7;                   % Number of antenna elements in the array

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

% Coordinate matrices (for full 2D field representation)
x = R.'*cos(theta);            % x-coordinates of the field
y = R.'*sin(theta);            % y-coordinates of the field

% Custom color for annotations
teal = [0 0.5 0.5];

% Beam angles and element spacings
angle_deg = [0, 30, 60, 90];
d_values = [lambda/2, lambda/4];

% Loop over element spacings
for d_idx = 1:length(d_values)
    d = d_values(d_idx);

    % Group angles two-by-two per figure
    for fig_group = 1:2
        figure('Color','w');

        for subidx = 1:2
            ang_idx = (fig_group - 1)*2 + subidx;
            theta0 = angle_deg(ang_idx) * pi/180;
            delta = -k * d * sin(theta0);
            it = 1;

            % Array geometry (centered, symmetric)
            r1x = -3*d; r2x = -2*d; r3x = -d; r4x = 0;
            r5x = d;    r6x = 2*d;  r7x = 3*d;
            r1y = 0; r2y = 0; r3y = 0; r4y = 0;
            r5y = 0; r6y = 0; r7y = 0;

            % Field computation at each (x, y)
            for ix = 1:length(R)
                for iy = 1:length(theta)
                    R1 = sqrt((x(ix,iy)-r1x)^2 + (y(ix,iy)-r1y)^2);
                    E1(ix,iy) = cos(omega*t(it) - k*R1 - delta*-3);

                    R2 = sqrt((x(ix,iy)-r2x)^2 + (y(ix,iy)-r2y)^2);
                    E2(ix,iy) = cos(omega*t(it) - k*R2 - delta*-2);

                    R3 = sqrt((x(ix,iy)-r3x)^2 + (y(ix,iy)-r3y)^2);
                    E3(ix,iy) = cos(omega*t(it) - k*R3 - delta*-1);

                    R4 = sqrt((x(ix,iy)-r4x)^2 + (y(ix,iy)-r4y)^2);
                    E4(ix,iy) = cos(omega*t(it) - k*R4 - delta*0);

                    R5 = sqrt((x(ix,iy)-r5x)^2 + (y(ix,iy)-r5y)^2);
                    E5(ix,iy) = cos(omega*t(it) - k*R5 - delta*1);

                    R6 = sqrt((x(ix,iy)-r6x)^2 + (y(ix,iy)-r6y)^2);
                    E6(ix,iy) = cos(omega*t(it) - k*R6 - delta*2);

                    R7 = sqrt((x(ix,iy)-r7x)^2 + (y(ix,iy)-r7y)^2);
                    E7(ix,iy) = cos(omega*t(it) - k*R7 - delta*3);
                end
            end

            E = E1 + E2 + E3 + E4 + E5 + E6 + E7;

            % Array factor calculation
            A = ones(1, az);
            Fa = zeros(1, length(theta));
            for i = 0:(length(A)-1)
                Fa = Fa + (A(i+1) * exp(1i*i*delta + 1i*k*(i*d - ((length(A)-1)/2)*d)*cos(theta)));
            end
            Fa = abs(Fa);
            Fa = Fa / max(Fa);
            
            %% Visualization
            subplot(2,2,(subidx-1)*2 + 1)
            polar(theta, -Fa, 'k'); hold on;
            title(['Array Factor: \theta_0 = ' num2str(angle_deg(ang_idx)) '°, d = \lambda/' num2str(lambda/d)]);

            subplot(2,2,(subidx-1)*2 + 2)
            pcolor(x/max(max(x)), y/max(max(y)), E); shading interp;
            pbaspect([1 1 1]);
            title(['E(x,y): \theta_0 = ' num2str(angle_deg(ang_idx)) '°, d = \lambda/' num2str(lambda/d)]);
            colormap jet;
            xlabel('x'); ylabel('y');
        end

        sgtitle(['7-element Array – d = \lambda/' num2str(lambda/d) ...
                 ' – Angles: ' num2str(angle_deg((fig_group-1)*2+1)) '° & ' ...
                 num2str(angle_deg((fig_group-1)*2+2)) '°']);
    end
end

%% Disk–Gap–Band (DGB) Parachute Visualization
clear; close all; clc;

%% Inputs (all in inches for plotting)
Rdisk = 20;     % disk radius
Rgap  = 10;     % gap width
Rband = 30;     % band width
N     = 35;     % number of contour levels
Zmax  = Rdisk + Rband;   % canopy height scale (roughly)

%% Derived dimensions
Rtot = Rdisk + Rgap + Rband;  % total radius

%% Generate heights
zvals = linspace(Zmax,0,N);
colors = lines(N);

figure(1); hold on; grid on;
axis equal vis3d; view(45,25);
title("Disk–Gap–Band Parachute (wireframe rings)");
xlabel("X (in)"); ylabel("Y (in)"); zlabel("Z (in)");

for k = 1:N
    z = zvals(k);

    % Dome scaling (spherical shrink)
    scale = sqrt(max(0, 1 - (z/Zmax)^2));

    % Effective radius at this height
    Reff = Rtot * scale;

    % If inside disk radius: plot solid disk
    if Reff <= Rdisk
        theta = linspace(0,2*pi,200);
        xring = Reff*cos(theta);
        yring = Reff*sin(theta);
        plot3(xring, yring, z*ones(size(xring)), 'Color', colors(k,:), 'LineWidth', 1.2);

    % If inside gap: skip plotting (open space)
    elseif Reff <= Rdisk + Rgap
        continue;

    % If inside band: plot annulus
    else
        theta = linspace(0,2*pi,200);
        xring_outer = Reff*cos(theta);
        yring_outer = Reff*sin(theta);
        plot3(xring_outer, yring_outer, z*ones(size(xring_outer)), ...
              'Color', colors(k,:), 'LineWidth', 1.2);

        % inner edge of the band (gap boundary)
        xring_inner = (Rdisk+Rgap)*scale*cos(theta);
        yring_inner = (Rdisk+Rgap)*scale*sin(theta);
        plot3(xring_inner, yring_inner, z*ones(size(xring_inner)), ...
              'Color', colors(k,:), 'LineWidth', 1.2);
    end
end

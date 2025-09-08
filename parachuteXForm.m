%% X-form Parachute Model with 3D Visualization (no shroud lines)
clear; close all; clc;

%% Inputs
Cd = .97;             % Drag coefficient
p  = 1.0;             % Air density kg/m^3 (~6000 ft)
m  = 36/2.2;          % Mass (lb to kg)
g  = 9.81;            % Gravity

v_target = 100 * 12*2.54/100;  % Target descent rate (ft/s → m/s)

%% Required Area
S = 2*m*g / (p*(v_target^2)*Cd);  % Required area (m^2)

%% X-form Geometry
AR = 2;                     % Aspect ratio L/w
w  = sqrt(S / (2*AR - 1));  % width (m)
L  = AR * w;                % length (m)

% Convert to inches for construction
w_in = w*100/2.54;
L_in = L*100/2.54;

fprintf("X-form parachute design:\n");
fprintf("Target descent rate: %.1f ft/s\n", v_target*100/(12*2.54));
fprintf("Projected area needed: %.2f m^2\n", S);
fprintf("Rectangle width w = %.1f in\n", w_in);
fprintf("Rectangle length L = %.1f in\n", L_in);

%% Shroud Line Length (calc only, not plotted)
numLines = 8; % typical (2 per arm)
armTip = L_in/2; % from center to tip
Ls = numLines*(armTip/36 + 0.25); % yards (simple approx)
fprintf("Shroud line length needed ≈ %.1f yds total\n", Ls);

%% Pattern Plot (flat cruciform)
figure(1); hold on; axis equal; grid on;
% Vertical arm
x = [-w_in/2 w_in/2 w_in/2 -w_in/2];
y = [-L_in/2 -L_in/2 L_in/2 L_in/2];
fill(x,y,[0.8 0.8 1]);

% Horizontal arm
x = [-L_in/2 -L_in/2 L_in/2 L_in/2];
y = [-w_in/2 w_in/2 w_in/2 -w_in/2];
fill(x,y,[0.6 0.6 1]);

title("X-form Parachute Flat Pattern");
xlabel("in"); ylabel("in");

%% X-form 3D Visualization (wireframe, smooth rounded apex)
figure(5); clf; hold on; grid on;
axis equal vis3d; view(45,25);
title("X-form Parachute Wireframe (rounded apex)");
xlabel("X (in)"); ylabel("Y (in)"); zlabel("Z (in)");

if ~exist('L_in','var') || ~exist('w_in','var')
    L_in = 28;   % total arm length (inches)
    w_in = 14;   % total arm width (inches)
end

N     = 35;                    % number of contour rings
Zmax  = L_in/2;                % dome "height"
zvals = linspace(Zmax,0,N);    % top → base
colors = lines(N);

for k = 1:N
    z = zvals(k);

    % Shrink factor from spherical dome
    scale = sqrt(max(0, 1 - (z/Zmax)^2));

    % Arm half-width (constant sideways)
    w_now = (w_in/2);

    % Arm half-length shrinks with height
    L_now = (L_in/2) * scale;

    % --- Smooth apex handling ---
    if L_now <= w_now
        % Collapse to a small circle cap instead of flat square
        theta = linspace(0,2*pi,100);
        rcap = w_now * (L_now / (w_in/2)); % taper radius to zero
        xring = rcap * cos(theta);
        yring = rcap * sin(theta);
    else
        % Normal cruciform outline
        [xring, yring] = cruciform_outline(w_now, L_now);
    end

    % Plot ring at height z
    plot3(xring, yring, z*ones(size(xring)), ...
          'Color', colors(k,:), 'LineWidth', 1.2);
end

%% Local function for cruciform perimeter
function [xring, yring] = cruciform_outline(w, L)
    % outline of cruciform = vertical arm + horizontal arm
    xring = [-w,  +w,  +w,  +L,  +L,  +w,  +w,  -w,  -w,  -L,  -L,  -w,  -w];
    yring = [-L,  -L,  -w,  -w,  +w,  +w,  +L,  +L,  +w,  +w,  -w,  -w,  -L];
end

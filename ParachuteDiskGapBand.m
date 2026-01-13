%% Disk–Gap–Band (DGB) Parachute Visualization + Sewing Pattern (One Gore)
% Includes a center spill hole (vent).
clear; close all; clc;

%% ===================== USER INPUTS (inches) =====================
Rdisk = 6;     % disk radius
Rgap  = 1;     % gap width (open space)
Rband = 1.5;     % band width (fabric from gap edge to skirt)

Rspill = 1;     % SPILL HOLE radius (in)  <-- set to 0 for no spillhole

numS = 16;      % number of gores / chute sections
seam = 0.5;     % seam allowance (in)

% Plot resolution
dr = 0.25;      % radial step for 3D rings (smaller = smoother)
M  = 1200;      % resolution for sewing pattern unwrapping

%% ===================== BASIC CHECKS =====================
if Rspill < 0
    error("Rspill must be >= 0.");
end
if Rspill >= Rdisk
    error("Rspill must be smaller than Rdisk (spill hole must fit inside the disk).");
end

%% ===================== DERIVED DIMENSIONS =====================
Rtot = Rdisk + Rgap + Rband;   % total radius at skirt
Zmax = Rdisk + Rband;          % height scale (kept from your Script 1 idea)
A    = Rtot;                   % base radius for profile
b    = Zmax;                   % max height for profile

% Canopy profile (similar look to your Script 2 rings):
% z(r) = b * sqrt( 1 - (r/A)^2 ),  r in [0, A]
z_of_r = @(r) b .* sqrt(max(0, 1 - (r./A).^2));

%% ===================== FIGURE 1: 3D DGB WIREFRAME RINGS =====================
figure(1); clf; hold on; grid on;
axis equal; view(45,25);
title("Disk–Gap–Band Parachute (wireframe rings) with spill hole");
xlabel("X (in)"); ylabel("Y (in)"); zlabel("Z (in)");

theta = linspace(0, 2*pi, 400);

for r = 0:dr:A
    z = z_of_r(r);

    % --- SPILLHOLE: skip anything inside the vent radius ---
    if r < Rspill
        continue;
    end

    if r <= Rdisk
        % Disk region: draw ring at radius r
        X = r*cos(theta);
        Y = r*sin(theta);
        plot3(X, Y, z*ones(size(X)), 'LineWidth', 1.2);

    elseif r <= (Rdisk + Rgap)
        % Gap region: skip drawing (open space)
        continue;

    else
        % Band region: draw ring at radius r
        Xo = r*cos(theta);
        Yo = r*sin(theta);
        plot3(Xo, Yo, z*ones(size(Xo)), 'LineWidth', 1.2);

        % Also draw the inner edge of the band (gap boundary) for reference
        rin = (Rdisk + Rgap);
        zi  = z_of_r(rin);
        Xi = rin*cos(theta);
        Yi = rin*sin(theta);
        plot3(Xi, Yi, zi*ones(size(Xi)), 'LineWidth', 1.2);
    end
end

% Draw the spillhole ring explicitly (nice visual)
if Rspill > 0
    zspill = z_of_r(Rspill);
    Xs = Rspill*cos(theta);
    Ys = Rspill*sin(theta);
    plot3(Xs, Ys, zspill*ones(size(Xs)), 'LineWidth', 2.0);
end

%% ===================== FIGURE 2: SEWING PATTERN (ONE GORE) =====================
% Unwrap into meridional distance y and gore half-width x(y).
% Gap means: no fabric between r=Rdisk and r=Rdisk+Rgap.
% Spill hole means: no fabric between r=0 and r=Rspill (disk starts at Rspill).

% Create r grid
r = linspace(0, A, M);
z = z_of_r(r);

% dz/dr numerically
dzdr = gradient(z, r);

% Meridional distance from apex
dsdr = sqrt(1 + dzdr.^2);
y = cumtrapz(r, dsdr);     % y(0)=0 at apex

% Gore half-width
xHalf = (pi/numS) .* r;

% Masks for fabric pieces
maskDisk = (r >= Rspill) & (r <= Rdisk);                 % disk fabric starts at spillhole
maskBand = (r >= (Rdisk + Rgap)) & (r <= Rtot);          % band fabric

% Extract curves
yDisk = y(maskDisk);  xDisk = xHalf(maskDisk);
yBand = y(maskBand);  xBand = xHalf(maskBand);

% Seam allowance (simple x expansion)
xDisk_seam = xDisk + seam;
xBand_seam = xBand + seam;

% Useful y-locations for cut boundaries
ySpill     = interp1(r, y, Rspill, 'linear', 'extrap');
yDiskEnd   = interp1(r, y, Rdisk, 'linear', 'extrap');
yBandStart = interp1(r, y, (Rdisk + Rgap), 'linear', 'extrap');
ySkirt     = y(end);

figure(2); clf; hold on; axis equal; grid on;
title("DGB Sewing Pattern (One Gore) with spill hole");
xlabel("Width (in)"); ylabel("Meridional length from apex (in)");

% ---- Disk gore outline ----
plot( xDisk, yDisk, 'k', 'LineWidth', 2);
plot(-xDisk, yDisk, 'k', 'LineWidth', 2);
plot([0 0], [ySpill yDiskEnd], '--k');   % centerline only across disk fabric region

% seam allowance (dotted)
plot( xDisk_seam, yDisk, ':k', 'LineWidth', 1.5);
plot(-xDisk_seam, yDisk, ':k', 'LineWidth', 1.5);

% ---- Band gore outline ----
plot( xBand, yBand, 'b', 'LineWidth', 2);
plot(-xBand, yBand, 'b', 'LineWidth', 2);
plot( xBand_seam, yBand, ':b', 'LineWidth', 1.5);
plot(-xBand_seam, yBand, ':b', 'LineWidth', 1.5);

% ---- Mark NO-FABRIC regions ----
% Spill hole region marker
yline(ySpill, '--m', 'Spill hole edge / cut');
text(0, ySpill/2, "SPILL HOLE (no fabric)", ...
    'HorizontalAlignment','center', 'Color','m', 'FontWeight','bold');

% Gap region marker
yline(yDiskEnd,   '--r', 'Disk edge / cut');
yline(yBandStart, '--r', 'Band starts / cut');
text(0, (yDiskEnd + yBandStart)/2, "GAP (no fabric)", ...
    'HorizontalAlignment','center', 'Color','r', 'FontWeight','bold');

% Annotation: skirt width for one gore
goreSkirtWidth = 2*(pi/numS)*Rtot;
text(0, ySkirt, sprintf("Skirt width of one gore: %.2f in", goreSkirtWidth), ...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom');

% Limits
xMax = max(xHalf) + 2*seam;
xlim([-xMax, xMax]);
ylim([0, ySkirt*1.05]);

%% ===================== OPTIONAL: EXPORT PATTERN TO PDF =====================
% figure(2);
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8.5 11]); % letter
% print('-dpdf','-r300','DGB_gore_pattern_with_spillhole_letter.pdf');

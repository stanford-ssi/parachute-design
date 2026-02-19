%% DGB Parachute: 3D Visualization + Disk/Band printable PDF pages (seam-length-controlled)
% Figure 1: 3D DGB visualization (disk-gap-band + spill hole)
% Disk pages: DISK_PAGES/DISK_page_XX.pdf
% Band pages: BAND_PAGES/BAND_page_XX.pdf
%
% CHANGE: You specify desired SEAM LENGTHS along the fabric (meridional arc length)
% for disk, gap, and band, and the script solves for radii that match those seam lengths.
%
% yComp scales ALL Y geometry INCLUDING the calibration square.
% Export/tiling windows are also in compensated coordinates.

clear; close all; clc;

%% ===================== USER INPUTS (inches) =====================
% Profile height scale (controls dome "depth")
Zmax = 14;              % inches

% Desired seam lengths along fabric (inches)
LdiskSeam = 12;         % spill edge -> disk edge (fabric)
LgapSeam  = 3;          % disk edge -> band start (no fabric)
LbandSeam = 4;          % band start -> skirt (fabric)

% Spill hole radius (plan-view inches)
Rspill = 2;

% Gores / panels
numS = 8;
seam = 0.5;

% Resolutions
M = 4000;               % arc-length mapping resolution
dr3D = 0.25;            % 3D ring step
nth3D = 400;            % 3D theta resolution

%% ===================== PRINT SETTINGS =====================
pageW  = 8.5;
pageH  = 11.0;
margin = 0.5;
dpi    = 300;

diskDir = 'DISK_PAGES';
bandDir = 'BAND_PAGES';

%% ===================== PRINT COMPENSATION =====================
yComp = 1/1.125;     % change if your printer stretches Y
xComp = 1;

%% ===================== SOLVE RADII FROM SEAM LENGTHS =====================
Rtot_guess = max(Rspill + (LdiskSeam + LgapSeam + LbandSeam), Rspill + 1);

for iter = 1:30
    A = Rtot_guess;
    z_of_r = @(rr) Zmax .* sqrt(max(0, 1 - (rr./A).^2));

    r = linspace(0, A, M);
    z = z_of_r(r);
    dzdr = gradient(z, r);
    dsdr = sqrt(1 + dzdr.^2);
    s = cumtrapz(r, dsdr);   % meridional distance from apex

    if Rspill >= A
        error('Rspill must be smaller than the final skirt radius Rtot. Increase guess/seam lengths.');
    end

    sAt = @(rq) interp1(r, s, rq, 'linear', 'extrap');
    s_spill = sAt(Rspill);

    s_disk_edge   = s_spill + LdiskSeam;
    s_band_start  = s_disk_edge + LgapSeam;
    s_skirt       = s_band_start + LbandSeam;

    Rdisk_new      = interp1(s, r, s_disk_edge,  'linear', 'extrap');
    RbandStart_new = interp1(s, r, s_band_start, 'linear', 'extrap');
    Rtot_new       = interp1(s, r, s_skirt,      'linear', 'extrap');

    if ~isfinite(Rtot_new) || Rtot_new <= 0
        error('Failed to solve radii. Check Zmax and seam lengths.');
    end

    if abs(Rtot_new - Rtot_guess) < 1e-4
        Rdisk = Rdisk_new;
        Rtot  = Rtot_new;
        Rgap  = max(0, RbandStart_new - Rdisk_new);
        Rband = max(0, Rtot_new - RbandStart_new);
        break;
    end

    Rtot_guess = 0.6*Rtot_guess + 0.4*Rtot_new;

    if iter == 30
        Rdisk = Rdisk_new;
        Rtot  = Rtot_new;
        Rgap  = max(0, RbandStart_new - Rdisk_new);
        Rband = max(0, Rtot_new - RbandStart_new);
    end
end

RbandStart = Rdisk + Rgap;

%% ===================== BASIC CHECKS =====================
if Rspill < 0, error('Rspill must be >= 0.'); end
if Rspill >= Rdisk, error('Rspill must be smaller than Rdisk.'); end
if RbandStart > Rtot, error('Band start exceeds skirt radius. Check seam lengths.'); end

fprintf('\nSolved radii from seam lengths:\n');
fprintf('  Rdisk      = %.3f in\n', Rdisk);
fprintf('  Rgap       = %.3f in\n', Rgap);
fprintf('  Rband      = %.3f in\n', Rband);
fprintf('  Rtot       = %.3f in\n', Rtot);
fprintf('  Band start = %.3f in\n', RbandStart);

%% ===================== 3D VISUALIZATION (FIGURE 1) =====================
A = Rtot;
z_of_r = @(rr) Zmax .* sqrt(max(0, 1 - (rr./A).^2));

figure(1); clf; hold on; grid on;
axis equal; view(45,25);
title('DGB Parachute 3D Visualization (disk-gap-band + spill hole)');
xlabel('X (in)'); ylabel('Y (in)'); zlabel('Z (in)');

theta = linspace(0, 2*pi, nth3D);

for rr = 0:dr3D:A
    if rr < Rspill
        continue; % spill hole (open)
    end

    if rr <= Rdisk
        % Disk fabric region
        zz = z_of_r(rr);
        plot3(rr*cos(theta), rr*sin(theta), zz*ones(size(theta)), 'k', 'LineWidth', 1.0);

    elseif rr <= (Rdisk + Rgap)
        % Gap region (open)
        continue;

    else
        % Band fabric region
        zz = z_of_r(rr);
        plot3(rr*cos(theta), rr*sin(theta), zz*ones(size(theta)), 'b', 'LineWidth', 1.0);

        % Draw inner band edge for reference (gap boundary)
        rin = (Rdisk + Rgap);
        zin = z_of_r(rin);
        plot3(rin*cos(theta), rin*sin(theta), zin*ones(size(theta)), 'b', 'LineWidth', 1.0);
    end
end

% Draw spill edge ring for reference
if Rspill > 0
    zsp = z_of_r(Rspill);
    plot3(Rspill*cos(theta), Rspill*sin(theta), zsp*ones(size(theta)), 'm', 'LineWidth', 1.5);
end

%% ===================== FINAL UNWRAP USING SOLVED RADII =====================
r = linspace(0, A, M);
z = z_of_r(r);

dzdr = gradient(z, r);
dsdr = sqrt(1 + dzdr.^2);
y = cumtrapz(r, dsdr);
xHalf = (pi/numS).*r;

% Disk (fabric)
maskDisk = (r >= Rspill) & (r <= Rdisk);
xDisk  = xHalf(maskDisk);
yDisk  = y(maskDisk);
xDiskS = xDisk + seam;

% Band (fabric)
maskBand = (r >= RbandStart) & (r <= Rtot);
xBand  = xHalf(maskBand);
yBand  = y(maskBand);
xBandS = xBand + seam;

ySpill     = interp1(r, y, Rspill,     'linear', 'extrap');
yBandStart = interp1(r, y, RbandStart, 'linear', 'extrap');

%% ===================== TILE SETTINGS (inches) =====================
tileW = pageW - 2*margin;
tileH = pageH - 2*margin;

diskXLim = [-max(xDiskS), +max(xDiskS)];
diskYLim = [0, max(yDisk)*1.05];

bandXLim = [-max(xBandS), +max(xBandS)];
bandYLim = [yBandStart, max(yBand)*1.05];

%% ===================== EXPORT DISK =====================
export_pages_newfig_TRUE_SCALE_compSpace( ...
    diskDir, 'DISK', pageW, pageH, margin, dpi, ...
    diskXLim, diskYLim, tileW, tileH, xComp, yComp, ...
    @(ax, xL, yL) draw_disk_comp(ax, xDisk, yDisk, xDiskS, yDisk, Rspill, numS, ySpill, xComp, yComp, xL, yL) );

%% ===================== EXPORT BAND =====================
export_pages_newfig_TRUE_SCALE_compSpace( ...
    bandDir, 'BAND', pageW, pageH, margin, dpi, ...
    bandXLim, bandYLim, tileW, tileH, xComp, yComp, ...
    @(ax, xL, yL) draw_band_comp(ax, xBand, yBand, xBandS, yBand, yBandStart, xComp, yComp, xL, yL) );

fprintf('\nDone.\nDisk pages: %s/\nBand pages: %s/\n', diskDir, bandDir);
fprintf('Using compensation: xComp=%.6f, yComp=%.6f\n', xComp, yComp);

%% ========================================================================
%% Export helper: NEW figure per page; tiling + axes limits in COMPENSATED space
function export_pages_newfig_TRUE_SCALE_compSpace(outDir, baseName, pageW, pageH, margin, dpi, ...
                                                  fullXLim, fullYLim, tileW, tileH, xComp, yComp, drawFcn)
    if ~exist(outDir,'dir'), mkdir(outDir); end

    fullXLimC = fullXLim * xComp;
    fullYLimC = fullYLim * yComp;

    fullW = diff(fullXLimC);
    fullH = diff(fullYLimC);
    nX = max(1, ceil(fullW / tileW));
    nY = max(1, ceil(fullH / tileH));

    pageIdx = 1;

    for iy = 1:nY
        for ix = 1:nX
            x0 = fullXLimC(1) + (ix-1)*tileW;
            y0 = fullYLimC(1) + (iy-1)*tileH;

            xL = [x0, x0 + tileW];
            yL = [y0, y0 + tileH];

            f = figure('Visible','off');
            set(f, 'Units','inches', 'Position',[1 1 pageW pageH]);
            set(f, 'PaperUnits','inches', 'PaperSize',[pageW pageH]);
            set(f, 'PaperPositionMode','manual', 'PaperPosition',[0 0 pageW pageH]);

            ax = axes(f, 'Units','inches', 'Position',[margin margin tileW tileH]);
            set(ax, 'ActivePositionProperty','position');

            xlim(ax, xL);
            ylim(ax, yL);

            axis(ax,'off');
            axis(ax,'normal');

            drawFcn(ax, xL, yL);
            drawnow;

            fname = fullfile(outDir, sprintf('%s_page_%02d.pdf', baseName, pageIdx));
            print(f, fname, '-dpdf', sprintf('-r%d', dpi), '-painters');
            close(f);

            pageIdx = pageIdx + 1;
        end
    end
end

%% ========================================================================
%% Disk draw in compensated space
function draw_disk_comp(ax, x, y, xS, yS, Rspill, numS, ySpill, xComp, yComp, xL, yL)
    hold(ax,'on');
    set(ax,'Visible','off');

    X  = x  * xComp;   Y  = y  * yComp;
    XS = xS * xComp;   YS = yS * yComp;

    plot(ax,  X,  Y,  'k', 'LineWidth', 2);
    plot(ax, -X,  Y,  'k', 'LineWidth', 2);
    plot(ax,  XS, YS, ':k','LineWidth', 1.5);
    plot(ax, -XS, YS, ':k','LineWidth', 1.5);

    if Rspill > 0
        xSp = (pi/numS)*Rspill * xComp;
        ySp = ySpill * yComp;
        plot(ax, [-xSp xSp], [ySp ySp], '--m', 'LineWidth', 1.5);
    end

    cx = xL(1) + 0.5;
    cy = yL(1) + 0.5;
    w  = 1 * xComp;
    h  = 1 * yComp;
    plot(ax, [cx cx+w cx+w cx cx], [cy cy cy+h cy+h cy], 'r-', 'LineWidth', 3);
end

%% ========================================================================
%% Band draw in compensated space
function draw_band_comp(ax, x, y, xS, yS, yBandStart, xComp, yComp, xL, yL)
    hold(ax,'on');
    set(ax,'Visible','off');

    X  = x  * xComp;   Y  = y  * yComp;
    XS = xS * xComp;   YS = yS * yComp;
    yStart = yBandStart * yComp;

    plot(ax,  X,  Y,  'b', 'LineWidth', 2);
    plot(ax, -X,  Y,  'b', 'LineWidth', 2);
    plot(ax,  XS, YS, ':b','LineWidth', 1.5);
    plot(ax, -XS, YS, ':b','LineWidth', 1.5);

    plot(ax, [-max(X) max(X)], [yStart yStart], '--k', 'LineWidth', 1.2);

    cx = xL(1) + 0.5;
    cy = yL(1) + 0.5;
    w  = 1 * xComp;
    h  = 1 * yComp;
    plot(ax, [cx cx+w cx+w cx cx], [cy cy cy+h cy+h cy], 'r-', 'LineWidth', 3);
end

%% ===================== DESCENT RATE ESTIMATE =====================
% Inputs
m_lb = 50;              % mass of rocket + recovery gear attached (lbm)
Cd   = .55;             % effective drag coefficient (typical ballpark: 1.5–2.2)
rho  = 1.00;            % air density (kg/m^3). ~1.225 sea level; ~0.95 around 6000 ft

% Unit conversions
m = m_lb * 0.45359237;  % kg
g = 9.80665;            % m/s^2

% Geometry (inches -> meters)
in2m = 0.0254;
Rtot_m   = Rtot   * in2m;
Rdisk_m  = Rdisk  * in2m;
Rspill_m = Rspill * in2m;
RbandStart_m = (Rdisk + Rgap) * in2m;   % start radius of band (outer edge of gap)

% Projected areas (m^2)
A_total = pi * (Rtot_m^2);

% "Open" areas: spill hole + gap annulus (no fabric)
A_spill = pi * (Rspill_m^2);
A_gap   = pi * (RbandStart_m^2 - Rdisk_m^2);  % annulus between disk edge and band start

% Effective projected area (simple first-order model)
% Clamp to avoid negative if you choose odd dimensions
A_eff = max(0.05*A_total, A_total - A_spill - A_gap);

% Terminal descent rate (m/s)
v_term = sqrt( (2*m*g) / (rho * Cd * A_eff) );

% Output
v_fts = v_term / 0.3048;

fprintf('\n=== Descent Rate Estimate ===\n');
fprintf('Mass: %.2f lbm (%.2f kg)\n', m_lb, m);
fprintf('Cd: %.2f, rho: %.3f kg/m^3\n', Cd, rho);
fprintf('Rtot: %.2f in, Rdisk: %.2f in, Rgap: %.2f in, Rspill: %.2f in\n', Rtot, Rdisk, Rgap, Rspill);
fprintf('A_total: %.4f m^2\n', A_total);
fprintf('A_spill: %.4f m^2, A_gap: %.4f m^2\n', A_spill, A_gap);
fprintf('A_eff:   %.4f m^2 (%.1f%% of total)\n', A_eff, 100*(A_eff/A_total));
fprintf('Terminal descent rate: %.2f m/s (%.1f ft/s)\n', v_term, v_fts);

clc;
clear;
close all;

%% 1) Paper dimensions (millimetres)
a       = 8.8;
b       = 8.8;
y       = 7.4;       % SRR outer width
z       = 7.7;       % SRR outer height
c       = 0.7;       % copper trace width
g       = 0.3;       % SRR inner–outer gap and split width
monoBody = 8.8;      % monopole height above ground
lg_mm    = 10;       % ground pad height (from paper)
monoW_mm = 3;        % monopole width
wg_mm    = 26.6;     % ground pad width (from paper)
boardW   = 26.6;     % PCB width
boardL   = 21.8;     % PCB length

monoLen_mm = monoBody + lg_mm;  % full monopole including base
gapX = monoW_mm/2 + boardW/2 + (a - y)/2;
gapY = lg_mm + (b - z)/2;

%% 2) Convert to metres
mm = 1e-3;
a        = a       * mm;
b        = b       * mm;
z        = z       * mm;
y        = y       * mm;
c        = c       * mm;
g        = g       * mm;
lg       = lg_mm   * mm;
monoLen  = monoLen_mm * mm;
monoW    = monoW_mm * mm;
wg       = wg_mm   * mm;
boardW   = boardW  * mm;
boardL   = boardL  * mm;
gapX     = gapX    * mm;
gapY     = gapY    * mm;

%% 3) Substrate
hsub = 1.6e-3;   % substrate thickness
epsR = 4.4;

%% 4) SRR construction
outerFull = antenna.Rectangle('Center',[y/2 z/2],'Length',y,'Width',z);
outerHole = antenna.Rectangle('Center',[y/2 z/2],'Length',y-2*c,'Width',z-2*c);
outerRing = outerFull - outerHole;
outerRing = outerRing - antenna.Rectangle('Center',[y/2  c/2],'Length',g,'Width',c);

ai = y - 2*(c+g);  bi = z - 2*(c+g);
innerFull = antenna.Rectangle('Center',[y/2 z/2],'Length',ai,'Width',bi);
innerHole = antenna.Rectangle('Center',[y/2 z/2],'Length',ai-2*c,'Width',bi-2*c);
innerRing = innerFull - innerHole;
yTopInner = z - (c+g) - c/2;
innerRing = innerRing - antenna.Rectangle('Center',[y/2 yTopInner],'Length',g,'Width',c);

SRR = outerRing + innerRing;

figure('Name','SRR unit cell');
show(SRR); view(0,90); axis equal tight;
title('Double SRR');

%% 5a) Twin ground pads (coplanar layout from paper)
% Left pad: from x = 0 to x = (boardW - monoW) / 2
leftPad = antenna.Rectangle('Length', (boardW - monoW)/2, 'Width', lg, ...
    'Center', [(boardW - monoW)/4 , lg/2]);

% Right pad: from x = (boardW + monoW)/2 to x = boardW
rightPad = antenna.Rectangle('Length', (boardW - monoW)/2, 'Width', lg, ...
    'Center', [boardW - (boardW - monoW)/4 , lg/2]);

groundPads = leftPad + rightPad;

%% 5b) Monopole strip (from top of ground up)
mono = antenna.Rectangle('Center',[boardW/2 , monoLen/2], ...
                         'Length',monoW,'Width',monoLen);

% Debugging view: Show ground pads + monopole strip
figure('Name','Ground Pads + Monopole');
show(groundPads + mono);
axis equal; view(0,90);
title('Monopole strip + Twin Ground Pads');

%% 6) Shift SRR vertically to start from top of ground
SRRshift = translate(SRR, [gapX , gapY , 0]);

%% 7) Board outline and stack
board = antenna.Rectangle('Center',[boardW/2 , boardL/2], ...
                          'Length',boardW,'Width',boardL);
FR4 = dielectric('FR4');  FR4.EpsilonR = epsR;  FR4.Thickness = hsub;

pcb = pcbStack;
pcb.BoardShape     = board;
pcb.BoardThickness = hsub;
pcb.Layers         = { mono + SRRshift + groundPads, FR4 , board };

%% 8) Feed location at base center of monopole
pcb.FeedLocations = [boardW/2 , 0 , 3 , 1];
pcb.FeedDiameter = 1e-3;


%% 9) Geometry view
figure('Name','Planar ESA Geometry');
show(pcb); axis equal; view(0,90);
title('Planar SRR-monopole geometry');

%% 10) Meshing
figure('Name','Mesh settings');
mesh(pcb,'MaxEdgeLength',0.1,'MinEdgeLength',0.01,'GrowthRate',0.7);
title('\lambda/10 mesh for 1 – 8 GHz simulation');

%% 11) S11 simulation
freq = linspace(1e9,8e9,300);
S = sparameters(pcb,freq,50);

figure('Name','S11 response');
rfplot(S,1,1); grid on; ylim([-35 0]);
xlabel('Frequency (GHz)'); ylabel('|S_{11}| (dB)');

% 11b) S11 of plain monopole (no SRR)
pcb_noSRR = pcbStack;
pcb_noSRR.BoardShape = board;
pcb_noSRR.BoardThickness = hsub;
pcb_noSRR.Layers = { mono + groundPads , FR4 , board };  % no SRR
pcb_noSRR.FeedLocations = [boardW/2 , 0 , 3 , 1];
pcb_noSRR.FeedDiameter = 1e-3;

S_noSRR = sparameters(pcb_noSRR,freq,50);

% 11c) Plot both |S11|
figure('Name','S11 Comparison: With and Without SRR');
rfplot(S,1,1); hold on;
rfplot(S_noSRR,1,1);
legend('|S_{11}| with SRR','|S_{11}| without SRR');
grid on; ylim([-35 0]);
xlabel('Frequency (GHz)'); ylabel('|S_{11}| (dB)');
title('|S_{11}| Comparison');


%% 12) 3-D pattern at 2.45 GHz
figure('Name','3D Radiation Pattern @ 2.45 GHz');
pattern(pcb,2.45e9);

%% 13) Elevation cut at phi = 0°
figure('Name','Elevation φ = 0°, 2.45 GHz');
patternElevation(pcb, 2.45e9, 0); hold on;
legend('With SRR');
title('Elevation pattern (XZ plane) @ 2.45 GHz');

figure('Name','Elevation φ = 0°, 2.45 GHz');
patternElevation(pcb_noSRR, 2.45e9, 0); hold on;
legend('Without SRR');
title('Elevation pattern (XZ plane) @ 2.45 GHz');


S11_dB = 20*log10(abs(squeeze(S.Parameters(1,1,:))));
fGHz = freq / 1e9;

% Find indices where |S11| <= -10 dB
idx_below_10dB = find(S11_dB <= -10);

if ~isempty(idx_below_10dB)
    bw_start = fGHz(idx_below_10dB(1));
    bw_end   = fGHz(idx_below_10dB(end));
    bw_MHz   = (bw_end - bw_start) * 1e3;
    fprintf('With SRR: -10 dB bandwidth = %.2f – %.2f GHz (%.1f MHz)\n', ...
        bw_start, bw_end, bw_MHz);
else
    fprintf('With SRR: No -10 dB bandwidth found.\n');
end

S11_noSRR_dB = 20*log10(abs(squeeze(S_noSRR.Parameters(1,1,:))));
idx_below_10dB_noSRR = find(S11_noSRR_dB <= -10);

if ~isempty(idx_below_10dB_noSRR)
    bw_start2 = fGHz(idx_below_10dB_noSRR(1));
    bw_end2   = fGHz(idx_below_10dB_noSRR(end));
    bw2_MHz   = (bw_end2 - bw_start2) * 1e3;
    fprintf('Without SRR: -10 dB bandwidth = %.2f – %.2f GHz (%.1f MHz)\n', ...
        bw_start2, bw_end2, bw2_MHz);
else
    fprintf('Without SRR: No -10 dB bandwidth found.\n');
end

function eta = estimateRadiationEfficiency(pcb, freq)
    % Estimate max gain and directivity at given frequency
    G_max = pattern(pcb, freq, 0, 0, 'Type', 'gain');
    D_max = pattern(pcb, freq, 0, 0, 'Type', 'directivity');

    % Convert from dB to linear
    G_lin = 10^(G_max / 10);
    D_lin = 10^(D_max / 10);

    % Radiation efficiency is gain/directivity
    eta = G_lin / D_lin;
end

eta_srr = estimateRadiationEfficiency(pcb, 2.45e9);
fprintf('Estimated Efficiency with SRR: %.2f%%\n', eta_srr * 100);

eta_plain = estimateRadiationEfficiency(pcb_noSRR, 2.45e9);
fprintf('Estimated Efficiency without SRR: %.2f%%\n', eta_plain * 100);

% Compactness based on Chu Limit
lambda = 3e8 / 2.45e9;
a_sphere = sqrt((boardW/2)^2 + monoLen^2);  % estimate radius in metres
ka = 2 * pi * a_sphere / lambda;
Q_chu = 1 / ka^3 + 1 / ka;

Q_srr = 2.45e9 / (bw_MHz * 1e6);
Q_plain = 2.45e9 / (bw2_MHz * 1e6);


if eta_srr <= 0.001
    Q_eff_srr = 0;
else
    Q_eff_srr = Q_srr / eta_srr;
end

if eta_plain <= 0.001
    Q_eff_plain = 0;
else
    Q_eff_plain = Q_plain / eta_plain;
end


fprintf('\n--- Compactness Analysis (Chu Limit) ---\n');
fprintf('ka = %.4f\n', ka);
fprintf('Chu Lower Bound Q = %.2f\n', Q_chu);
fprintf('Effective Q (SRR) = %.2f\n', Q_eff_srr);
fprintf('Effective Q (No SRR) = %.2f\n', Q_eff_plain);
clear all; clc;

%% Simulation parameters
freq   = 1e9;      % Hz
c      = 3e8;      % speed of light
lambda = c/freq;   
T      = 1/freq;   
omega  = 2*pi*freq;
k      = 2*pi/lambda;

Ns     = 30;       % spatial samples per λ
ds     = lambda/Ns;
Nt     = 35;       % temporal samples per period
dt     = T/Nt;
t      = 0:dt:T;   % one period

R      = 0:ds:8*lambda;     % radial domain
Ntheta = 240;               % number of angles
dtheta = 2*pi/Ntheta;
theta  = 0:dtheta:2*pi;     % 0…360°

%% Build spatial grid
x = R.' * cos(theta);
y = R.' * sin(theta);

%% Plot colors
teal = [0 0.5 0.5];

%% Part Γ – Binomial taper on a 7-element array
az    = 7;              % choose 5, 7 or 9
d     = lambda/2;       % element spacing
mid   = (az - 1)/2;     

% Compute binomial weights (normalized)
idx = 0:(az-1);
w   = factorial(az-1) ./ ( factorial(idx) .* factorial((az-1)-idx) );
A   = w / sum(w);       

% Steering angles vector
deltas = [0, pi/3, pi/2, 2*pi/3];  % in radians

% Single time snapshot
it = 1;                 % t = 0

%% Loop over deltas
for n = 1:length(deltas)
    delta = deltas(n);

    % Compute total field E at each (R,θ)
    E = zeros(length(R), length(theta));
    for e = 1:az
        rx = (e-1-mid) * d;              % x-position of element e
        ph = delta*(e-1-mid);           % steering phase
        for ix = 1:length(R)
            for iy = 1:length(theta)
                Re = sqrt((x(ix,iy) - rx)^2 + y(ix,iy)^2);
                E(ix,iy) = E(ix,iy) + ...
                           A(e) * cos(omega*t(it) - k*Re + ph);
            end
        end
    end

    % Compute array factor Fa(θ)
    Fa = zeros(1, length(theta));
    for i = 0:az-1
        Fa = Fa + A(i+1) * exp(-1i*i*delta + 1i*k*(i-mid)*d * cos(theta));
    end
    Fa = abs(Fa) / max(abs(Fa));

    %% Plotting for this delta
    figure('Color','w','Position',[100 100 900 400]);
    Fs = 12;

    % Array Factor
    subplot(1,2,1);
    h = polarplot(theta, Fa, 'LineWidth', 2, 'Color', teal);
    title(sprintf('%d-El AF, δ=%.0f° (tapered)', az, delta*180/pi), 'FontSize', Fs);

    % Near-Field Spatial Pattern
    subplot(1,2,2);
    pcolor(x / max(x(:)), y / max(y(:)), E);
    shading interp; axis equal off;
    title(sprintf('Near-Field E, δ=%.0f°', delta*180/pi), 'FontSize', Fs);
    xlabel(sprintf('%d elems, d=λ/2', az), 'FontSize', Fs);
end

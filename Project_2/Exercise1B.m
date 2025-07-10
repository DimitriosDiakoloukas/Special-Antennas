clear all; clc;
%% Simulation parameters
freq = 1e9; % Hz
c = 3e8; % free space speed
lambda = c/freq;
T = 1/freq;
omega = 2*pi*freq;
k = 2*pi/lambda;
Ns = 30; % Number of samples per wavelength
ds = lambda/Ns; % Spatial Discretization
Nt = 35; % Number of time samples per period
dt = T/Nt; % Temporal discretization
t = 0:dt:(1*T); % Increase the number of periods here for longer simulations
R = (0*lambda):ds:(8*lambda);
Ntheta = 240; % Number of angular discretization
dtheta = 2*pi/Ntheta;
theta = 0:dtheta:(2*pi);


%% Generate Domain
x=R.'*cos(theta);
y=R.'*sin(theta);
%% Output properties
teal = [ 0 0.5 0.5]; % maps for unconventional coloring
origBrownColor=[114/256 70/256 43/256];
 
%% Question B - Analysis of linear arrays with 5 and 9 elements at two phase shift values -->%% Animate   
array_sizes = [5, 9];      % number of elements in each scenario
deltaAll   = [0, pi/3, pi/2, 2*pi/3];    % 0°, 60°, 90° and 120° steering
d          = lambda/2;     % spacing = λ/2

for n = 1:length(array_sizes)
    az  = array_sizes(n);       % 5 or 9
    mid = (az-1)/2;             % center‐index for symmetric layout

    for ps = 1:length(deltaAll)
        delta = deltaAll(ps);   % current steering angle
        it    = 1;              % one time snapshot (t=0)

        % build element positions automatically
        x_elem = (-mid:mid)*d;  
        y_elem = zeros(1,az);

        % compute total field E
        E = zeros(length(R),length(theta));
        for e = 1:az
          rx = x_elem(e);
          ry = 0;
          for ix = 1:length(R)
            for iy = 1:length(theta)
              R_e = sqrt((x(ix,iy)-rx)^2 + (y(ix,iy)-ry)^2);
              E(ix,iy) = E(ix,iy) + ...
                cos(omega*t(it) - k*R_e + delta*(e-1-mid));
            end
          end
        end

        % === plot ===
        figure; set(gcf,'Color',[1 1 1]); Fs=10;

        % 1) Array factor (polar)
        subplot(1,2,1); set(gca,'FontSize',Fs);
        A  = ones(1,az);
        Fa=zeros(1,length(theta));
        for i=0:(length(A)-1)
            temp = (A(i+1) * exp(-1i*i*delta + 1i*k*(i*d-((length(A)-1)/2)*d)*cos(theta)));
            Fa = Fa + temp;
        end
        Fa=abs(Fa);
        polar(theta, -Fa/max(Fa)); axis off;

        % 2) Spatial field distribution
        subplot(1,2,2);
        pcolor(x/max(x(:)), y/max(y(:)), E); shading interp;
        pbaspect([1 1 1]); ylim([0 1])

        % labels
        xlabel(sprintf('Array: %d elements, d=\\lambda/2',az),...
               'FontSize',Fs+2);
        title(sprintf('\\delta = %d°', round(delta*180/pi)),...
              'FontSize',Fs+2,'Color',teal);
    end
end

%% PART A
clearvars; close all; clc;

c = 3e8;                     % m/s (included for completeness, though optical frequencies use wavelength directly)
lambda = 1550e-9;            % m
k = 2*pi/lambda;             % wave number
N = 8;                       % number of elements in the array
d = lambda/2;                % spacing between elements

% Steering angles (in degrees)
steer_angles = [30, 60, 90];

% Discretize observation angles for the array factor
theta = linspace(0,2*pi,720);  % 0–360° in 0.5° steps

%% Phase calculation and array factor (AF)
figure('Color','w','Position',[100 100 900 300]);
for idx = 1:length(steer_angles)
    theta0 = steer_angles(idx) * pi/180;  % convert to radians
    
    % Phase for each element: phi_n = -n * k * d * sin(theta0)
    n = 0:(N-1);
    phi_n = -n * k * d * sin(theta0);     % in radians
    
    % Print phase values
    fprintf('Steering = %2d°:\n', steer_angles(idx));
    fprintf(' n    phi_n [rad]    phi_n [deg]\n');
    for m = 1:N
        fprintf('%2d   %10.4f    %10.2f\n', n(m), phi_n(m), phi_n(m)*180/pi);
    end
    fprintf('\n');
    
    % Compute AF(theta)
    % vector version: n = 0:N-1
    % an = 1
    AF = exp(1j*( (0:N-1)'*k*d * (sin(theta) - sin(theta0)) ) );
    AF = sum(AF,1);
    AF_mag = abs(AF)/N;
   
    % Polar plot of the normalized AF
    subplot(1,3,idx);
    polarplot(theta, AF_mag, 'LineWidth',1.8);
    title(sprintf('Steering %d°', steer_angles(idx)));
    thetalim([0 360]);
    rlim([0 1]);
    grid on;
end


%% PART B: Διωνυμική κατανομή βαρών
% Για εφαρμογή Binomial weighting:
% 1. Επιλέγουμε N=8, υπολογίζουμε διωνυμικούς συντελεστές B(n) = nchoosek(N-1, n-1).
% 2. Κανονικοποιούμε ώστε max|B|=1: b_n = B(n)/max(B).
% 3. Αντικαθιστούμε A(n)=b_n στο AF, αφήνουμε φ_n όπως στο μέρος A για την κατεύθυνση.
% Πρακτική υλοποίηση:
% - Ρυθμιστές πλάτους (MZM couplers) σε κάθε κλάδο πριν το phase shifter.
% - Ρυθμίζουν την ένταση κατά αναλογία b_n μέσω ηλεκτρικού σήματος ελέγχου.


clc; clear; close all;

%% ------------------------------------------------------------------------
% Task 1 – just plotting airfoils (yours was fine)
% -------------------------------------------------------------------------

N = 1000;

[naca0018x, naca0018y, ~, ~] = airfoilgen(0,0,18, N);         % NACA 0018
[naca2418x, naca2418y, xcamber2418, ycamber2418] = airfoilgen(2,4,18, N); % NACA 2418
[naca4412x, naca4412y, xcamber4412, ycamber4412] = airfoilgen(4,4,12,N); % NACA 4412

tiledlayout(1,2)

nexttile; % Plot airfoils
plot(naca0018x,naca0018y,'linewidth',2); 
hold on
grid on
yline(0, 'b--', 'LineWidth', 1);
title('NACA 0018')
axis equal
legend('Airfoil', 'Camberline')
xlabel('Normalized Chord')
ylabel('Normalized Thickness')

nexttile
plot(naca2418x,naca2418y,'linewidth',2); 
hold on;
plot(xcamber2418, ycamber2418, 'b--', 'LineWidth', 1);
grid on
title('NACA 2418')
axis equal
legend('Airfoil', 'Camberline')
xlabel('Normalized Chord')
ylabel('Normalized Thickness')

% print('NACA0018_vs_NACA2418', '-dpng', '-r300')

%% Task 2 Deliverable 1

N = linspace(0, 1, 500);

CLtot = zeros(size(N)); % preallocate for the for loop
Numtot = zeros(size(N));

for i = 5:length(N) % start the loop at 5 panels since the method doesn't work with 1

[naca0012x, naca0012y] = airfoilgen(0,0,12, i);

ALPHA = 5; % set AOA to 5 deg

[CL] = Vortex_Panel(naca0012x,naca0012y,ALPHA);

Num = numel(naca0012x) - 1; % total number of panels without repeated last point

CLtot(i) = CL;
Numtot(i) = Num;

end

% exact value from highest panel count
CL_exact = CLtot(end);

% relative error
rel_err = abs(CLtot - CL_exact) ./ CL_exact;

% find first panel count meeting 1% error
idx = find(rel_err < 0.01, 1, 'first');
panels_1percent = Numtot(idx);

figure;
plot(Numtot, CLtot, 'LineWidth', 2)
ylabel('C_l values')
xlabel('Number of Panels')
title('C_l Convergence Study')
ylim([0,1])
print('ConvStudy', '-dpng', '-r300')

%% Task 2 Deliverable 2
N = 1000;

imgPath   = 'NACA0018.png';
[naca0018x , naca0018y] = digitize2D(imgPath, [0 1], 0.18, N);

[naca0006x, naca0006y] = airfoilgen(0,0,06, N);
[naca0012x, naca0012y] = airfoilgen(0,0,12, N);

alphas = -5:20;
cl = zeros(length(alphas),3);

for i = alphas
    cl(i+6,1) = Vortex_Panel(naca0006x,naca0006y,i);
    cl(i+6,2) = Vortex_Panel(naca0012x,naca0012y,i);
    cl(i+6,3) = Vortex_Panel(naca0018x,naca0018y,i);
end

figure()
plot(alphas,cl(:,1))
ylim([-1 3])
title('Cl vs. Angle of Attack for Varied Thickness')
xlabel('Angle of Attack (deg)')
ylabel('c_l')
hold on;
plot(alphas,cl(:,2))
hold on;
plot(alphas,cl(:,3))
legend('NACA 0006', 'NACA 0012', 'NACA 0018');
print('naca_thickness','-dpng')

%% zero-lift angle and lift slope
data = zeros(2,3);

for i = 1:3
    % estimate zero-lift angle and slope from a linear fit
    p = polyfit(alphas,cl(:,i),1);
    y = @(x) p(1)*x + p(2);
    data(1,i) = fzero(y,0);
    data(2,i) = p(1);
end

%% ------------------------------------------------------------------------
% Task 3 – EFFECT OF CAMBER on lift
% NACA 0012 (sym), 2412 (moderate), 4412 (high)
% -------------------------------------------------------------------------
N = 1000;

[naca0012x, naca0012y] = airfoilgen(0,0,12, N);
[naca2418x, naca2418y, xcamber2418, ycamber2418] = airfoilgen(2,4,18, N); % NACA 2418
[naca4412x, naca4412y, xcamber4412, ycamber4412] = airfoilgen(4,4,12,N); % NACA 4412

alphas = -5:20;
nA = numel(alphas);
cl_cam = zeros(nA,3);

for i = 1:nA
    a = alphas(i);
    cl_cam(i,1) = Vortex_Panel(naca0012x, naca0012y, a);
    cl_cam(i,2) = Vortex_Panel(naca2418x, naca2418y, a);
    cl_cam(i,3) = Vortex_Panel(naca4412x, naca4412y, a);
end

% plot all three on the SAME figure (assignment wants this)
figure()
plot(alphas, cl_cam(:,1),'o-','LineWidth',1.2); 
hold on;
plot(alphas, cl_cam(:,2),'s-','LineWidth',1.2);
plot(alphas, cl_cam(:,3),'^-','LineWidth',1.2);
grid on
xlabel('\alpha (deg)'); ylabel('c_l');
title('Effect of Camber on c_l vs \alpha');
legend('NACA 0012','NACA 2412','NACA 4412','Location','NorthWest');

% find zero-lift and slope for the three cambered sections
data_cam = zeros(2,3);  % row1 = alpha_L0, row2 = a0
for j = 1:3
    % only fit near alpha = 0 to get cleaner slope
    idxfit = (alphas >= -2) & (alphas <= 4);
    p = polyfit(alphas(idxfit), cl_cam(idxfit,j), 1);
    f = @(x) p(1)*x + p(2);
    data_cam(1,j) = fzero(f, 0);    % zero-lift AoA
    data_cam(2,j) = p(1);           % slope per deg
end

% show a small table in the command window
fprintf('\nTask 3: Cambered Airfoils – Vortex Panel Results\n');
fprintf('Airfoil\t\talpha_L0 (deg)\ta0 (1/deg)\n');
names = {'NACA 0012','NACA 2412','NACA 4412'};
for j = 1:3
    fprintf('%s\t%10.3f\t%10.4f\n', names{j}, data_cam(1,j), data_cam(2,j));
end

%% Thin Airfoil Theory (TAT) for the same three airfoils

data_cam_thin = zeros(2,3);   % row1 = alpha_L0 (deg), row2 = a0 (1/deg)

% NACA 0012: symmetric, so alpha_L0 = 0, a0 ≈ 2π per rad
a0_rad  = 2*pi;                      % per rad
a0_deg  = a0_rad * (pi/180);         % per deg
data_cam_thin(1,1) = 0.0;            % alpha_L0 (deg)
data_cam_thin(2,1) = a0_deg;         % slope (1/deg)


% NACA 2412: m=2, p=4
[alphaL0_2412_deg, a0_2412_deg] = thin_airfoil_NACA4(2,4);
data_cam_thin(1,2) = alphaL0_2412_deg;
data_cam_thin(2,2) = a0_2412_deg;

% NACA 4412: m=4, p=4
[alphaL0_4412_deg, a0_4412_deg] = thin_airfoil_NACA4(4,4);
data_cam_thin(1,3) = alphaL0_4412_deg;
data_cam_thin(2,3) = a0_4412_deg;

% print in the same format as the vortex panel results
fprintf('\nTask 3: Thin Airfoil Theory Results\n');
fprintf('Airfoil\t\talpha_L0 (deg)\ta0 (1/deg)\n');
for j = 1:3
    fprintf('%s\t%10.3f\t%10.4f\n', names{j}, data_cam_thin(1,j), data_cam_thin(2,j));
end

function [alphaL0_deg, a0_per_deg] = thin_airfoil_NACA4(m_0, p_0)
% Thin airfoil helper
% Thin airfoil theory for 4-digit NACA: NACA mp12
% Inputs:
%   m_0 = first digit (max camber, percent of chord, e.g. 2 for NACA 2412)
%   p_0 = second digit (location of max camber, tenths of chord, e.g. 4 for NACA 2412)
% Outputs:
%   alphaL0_deg : zero-lift angle of attack in degrees
%   a0_per_deg  : lift-curve slope in 1/deg (thin airfoil → ~2π per rad)

m = m_0 / 100;   % max camber as fraction of chord
p = p_0 / 10;    % location of max camber as fraction of chord

% --- Symmetric case: NACA 00xx ---
if m == 0
    alphaL0_deg = 0.0;
    a0_rad      = 2*pi;           % 2π per rad
    a0_per_deg  = a0_rad * (pi/180);  % *** correct: per degree ***
    return;
end

% --- Cambered case: numerical integral in theta-space ---
Nint  = 1000;
theta = linspace(0, pi, Nint);
x     = 0.5*(1 - cos(theta));   % x/c

dyc_dx = zeros(size(x));
for k = 1:Nint
    if x(k) < p
        dyc_dx(k) = (2*m/(p^2))*(p - x(k));
    else
        dyc_dx(k) = (2*m/((1-p)^2))*(p - x(k));
    end
end

% alpha_L0 (rad) = -(1/π) ∫_0^π (dyc/dx)(cosθ - 1) dθ
integrand   = dyc_dx .* (cos(theta) - 1);
alphaL0_rad = -(1/pi) * trapz(theta, integrand);

% Lift-curve slope is still 2π per rad
a0_rad      = 2*pi;
a0_per_deg  = a0_rad * (pi/180);    
alphaL0_deg = alphaL0_rad * 180/pi;
end




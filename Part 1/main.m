clc;clear;close all;

%% Plots for Task 1, Deliverable 1

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
panels_1percent = Numtot(idx)

figure;
plot(Numtot, CLtot, 'LineWidth', 2)
ylabel('C_l values')
xlabel('Number of Panels')
title('C_l Convergence Study')
ylim([0,1])
print('ConvStudy', '-dpng', '-r300')

%% Task 2 Deliverable 2

N = 1000;
[naca0006x, naca0006y] = airfoilgen(0,0,06, N);
[naca0012x, naca0012y] = airfoilgen(0,0,12, N);
[naca0018x, naca0018y] = airfoilgen(0,0,18, N);


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

% alphas = -5:30;
% nA = numel(alphas);
% cl_cam = zeros(nA,3);
% 
% for i = 1:nA
%     a = alphas(i);
%     cl_cam(i,1) = Vortex_Panel(naca0018x, naca0018y, a);
%     cl_cam(i,2) = Vortex_Panel(naca2418x, naca2418y, a);
%     cl_cam(i,3) = Vortex_Panel(naca4412x, naca4412y, a);
% end
% 
% % plot all three on the SAME figure (assignment wants this)
% figure;
% plot(alphas, cl_cam(:,1),'o-','LineWidth',1.2); hold on;
% plot(alphas, cl_cam(:,2),'s-','LineWidth',1.2);
% plot(alphas, cl_cam(:,3),'^-','LineWidth',1.2);
% grid on
% xlabel('\alpha (deg)'); ylabel('c_l');
% title('Effect of Camber on c_l vs \alpha');
% legend('NACA 0012','NACA 2412','NACA 4412','Location','NorthWest');
% 
% % find zero-lift and slope for the three cambered sections
% data_cam = zeros(2,3);  % row1 = alpha_L0, row2 = a0
% for j = 1:3
%     % only fit near alpha = 0 to get cleaner slope
%     idxfit = (alphas >= -2) & (alphas <= 4);
%     p = polyfit(alphas(idxfit), cl_cam(idxfit,j), 1);
%     f = @(x) p(1)*x + p(2);
%     data_cam(1,j) = fzero(f, 0);    % zero-lift AoA
%     data_cam(2,j) = p(1);           % slope per deg
% end
% 
% % show a small table in the command window
% fprintf('\nTask 3: Cambered Airfoils – Vortex Panel Results\n');
% fprintf('Airfoil\t\talpha_L0 (deg)\ta0 (1/deg)\n');
% names = {'NACA 0012','NACA 2412','NACA 4412'};
% for j = 1:3
%     fprintf('%s\t%10.3f\t%10.4f\n', names{j}, data_cam(1,j), data_cam(2,j));
% end
% 

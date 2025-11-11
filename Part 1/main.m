clc;clear;close all;

%% Plots for Task 1, Deliverable 1

N = 1000;

[naca0018x, naca0018y, ~, ~] = airfoilgen(0,0,18, N); % generate airfoils to plot
[naca2418x, naca2418y, xcamber2418, ycamber2418] = airfoilgen(2,4,18, N);


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

% N = 1000;
% [naca0006x, naca0006y] = airfoilgen(0,0,06, N);
% [naca0012x, naca0012y] = airfoilgen(0,0,12, N);
% [naca0018x, naca0018y] = airfoilgen(0,0,18, N);
% 
% 
% alphas = -5:20;
% cl = zeros(length(alphas),3);
% for i = alphas
%     cl(i+6,1) = Vortex_Panel(naca0006x,naca0006y,i);
%     cl(i+6,2) = Vortex_Panel(naca0012x,naca0012y,i);
%     cl(i+6,3) = Vortex_Panel(naca0018x,naca0018y,i);
% end
% %{
% figure()
% plot(alphas,cl(:,1))
% ylim([-1 3])
% title('NACA 0006')
% xlabel('Angle of Attack (deg)')
% ylabel('c_l')
% print('naca0006','-dpng')
% %}
% figure()
% plot(alphas,cl(:,2))
% ylim([-1 3])
% title('NACA 0012')
% xlabel('Angle of Attack (deg)')
% ylabel('c_l')
% print('naca0012','-dpng')
% 
% figure()
% plot(alphas,cl(:,3))
% ylim([-1 3])
% title('NACA 0018')
% xlabel('Angle of Attack (deg)')
% ylabel('c_l')
% print('naca0018','-dpng')
% 
% %% zero-lift angle and lift slope
% data = zeros(2,3);
% 
% for i = 1:3
%     % estimate zero-lift angle and slope from a linear fit
%     p = polyfit(alphas,cl(:,i),1);
%     y = @(x) p(1)*x + p(2);
%     data(1,i) = fzero(y,0);
%     data(2,i) = p(1);
% end
% 
% 
% %% task 3
% [naca0012x, naca0012y] = airfoilgen(0,0,12, N);
% [naca2412x, naca2412y] = airfoilgen(2,4,12, N);
% [naca4412x, naca4412y] = airfoilgen(4,4,12, N);
% 
% alphas = -5:20;
% cl2 = zeros(length(alphas),3);
% for i = alphas
%     cl2(i+6,1) = Vortex_Panel(naca0012x,naca0012y,i);
%     cl2(i+6,2) = Vortex_Panel(naca2412x,naca2412y,i);
%     cl2(i+6,3) = Vortex_Panel(naca4412x,naca4412y,i);
% end
% 
% %{
% figure()
% plot(alphas,cl2(:,1))
% ylim([-1 3])
% title('NACA 0012')
% xlabel('Angle of Attack (deg)')
% ylabel('c_l')
% %}
% 
% 
% figure()
% plot(alphas,cl2(:,2))
% ylim([-1 3])
% title('NACA 2412')
% xlabel('Angle of Attack (deg)')
% ylabel('c_l')
% print('naca2412','-dpng')
% 
% figure()
% plot(alphas,cl2(:,3))
% ylim([-1 3])
% title('NACA 4412')
% xlabel('Angle of Attack (deg)')
% ylabel('c_l')
% print('naca4412','-dpng')
% 
% %% zero-lift angle and lift slope
% data2 = zeros(2,3);
% 
% for i = 1:3
%     % estimate zero-lift angle and slope from a linear fit
%     p = polyfit(alphas,cl2(:,i),1);
%     y = @(x) p(1)*x + p(2);
%     data2(1,i) = fzero(y,0);
%     data2(2,i) = p(1);
% end

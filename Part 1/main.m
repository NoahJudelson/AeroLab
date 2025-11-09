clc;clear;close all;

%% Plots for Task 1, Deliverable 1

[naca0018x, naca0018y, ~, ~] = airfoilgen(0,0,18); % generate airfoils to plot
[naca2418x, naca2418y, xcamber2418, ycamber2418] = airfoilgen(2,4,18);


tiledlayout(1,2)

nexttile; % Plot airfoils
plot(naca0018x,naca0018y,'linewidth',2); 
grid on
title('NACA 0018')
ylim([-0.5,0.5]);
legend('Airfoil', 'Camberline')

nexttile
plot(naca2418x,naca2418y,'linewidth',2); 
hold on;
plot(xcamber2418, ycamber2418, 'b--', 'LineWidth', 1);
grid on
title('NACA 2418')
ylim([-0.5,0.5]);
legend('Airfoil', 'Camberline')

print('NACA 0018 vs. NACA 2418', '-dpng', '-r300')

%% Task 2 Deliverable 2
[naca0006x, naca0006y] = airfoilgen(0,0,06);
[naca0012x, naca0012y] = airfoilgen(0,0,12);
[naca0018x, naca0018y] = airfoilgen(0,0,18);


alphas = -5:20;
cl = zeros(length(alphas),3);
for i = alphas
    cl(i+5,1) = Vortex_Panel(naca0006x,naca0006y,i);
    cl(i+5,2) = Vortex_Panel(naca0012x,naca0012y,i);
    cl(i+5,3) = Vortex_Panel(naca0018x,naca0018y,i);
end

figure()
plot(alphas,cl(:,1))


figure()
plot(alphas,cl(:,2))

figure()
plot(alphas,cl(:,3))

%% zero-lift angle and lift slope
data = zeros(2,3);

for i = 1:3
    % estimate zero-lift angle and slope from a linear fit
    p = polyfit(alphas,cl(:,i),1);
    y = @(x) p(1)*x + p(2);
    data(1,i) = fzero(y,0);
    data(2,i) = p(1);
end


%% task 3
[naca0012x, naca0012y] = airfoilgen(0,0,12);
[naca02412x, naca02412y] = airfoilgen(2,4,12);
[naca4412x, naca4412y] = airfoilgen(4,4,12);

alphas = -5:20;
cl2 = zeros(length(alphas),3);
for i = alphas
    cl2(i+5,1) = Vortex_Panel(naca0012x,naca0012y,i);
    cl2(i+5,2) = Vortex_Panel(naca2412x,naca2412y,i);
    cl2(i+5,3) = Vortex_Panel(naca4412x,naca4412y,i);
end

figure()
plot(alphas,cl2(:,1))


figure()
plot(alphas,cl2(:,2))

figure()
plot(alphas,cl2(:,3))

%% zero-lift angle and lift slope
data2 = zeros(2,3);

for i = 1:3
    % estimate zero-lift angle and slope from a linear fit
    p = polyfit(alphas,cl2(:,i),1);
    y = @(x) p(1)*x + p(2);
    data2(1,i) = fzero(y,0);
    data2(2,i) = p(1);
end
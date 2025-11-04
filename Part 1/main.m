clc;clear;close all;
%% task 2 deliverable 2
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

%{
data = zeros(4,3);
for i = 1:3
    negativeAlpha = find((cl(:,i)<0),1,"first");
din


end
%}



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


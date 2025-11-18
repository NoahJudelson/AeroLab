clc;clear;close all;

a0_t = 2*pi;
a0_r = 2*pi;
c_t = linspace(0,1);
c_r = ones(size(c_t));
aero_t = 0; 
aero_r = 0;
geo_t = 5*pi/180;
geo_r = 5*pi/180;
N = 50;
e = ones(size(c_t));
figure()
for i = 4:2:10

    for j = 1:length(c_t)
        b = i*(c_r(j)+c_t(j))/2;
        [e(j),~,~] = PLLT(b,a0_t,a0_r,c_t(j),1,aero_t,aero_r,geo_t,geo_r,N);
    end

    delta = (1./e) - 1;
    hold on
    plot(c_t,delta);
    hold on
    xlabel('Taper Ratio (c_t/c_r)')
    ylabel('Induced Drag Factor (\delta)')

end

legend('AR = 4','AR = 6','AR = 8','AR = 10')
print('Part 2 plot','-dpng')
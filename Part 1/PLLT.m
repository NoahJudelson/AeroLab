function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%{  
t = tips
r = root

b = span
a0 = lift slope (per rad)
c = chord
aero = zero-lift angle of attack (rad)
geo = geometric angle of attack (rad)

N = number of terms in series
%}

theta = zeros([N 1]);
p = zeros([N 1]);
coeff_A = zeros([N N]);
for i = 1:N
    % define thetas and y location
    theta(i) = (i*pi)/(2*N);
    y = (-b/2)*cos(theta(i));
    % find a_0, zero-lift angle at each theta to get RHS
    aero_y = (2*(aero_r-aero_t)/b)*y + aero_r;
    alpha_y = (2*(geo_r-geo_t)/b)*y + geo_r;
    p(i) = alpha_y-aero_y;

    % find a0, c at each theta to get coefficient matrix
    a0_y = (2*(a0_r-a0_t)/b)*y + a0_r;
    c_y = (2*(c_r-c_t)/b)*y + c_r;
    for j = 1:N
        n = 2*j-1;
        coeff_A(i,j) = (4*b/(a0_y*c_y))*sin(n*theta(i)) + n*sin(n*theta(i))/sin(theta(i));
    end
end
% solve for coefficients
A = coeff_A\p;

%calculate CL
area = ((c_r+c_t)/2)*b;
AR = b*b/area;
c_L = A(1)*pi*AR;

%calculate CD_i
delta = 0;
for i = 2:N
    n = 2*i-1;
    delta = delta + n*((A(i)/A(1))^2);
end

c_Di = ((c_L*c_L)/(pi*AR)) * (1+delta);
e = 1/(1+delta);

end


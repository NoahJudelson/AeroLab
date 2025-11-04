function [x, y] = airfoilgen(m_0, p_0, t_0)

% This function generates a plot of the location of points on a 4 digit
% NACA airfoil, inputting all geometric values from the 4 digit code.

m_use = m_0 / 100; % calculate useful values of m,p,t
p_use = p_0 / 10;
t_use = t_0 / 100;

c = t_use / (t_0 / 100); % calculate chord length


x = linspace(0,c,1000); % get x vector


% y_t equation
y_t = (t_use / 0.2) * c * ((0.2969*sqrt(x./c)) - (0.1260 * (x ./ c)) ...
    - (0.3516 * (x ./ c).^2) + (0.2843 * (x ./ c).^3) - (0.1036 * (x ./ c).^4));


% mean camber line function
if x < (p_use * c)
    y_c1 = m_use * (x ./ (p_use^2)) * ((2 * p_use) - (x ./ c));

else
    y_c2 = m_use * ((c-x) ./ ((1-p_use)^2)) .* (1 + (x/c) - (2 * p_use));

y_c = [y_c1 y_c2];


xi = ;


x_u = x - (y_t .* sin(xi));
x_l = x + (y_t .* sin(xi));
y_u = y_c + (y_t .* cos(xi));
y_l = y_c - (y_t .* cos(xi));







end


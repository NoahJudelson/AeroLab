clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3 Task 1: Effect of Angle of Attack on C_L (Wing) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_panels = 400;
alpha_range_deg = -5:1:10;
alpha_range_rad = deg2rad(alpha_range_deg);

%% Root: NACA 2412
[x_root, y_root, ~, ~] = airfoilgen(2,4,12, N_panels);
CL_root = zeros(size(alpha_range_deg));
for k = 1:numel(alpha_range_deg)
    CL_root(k) = Vortex_Panel(x_root, y_root, alpha_range_deg(k));
end
p_root     = polyfit(alpha_range_rad, CL_root, 1);
a0_r       = p_root(1);                 % lift-curve slope root (1/rad)
alphaL0_r  = -p_root(2)/p_root(1);      % zero-lift AoA root (rad)

%% Tip: NACA 0012
[x_tip, y_tip, ~, ~] = airfoilgen(0,0,12, N_panels);
CL_tip = zeros(size(alpha_range_deg));
for k = 1:numel(alpha_range_deg)
    CL_tip(k) = Vortex_Panel(x_tip, y_tip, alpha_range_deg(k));
end
p_tip      = polyfit(alpha_range_rad, CL_tip, 1);
a0_t       = p_tip(1);                  % lift-curve slope tip (1/rad)
alphaL0_t  = -p_tip(2)/p_tip(1);        % zero-lift AoA tip (rad)

%% Wing geometry and PLLT
b   = 36;                               % span [ft]
c_r = 5 + 4/12;                         % root chord [ft]
c_t = 3 + 7/12;                         % tip chord [ft]
twist_root_minus_tip = 2;              % root is 2 deg higher AoA than tip
N_terms = 50;                           % PLLT series terms

alpha_tip_deg = -5:0.5:15;              % geometric AoA at the tip (deg)
CL_wing = zeros(size(alpha_tip_deg));

for k = 1:numel(alpha_tip_deg)
    geo_t = deg2rad(alpha_tip_deg(k));
    geo_r = geo_t + deg2rad(twist_root_minus_tip);
    [~, CL_wing(k), ~] = PLLT(b, a0_t, a0_r, c_t, c_r, ...
                              alphaL0_t, alphaL0_r, ...
                              geo_t, geo_r, N_terms);
end

figure;
plot(alpha_tip_deg, CL_wing, 'LineWidth', 1.5);
grid on; xlabel('\alpha_{tip} (deg)'); ylabel('C_L');
title('Cessna 180: C_L vs Tip Angle of Attack');
print('part3task1','-dpng');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3 Task 2: Estimation of Profile Drag c_D0 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model sectional drag coefficient c_d for the NACA 0012 tip airfoil
% as a function of angle of attack, using Theory of Wing Sections data.

%% 2.1 Experimental drag-polar data for NACA 0012 (TIP)
% Replace with data extracted from the NACA 0012 drag polar figure.

CL_exp = [-0.4  -0.2  0.0  0.2  0.4  0.6  0.8  1.0  1.2];
CD_exp = [0.0065  0.006  0.006  0.0061  0.007  0.0075  0.008  0.01 0.012];

CL_exp = CL_exp(:);    % ensure column
CD_exp = CD_exp(:);

%% 2.2 Fit quadratic drag polar: c_d = c_d0 + k * C_l^2

X = [ones(size(CL_exp)), CL_exp.^2];   % design matrix [1, C_l^2]
beta = X \ CD_exp;                     % least-squares fit

cd0 = beta(1);                         % c_d at zero lift
k    = beta(2);                        % quadratic coefficient

fprintf('Tip drag polar fit: c_d = %.5f + %.5f * C_l^2\n', cd0, k);

%% 2.3 Build c_d(alpha) MODEL using the lift curve from Task 1

alpha_model_deg = -5:0.5:15;           % angle of attack range
alpha_model_rad = deg2rad(alpha_model_deg);

% a0_t and alphaL0_t were computed in Task 1
CL_model = a0_t * (alpha_model_rad - alphaL0_t);   % predicted sectional C_l(α)
CD_model = cd0 + k * CL_model.^2;                  % quadratic c_d(α) model

%% 2.4 Convert EXPERIMENTAL C_l to EXPERIMENTAL α using same lift curve

alpha_exp_rad = CL_exp ./ a0_t + alphaL0_t;
alpha_exp_deg = rad2deg(alpha_exp_rad);

%% 2.5 Plot c_d vs alpha (MODEL + EXPERIMENTAL)

figure;
hold on;
plot(alpha_model_deg, CD_model, 'LineWidth',1.6, 'DisplayName','Model c_d(\alpha)');
plot(alpha_exp_deg, CD_exp, 'o', 'MarkerSize',6, ...
     'DisplayName','Experimental c_d(\alpha)');
grid on; box on;
xlabel('\alpha_{tip} (deg)');
ylabel('c_d');
title('NACA 0012 Tip Section: Profile Drag Coefficient vs Angle of Attack');
legend('Location','best');
print('part3task2','-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3 Task 3: Effect of Angle of Attack on C_D (Wing) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Angle of attack range (same as Task 1)
alpha_deg = alpha_tip_deg;   % uses -5 : 0.5 : 15 from Task 1
alpha_rad = deg2rad(alpha_deg);

% Preallocate
CDi_wing = zeros(size(alpha_deg));
CD0_wing = zeros(size(alpha_deg));
CD_total = zeros(size(alpha_deg));

% Wing planform area (ft^2 → convert to m^2 only if needed later)
S = (c_r + c_t)/2 * b;

for n = 1:numel(alpha_deg)

    % ---- 1. Induced drag coefficient from PLLT ----
    geo_t = deg2rad(alpha_deg(n));
    geo_r = geo_t + deg2rad(twist_root_minus_tip);

    [~, CL_here, CDi_here] = PLLT(b, a0_t, a0_r, c_t, c_r, ...
                                  alphaL0_t, alphaL0_r, ...
                                  geo_t, geo_r, N_terms);

    CDi_wing(n) = CDi_here;

    % ---- 2. Profile drag coefficient using tip-section model ----
    % From Task 2: c_d = cd0 + k * c_l^2
    CL_tip_local = a0_t * (geo_t - alphaL0_t);   % sectional CL at tip
    CD0_wing(n) = cd0 + k * CL_tip_local.^2;

    % ---- 3. Total drag coefficient ----
    CD_total(n) = CDi_wing(n) + CD0_wing(n);
end

% ---- Plot for Task 3 ----
figure;
plot(alpha_deg, CD_total, 'k-', 'LineWidth',1.5, 'DisplayName','Total C_D'); hold on;
plot(alpha_deg, CDi_wing, 'r--', 'LineWidth',1.5, 'DisplayName','Induced C_{D,i}');
plot(alpha_deg, CD0_wing, 'b-.', 'LineWidth',1.5, 'DisplayName','Profile C_{D,0}');
grid on; xlabel('\alpha (deg)'); ylabel('C_D');
title('Cessna 180: Drag Coefficients vs Angle of Attack');
legend('Location','best');
print('part3task3','-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 3 Task 4: Effect of Airspeed on Thrust Required (Wing) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = 2500;                     % weight (lb)
rho = 0.001756;               % slug/ft^3 at 10,000 ft (standard)
V_knots = 40:1:140;           % airspeed range
V_ft_s = V_knots * 1.68781;   % convert to ft/s

T_required = zeros(size(V_ft_s));

for i = 1:numel(V_ft_s)
    V = V_ft_s(i);

    % ---- 1. Find CL required for steady level flight ----
    CL_req = 2*W / (rho * V^2 * S);

    % ---- 2. Convert CL_req to geometric AoA at the tip ----
    alpha_req_rad = CL_req / a0_t + alphaL0_t;
    alpha_req_deg = rad2deg(alpha_req_rad);

    % Bound angle to model limits
    alpha_req_deg = max(min(alpha_req_deg, max(alpha_tip_deg)), min(alpha_tip_deg));

    % ---- 3. Interpolate drag coefficients from Task 3 ----
    CD_now = interp1(alpha_deg, CD_total, alpha_req_deg, 'linear', 'extrap');

    % ---- 4. Compute thrust required (= drag) ----
    T_required(i) = 0.5 * rho * V^2 * S * CD_now;   % pounds
end

%% Plot required thrust vs airspeed
figure;
plot(V_knots, T_required, 'LineWidth',1.5);
grid on; xlabel('Airspeed (knots)'); ylabel('Thrust Required (lb)');
title('Cessna 180: Thrust Required for Steady, Level Flight');
print('part3task4','-dpng');


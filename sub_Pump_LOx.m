% sub_Pump_LOx.m
function [properties, key_values] = sub_Pump_LOx(inputs, properties)

%% 1. Thermodynamics & Hydraulics
p_in = properties.p;
T_in = properties.T;

properties.warnings = struct('code', {}, 'message', {}, 'value', {});

rho_in = double(py.CoolProp.CoolProp.PropsSI('D', 'T', T_in, 'P', p_in, 'Oxygen'));
mu_in  = double(py.CoolProp.CoolProp.PropsSI('VISCOSITY', 'T', T_in, 'P', p_in, 'Oxygen'));
h_in   = double(py.CoolProp.CoolProp.PropsSI('H', 'T', T_in, 'P', p_in, 'Oxygen'));
s_in   = double(py.CoolProp.CoolProp.PropsSI('S', 'T', T_in, 'P', p_in, 'Oxygen'));

p_intermediate = p_in + inputs.delta_p_pump_LOx;
p_final = p_intermediate - inputs.delta_p_partial;

h_out_s = double(py.CoolProp.CoolProp.PropsSI('H', 'P', p_intermediate, 'S', s_in, 'Oxygen'));
h_out_intermediate = h_in + (h_out_s - h_in) / inputs.eta_pump_LOx;

P_pump = inputs.m_dot_oxidizer * (h_out_intermediate - h_in);

p_vapor = double(py.CoolProp.CoolProp.PropsSI('P', 'T', T_in, 'Q', 0, 'Oxygen'));

Head_Total = inputs.delta_p_pump_LOx / (rho_in * inputs.g0);
Q = inputs.m_dot_oxidizer / rho_in;

% Turbine team forced RPM
N_max_rpm = 26500;
% With this RPM and a lower total head for LOx, 1 stage is sufficient
n_stages = 1;
Head_per_stage = Head_Total / n_stages;

%% 2. Design inputs

use_inducer = true; % Required due to high 20k RPM

eta_inducer = 0.60; 
delta_p_inducer = 1.0e5;   % [Pa] (1 bar boost from the inducer)
alpha1_inducer = 0.0;      % [rad]
w1u_factor_inducer = 0.50; 

% -------- Impeller design inputs (SET B inspired) --------
eta_hyd_guess = 0.80;      % Match realistic correlation
alpha1_base = 0.0;
delta_h = 0.30;
b2_over_D2 = 0.09;

N_ss = 410;

%% 3. Effective inlet conditions after inducer

p_in_eff = p_in;
alpha1_eff = alpha1_base;

if use_inducer
    p_in_eff = p_in + delta_p_inducer;
    alpha1_eff = alpha1_inducer;
end

D_pipe_suction = 0.08; % [m] Suction pipe diameter (LOx requires slightly larger than LCH4)
A_pipe_suction = pi/4 * D_pipe_suction^2; % [m^2] Suction pipe area
c_abs_suction = Q / A_pipe_suction; % [m/s] Absolute velocity at suction pipe

NPSH_system = (p_in - p_vapor) / (rho_in * inputs.g0) + c_abs_suction^2 / (2 * inputs.g0);

rpm_limit = N_ss * (NPSH_system^0.75) / sqrt(Q);
Security_factor = 0.9;
N_max_rpm = Security_factor * rpm_limit;

Omega = N_max_rpm * (pi / 30);

%% 4. Specific speed and preliminary coefficients

N_s_universal = (Omega * sqrt(Q)) / (inputs.g0 * Head_per_stage)^0.75;

% Hydraulic efficiency correlation used as estimate
m = 0.08 * (1 / Q)^0.15 * (0.85 / N_s_universal)^0.06;
eta_hyd_corr = 1 - 0.065 * (1 / Q)^m - 0.23 * (0.3 - log10(2.3 * N_s_universal))^2 * (1 / Q)^0.05;
    
% Assume standard aerospace losses to bridge the gap
eta_vol  = 0.91;  % 9% lost to internal leakage
eta_mech = 0.90;  % 10% lost to bearing/seal friction
    
% Calculate Total Pump Efficiency
eta_total_calculated = eta_hyd_corr * eta_vol * eta_mech; 

% Use guessed efficiency for geometry closure, keep correlated one as reference
eta_hyd = eta_hyd_guess;

% SET B-style guided coefficients
psi_id = 0.605 * exp(-0.408 * N_s_universal);
psi = psi_id * eta_hyd;

lambda_c = -0.103 * log10(N_s_universal) + 1.1855;
lambda_w = 0.2144 * exp(0.1745 * N_s_universal);

%% 5. Impeller sizing

u2 = sqrt(inputs.g0 * Head_per_stage / psi_id);
D2 = 60 * u2 / (pi * N_max_rpm);

% Rotor tip diameter ratio guided by cavitation correlation
delta_t = sqrt(delta_h^2 + ...
    0.588 * (N_s_universal^(4/3)) * psi_id * (((lambda_c + lambda_w) / lambda_w)^(1/3)));

D1h = delta_h * D2;
D1t = delta_t * D2;
b1 = (D1t - D1h) / 2;

if b1 <= 0
    error('b1 <= 0. Revisa delta_t, delta_h o la correlación de entrada.');
end

% Inlet meridional velocity from geometry
c1m_base = 2 * Q / (pi * D2 * (delta_t + delta_h) * b1);
c1m_eff = c1m_base;

phi_base = c1m_base / u2;
phi_eff = c1m_eff / u2;

% Rotor meridional velocity ratio from geometry
xi = (delta_t^2 - delta_h^2) / (4 * b2_over_D2);

if xi <= 0
    error('xi <= 0. Revisa delta_t, delta_h y b2_over_D2.');
end

% Degree of reaction becomes an output
R = 1 ...
    - psi / 2 ...
    + (phi_eff^2 / (2 * psi)) * ((1 - xi^2) + tan(alpha1_eff)^2 * (1 - delta_t^2)) ...
    - phi_eff * delta_t * tan(alpha1_eff);

% Flow angles
alpha2 = atan( psi / (xi * phi_eff) + (delta_t / xi) * tan(alpha1_eff) );
beta1  = atan( delta_t / phi_eff - tan(alpha1_eff) );
beta2  = atan( (1 / (phi_eff * xi)) * (1 - psi) - (delta_t / xi) * tan(alpha1_eff) );

% Velocities
c1m = c1m_eff;
c1u = c1m * tan(alpha1_eff);

c2m = u2 * xi * phi_eff;
c2u = c2m * tan(alpha2);

A1 = Q / c1m;
A2 = Q / c2m;

b2 = b2_over_D2 * D2;

w1u = c1m * tan(beta1);
w2u = c2m * tan(beta2);

w1 = c1m / cos(beta1);
w2 = c2m / cos(beta2);
c1 = c1m / cos(alpha1_eff);
c2 = c2m / cos(alpha2);

%% 6. Cavitation check & Geometry Validation

g = inputs.g0;

% 1. Required input parameters for NPSH calculations
sigma_inducer = 0.065;
sigma_impeller = 0.2;
B_inducer = 0.15;
B_impeller = 0.12;
delta_p_between_inducer_impeller = 0.01 * delta_p_inducer; % Assume ~1% of total head loss

% 2. Variables reutilizadas
delta_p_inducer_static = delta_p_inducer; 
D_ind_tip_in = D1t;                       
D_ind_hub_in = D1h;                       
c_u_inducer_in = 0;                       
p_suction = p_in;                         

%% 6.1 Inducer NPSH & Geometry

A_inducer_in_geom = pi/4 * (D_ind_tip_in^2 - D_ind_hub_in^2);
A_inducer_in_eff = A_inducer_in_geom * (1 - B_inducer);

if A_inducer_in_eff <= 0
    error('The effective inducer inlet area is not positive. Check blockage factor.');
end

c_m_inducer_LE = Q / A_inducer_in_eff;

% NPSH disponible en la entrada del inductor
NPSH_A_inducer = (p_suction - p_vapor) / (rho_in * g) + c_abs_suction^2 / (2*g);

% Triángulo de velocidades y NPSH requerido del inducer
u_inducer_tip = Omega * D_ind_tip_in / 2;
u_ind_hub     = Omega * D_ind_hub_in / 2;

w_inducer_tip_LE = hypot(c_m_inducer_LE, u_inducer_tip - c_u_inducer_in);

NPSH_R_inducer = sigma_inducer * w_inducer_tip_LE^2 / (2*g);
margin_inducer = NPSH_A_inducer - NPSH_R_inducer;

%% 6.1b Inducer exit swirl 

Delta_h_ind = delta_p_inducer_static / (rho_in * g);
c_u_ind_out_tip = Delta_h_ind * g / u_inducer_tip;
c_u_ind_out_hub = Delta_h_ind * g / u_ind_hub;
c_u_ind_out_mean = 0.5 * (c_u_ind_out_tip + c_u_ind_out_hub);

%% 6.2 Main Impeller NPSH & Geometry

p_impeller_eye = p_suction + delta_p_inducer_static - delta_p_between_inducer_impeller;

D_impeller_tip_in = D1t;
D_impeller_hub_in = D1h;

A_impeller_in_geom = pi/4 * (D_impeller_tip_in^2 - D_impeller_hub_in^2);
A_impeller_in_eff = A_impeller_in_geom * (1 - B_impeller);

if A_impeller_in_eff <= 0
    error('The effective impeller inlet area is not positive.');
end

c_m_impeller_LE = Q / A_impeller_in_eff;
u_impeller_tip_in = Omega * D_impeller_tip_in / 2;

c_u_impeller_in = c_u_ind_out_mean;

c_abs_impeller_eye = hypot(c_m_impeller_LE, c_u_impeller_in);
w_impeller_LE = hypot(c_m_impeller_LE, u_impeller_tip_in - c_u_impeller_in);

NPSH_A_impeller = (p_impeller_eye - p_vapor) / (rho_in * g) + c_abs_impeller_eye^2 / (2*g);
NPSH_R_impeller = sigma_impeller * w_impeller_LE^2 / (2*g);
margin_impeller = NPSH_A_impeller - NPSH_R_impeller;

%% 6.3 Inducer flow angles for CAD export

beta_ind_out_flujo = atan2(u_inducer_tip - c_u_ind_out_tip, c_m_inducer_LE);
beta_ind_hub_flujo = atan2(u_ind_hub - c_u_ind_out_hub, c_m_inducer_LE);

i_inc = 3 * (pi/180);
beta_pala_tip = beta_ind_out_flujo + i_inc;
beta_pala_hub = beta_ind_hub_flujo + i_inc;

%% 6.4 Slip Factor & Ángulo Físico de la Pala (Impulsor)

N_B_R = 6; 
delta_M = (delta_t + delta_h) / 2;

beta2B_guess = beta2; 
error_beta = 1;
iter = 0;

while error_beta > 1e-5 && iter < 100
    iter = iter + 1;
    SF_wiesner = 1 - (sqrt(cos(beta2B_guess)) / (N_B_R^0.7));
    delta_M_lim = exp(-8.16 * cos(beta2B_guess) / N_B_R);
    
    if delta_M > delta_M_lim
        SF = SF_wiesner * (1 - ((delta_M - delta_M_lim) / (1 - delta_M_lim))^3);
    else
        SF = SF_wiesner;
    end
    
    c2u_inf = c2u / SF;
    beta2B_new = atan((u2 - c2u_inf) / c2m);
    error_beta = abs(beta2B_new - beta2B_guess);
    beta2B_guess = beta2B_new;
end

beta2B = beta2B_guess;
properties.beta2B = rad2deg(beta2B);
properties.Slip_Factor = SF;
properties.N_blades = N_B_R;

% Guardar propiedades del Inductor para CAD
properties.Inducer_D_tip = D_ind_tip_in;
properties.Inducer_D_hub = D_ind_hub_in;
properties.Inducer_c_m = c_m_inducer_LE;
properties.Inducer_beta_pala_tip_deg = rad2deg(beta_pala_tip);
properties.Inducer_beta_pala_hub_deg = rad2deg(beta_pala_hub);
properties.Inducer_N_blades = 3;

%% 7. Warnings

if delta_t <= delta_h
    msg = 'Geometría no válida: delta_t <= delta_h, lo que implica b1 <= 0.';
    warning(msg);
    properties.warnings(end+1).code = 'INVALID_B1';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = delta_t - delta_h;
end

if margin_inducer < 0
    msg = sprintf('¡Peligro! Margen de cavitación del Inductor negativo (%.4f m)', margin_inducer);
    warning(msg);
    properties.warnings(end+1).code = 'INDUCER_CAVITATION';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = margin_inducer;
end

if margin_impeller < 0
    msg = sprintf('¡Peligro! Margen de cavitación del Impulsor negativo (%.4f m).', margin_impeller);
    warning(msg);
    properties.warnings(end+1).code = 'IMPELLER_CAVITATION';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = margin_impeller;
end

if b2_over_D2 < 0.04 || b2_over_D2 > 0.20
    msg = sprintf('b2/D2 = %.4f está fuera del rango preliminar recomendado (~0.04-0.20).', b2_over_D2);
    warning(msg);
    properties.warnings(end+1).code = 'B2_RATIO_WARNING';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = b2_over_D2;
end

if b2 < 0.5e-3
    msg = sprintf('b2 = %.6f m es demasiado pequeño desde el punto de vista de fabricación.', b2);
    warning(msg);
    properties.warnings(end+1).code = 'B2_TOO_SMALL';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = b2;
end

if phi_eff < 0.1 || phi_eff > 0.25
    msg = sprintf('phi = %.4f está fuera del rango preliminar típico (~0.10-0.25).', phi_eff);
    warning(msg);
    properties.warnings(end+1).code = 'PHI_OUT_OF_RANGE';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = phi_eff;
end

if psi < 0.4 || psi > 0.6
    msg = sprintf('psi = %.4f está fuera del rango preliminar típico (~0.40-0.60).', psi);
    warning(msg);
    properties.warnings(end+1).code = 'PSI_OUT_OF_RANGE';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = psi;
end

if R < 0.5 || R > 0.75
    msg = sprintf('R = %.4f está fuera del rango preliminar típico (~0.50-0.75).', R);
    warning(msg);
    properties.warnings(end+1).code = 'R_OUT_OF_RANGE';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = R;
end

%% 8. Outputs

properties.p = p_final;
properties.T = py.CoolProp.CoolProp.PropsSI('T', 'P', properties.p, 'H', h_out_intermediate, 'Oxygen');
properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Oxygen');
properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Oxygen');

properties.NPSH_system = NPSH_system;
properties.NPSH_A_inducer = NPSH_A_inducer;
properties.NPSH_R_inducer = NPSH_R_inducer;
properties.NPSH_A_impeller = NPSH_A_impeller;
properties.NPSH_R_impeller = NPSH_R_impeller;

properties.Margin_Inducer = margin_inducer;
properties.Margin_Impeller = margin_impeller;
properties.Head_Total = Head_Total;
properties.Head_Per_Stage = Head_per_stage;
properties.Volumetric_Q = Q;
properties.RPM = N_max_rpm;
properties.Specific_Speed_Ns_Universal = N_s_universal;
properties.P_Pump_LOx = P_pump;

properties.n_stages = n_stages;
properties.u2 = u2;
properties.D2 = D2;
properties.D1h = D1h;
properties.D1t = D1t;
properties.b1 = b1;
properties.b2 = b2;
properties.delta_h = delta_h;
properties.delta_t = delta_t;
properties.xi = xi;
properties.R = R;

properties.A1 = A1;
properties.A2 = A2;

properties.c1m = c1m;
properties.c1u = c1u;
properties.c2m = c2m;
properties.c2u = c2u;
properties.w1 = w1;
properties.w2 = w2;
properties.w1u = w1u;
properties.w2u = w2u;
properties.c1 = c1;
properties.c2 = c2;

properties.b2_over_D2 = b2_over_D2;

properties.eta_hyd = eta_hyd;
properties.eta_total_calculated = eta_total_calculated;
properties.eta_hyd_corr = eta_hyd_corr;
properties.psi_id = psi_id;
properties.psi = psi;
properties.mu_in = mu_in;

properties.alpha1_base = alpha1_base;
properties.alpha1_eff = alpha1_eff;
properties.phi_base = phi_base;
properties.phi_eff = phi_eff;
properties.c1m_base = c1m_base;
properties.c1m_eff = c1m_eff;
properties.p_in = p_in;
properties.p_in_eff = p_in_eff;
properties.delta_p_inducer = delta_p_inducer;

properties.alpha1 = alpha1_eff;
properties.beta1 = beta1;
properties.alpha2 = alpha2;

properties.warning_messages = {properties.warnings.message};

key_values = struct();
key_values.P_pump           = P_pump;
key_values.N_rpm            = N_max_rpm;
key_values.D_impeller       = D2;
key_values.u2               = u2;
key_values.psi              = psi;
key_values.phi              = phi_eff;
key_values.n_stages         = n_stages;
key_values.Head             = Head_Total; 
key_values.Margin_Inducer   = margin_inducer;
key_values.Margin_Impeller  = margin_impeller;
key_values.NPSH_system      = NPSH_system;
key_values.N_max_cavitation = rpm_limit;
key_values.Ns_pump          = N_s_universal;

end
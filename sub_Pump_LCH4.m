% sub_Pump_LCH4.m
function [properties, key_values] = sub_Pump_LCH4(inputs, properties)

%% 1. Thermodynamics & Hydraulics
p_in = properties.p;
T_in = properties.T;

rho_in = double(py.CoolProp.CoolProp.PropsSI('D', 'T', T_in, 'P', p_in, 'Methane'));
mu_in  = double(py.CoolProp.CoolProp.PropsSI('VISCOSITY', 'T', T_in, 'P', p_in, 'Methane'));
h_in   = double(py.CoolProp.CoolProp.PropsSI('H', 'T', T_in, 'P', p_in, 'Methane'));
s_in   = double(py.CoolProp.CoolProp.PropsSI('S', 'T', T_in, 'P', p_in, 'Methane'));

p_intermediate = p_in + inputs.delta_p_pump_LCH4;
p_final = p_intermediate - inputs.delta_p_partial;

h_out_s = double(py.CoolProp.CoolProp.PropsSI('H', 'P', p_intermediate, 'S', s_in, 'Methane'));
h_out_intermediate = h_in + (h_out_s - h_in) / inputs.eta_pump_LCH4;

P_pump = inputs.m_dot_fuel * (h_out_intermediate - h_in);

p_vapor = double(py.CoolProp.CoolProp.PropsSI('P', 'T', T_in, 'Q', 0, 'Methane'));

Head_Total = inputs.delta_p_pump_LCH4 / (rho_in * inputs.g0);
Q = inputs.m_dot_fuel / rho_in;

n_stages = 4;
Head_per_stage = Head_Total / n_stages;

%% 2. Design inputs

use_inducer = true;

eta_inducer = 0.60; %#ok
delta_p_inducer = 1.2e5;   % [Pa]
alpha1_inducer = 0.0;      % [rad]
w1u_factor_inducer = 0.50; % only used for NPSH_R correction if desired

% -------- Impeller design inputs (SET B inspired) --------
eta_hyd_guess = 0.7;
alpha1_base = 0.0;
delta_h = 0.30;

% New design input: outlet width ratio
b2_over_D2 = 0.09;

% Suction specific speed parameter
N_ss = 280;

%% 3. Effective inlet conditions after inducer

p_in_eff = p_in;
alpha1_eff = alpha1_base;

if use_inducer
    p_in_eff = p_in + delta_p_inducer;
    alpha1_eff = alpha1_inducer;
end

NPSH = (p_in_eff - p_vapor) / (rho_in * inputs.g0);
N_max_rpm = N_ss * (NPSH^0.75) / sqrt(Q);
Omega = N_max_rpm * (pi / 30);

%% 4. Specific speed and preliminary coefficients

N_s_universal = (Omega * sqrt(Q)) / (inputs.g0 * Head_per_stage)^0.75;

% Hydraulic efficiency correlation used as estimate
m = 0.08 * (1 / Q)^0.15 * (0.85 / N_s_universal)^0.06;
eta_hyd_corr = 1 - 0.065 * (1 / Q)^m - 0.23 * (0.3 - log10(2.3 * N_s_universal))^2 * (1 / Q)^0.05;

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

% Keep inducer influence on pressure/NPSH, not on c1m unless model is upgraded
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

%% 6. Cavitation check

w1u_eff = w1u;
if use_inducer
    w1u_eff = w1u * w1u_factor_inducer;
end

NPSH_R = lambda_c * (c1m^2 / (2 * inputs.g0)) + lambda_w * (w1u_eff^2 / (2 * inputs.g0));

%% Intento calcular n_ss

dn = D1h;
d1 = D1t;
kn = 1 - (dn/d1)^2;

denom_nss = ((lambda_c + lambda_w)*phi_eff^2 + lambda_w*(1 - phi_eff/tan(alpha1_eff))^2)^0.75;

nss_formula = 158*sqrt(phi_eff*kn)/denom_nss;



%% 7. Warnings

properties.warnings = struct('code', {}, 'message', {}, 'value', {});

if delta_t <= delta_h
    msg = 'Geometría no válida: delta_t <= delta_h, lo que implica b1 <= 0.';
    warning(msg);
    properties.warnings(end+1).code = 'INVALID_B1';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = delta_t - delta_h;
end

if b2_over_D2 < 0.04
    msg = sprintf('b2/D2 = %.4f está por debajo del rango preliminar recomendado (~0.04-0.20).', b2_over_D2);
    warning(msg);
    properties.warnings(end+1).code = 'B2_RATIO_LOW';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = b2_over_D2;
elseif b2_over_D2 > 0.20
    msg = sprintf('b2/D2 = %.4f está por encima del rango preliminar recomendado (~0.04-0.20).', b2_over_D2);
    warning(msg);
    properties.warnings(end+1).code = 'B2_RATIO_HIGH';
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

if ~isreal(alpha2) || ~isreal(beta1) || ~isreal(beta2)
    msg = 'Alguno de los ángulos de velocidad no es real. Revisa cierre cinemático.';
    warning(msg);
    properties.warnings(end+1).code = 'INVALID_ANGLES';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = NaN;
end

if NPSH_R > NPSH
    msg = sprintf('NPSH_R = %.4f m es mayor que NPSH = %.4f m. Margen de cavitación insuficiente.', NPSH_R, NPSH);
    warning(msg);
    properties.warnings(end+1).code = 'NPSH_INSUFFICIENT';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = NPSH_R - NPSH;
end

if abs(eta_hyd_corr - eta_hyd_guess) > 0.05
    msg = sprintf('La eficiencia asumida (%.4f) difiere bastante de la correlación preliminar (%.4f).', eta_hyd_guess, eta_hyd_corr);
    warning(msg);
    properties.warnings(end+1).code = 'ETA_HYD_MISMATCH';
    properties.warnings(end).message = msg;
    properties.warnings(end).value = eta_hyd_corr - eta_hyd_guess;
end

%% 8. Outputs

properties.p = p_final;
properties.T = py.CoolProp.CoolProp.PropsSI('T', 'P', properties.p, 'H', h_out_intermediate, 'Methane');
properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Methane');
properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Methane');

properties.NPSH = NPSH;
properties.NPSH_R = NPSH_R;
properties.Head_Total = Head_Total;
properties.Head_Per_Stage = Head_per_stage;
properties.Volumetric_Q = Q;
properties.RPM = N_max_rpm;
properties.Specific_Speed_Ns_Universal = N_s_universal;

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
properties.eta_hyd_corr = eta_hyd_corr;
properties.psi_id = psi_id;
properties.psi = psi;
properties.mu_in = mu_in;

% Store both base and effective inlet data
properties.alpha1_base = alpha1_base;
properties.alpha1_eff = alpha1_eff;
properties.phi_base = phi_base;
properties.phi_eff = phi_eff;
properties.c1m_base = c1m_base;
properties.c1m_eff = c1m_eff;
properties.p_in = p_in;
properties.p_in_eff = p_in_eff;
properties.delta_p_inducer = delta_p_inducer;

properties.nss_formula = nss_formula;

% Angles for printing / postprocess
properties.alpha1 = alpha1_eff;
properties.beta1 = beta1;
properties.alpha2 = alpha2;
properties.beta2 = beta2;
properties.beta2B = "inutil";

% Warnings
properties.warning_messages = {properties.warnings.message};

key_values = P_pump;

end
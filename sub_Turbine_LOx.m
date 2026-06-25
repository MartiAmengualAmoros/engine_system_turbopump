function [properties, key_values] = sub_Turbine_LOx(inputs, properties)
% SUB_TURBINE_LOx  Turbine 2 — drives the LOx pump shaft.
% Working fluid : methane (fuel stream continues from Turbine 1)
% Expansion     : p_mid -> p_turbine_LOx_out (= p_turbine_exit)

% ================================================================
% 0. READ INLET STATE
% ================================================================
p_in  = properties.p;
T_in  = properties.T;
p_out = inputs.p_turbine_LOx_out;

% ================================================================
% 1. THERMODYNAMIC MODEL
% ================================================================
h_in    = py.CoolProp.CoolProp.PropsSI('H', 'T', T_in, 'P', p_in,  'Methane');
s_in    = py.CoolProp.CoolProp.PropsSI('S', 'T', T_in, 'P', p_in,  'Methane');
h_out_s = py.CoolProp.CoolProp.PropsSI('H', 'P', p_out, 'S', s_in, 'Methane');

h_out       = h_in - inputs.eta_turbine_LOx * (h_in - h_out_s);
W_available = inputs.m_dot_fuel * (h_in - h_out);  % [W]
delta_h_is  = h_in - h_out_s;                      % [J/kg]

% ================================================================
% 2. STAGE COUNT
% ================================================================
C0_full   = sqrt(2 * delta_h_is);
u_blade   = pi * inputs.D_mean_LCH4 * inputs.N_shaft_LOx / 60;
nu_single = u_blade / C0_full;
nu_target = inputs.nu_target_LOx;

if nu_single >= 0.40
    n_stages      = 1;
    architecture  = 'Single-stage impulse';
    delta_h_stage = delta_h_is;
elseif nu_single >= 0.22
    n_stages      = 1;
    architecture  = 'Curtis (2-row velocity-compounded)';
    delta_h_stage = delta_h_is;
else
    delta_h_stage = 0.5 * (u_blade / nu_target)^2;
    n_stages      = ceil(delta_h_is / delta_h_stage);
    delta_h_stage = delta_h_is / n_stages;
    architecture  = sprintf('Rateau pressure-compounded (%d stages)', n_stages);
end

C0_stage  = sqrt(2 * delta_h_stage);
nu_actual = u_blade / C0_stage;

% ================================================================
% 3. VELOCITY TRIANGLES
% ================================================================
C1         = 0.96 * C0_stage;
alpha1_rad = deg2rad(18);
alpha1_deg = 18;

C1u  = C1 * cos(alpha1_rad);
C1ax = C1 * sin(alpha1_rad);

W1u = C1u - u_blade;
W1  = sqrt(W1u^2 + C1ax^2);
beta1_deg = rad2deg(atan2(C1ax, W1u));

C2u  = 0;
C2ax = C1ax;
W2u  = C2u - u_blade;
W2   = sqrt(W2u^2 + C2ax^2);
beta2_deg = rad2deg(atan2(C2ax, abs(W2u)));

psi      = u_blade * (C1u - C2u) / u_blade^2;
phi_flow = C1ax / u_blade;
R_deg    = 1 - (C1u + C2u) / (2 * u_blade);

% ================================================================
% 4. BLADE GEOMETRY
% ================================================================
rho_in  = py.CoolProp.CoolProp.PropsSI('D', 'T', T_in, 'P', p_in, 'Methane');
A_flow  = inputs.m_dot_fuel / (rho_in * C1ax);
h_blade = A_flow / (pi * inputs.D_mean_LOx);
D_tip   = inputs.D_mean_LOx + h_blade;
D_hub   = inputs.D_mean_LOx - h_blade;
HTR     = D_hub / D_tip;

% ================================================================
% 5. EXIT STATE
% ================================================================
T_out   = py.CoolProp.CoolProp.PropsSI('T',      'P', p_out, 'H', h_out, 'Methane');
rho_out = py.CoolProp.CoolProp.PropsSI('D',      'P', p_out, 'T', T_out, 'Methane');
cp_out  = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', T_out, 'P', p_out, 'Methane');

properties.T   = T_out;
properties.p   = p_out;
properties.rho = rho_out;
properties.c_p = cp_out;

% ================================================================
% 6. FEASIBILITY
% ================================================================
surplus = W_available - inputs.P_turbine_LOx_needed;

% ================================================================
% 7. KEY VALUES
% ================================================================
key_values.delta_p      = p_in - p_out;
key_values.W_available  = W_available;
key_values.surplus      = surplus;
key_values.delta_h_is   = delta_h_is;
key_values.PR           = p_in / p_out;
key_values.T_out        = T_out;
key_values.n_stages     = n_stages;
key_values.architecture = architecture;
key_values.delta_h_stage = delta_h_stage;
key_values.C0_full      = C0_full;
key_values.C0_stage     = C0_stage;
key_values.u_blade      = u_blade;
key_values.nu_actual    = nu_actual;
key_values.C1           = C1;
key_values.C1u          = C1u;
key_values.C1ax         = C1ax;
key_values.W1           = W1;
key_values.W2           = W2;
key_values.alpha1_deg   = alpha1_deg;
key_values.beta1_deg    = beta1_deg;
key_values.beta2_deg    = beta2_deg;
key_values.psi          = psi;
key_values.phi_flow     = phi_flow;
key_values.R_deg        = R_deg;
key_values.h_blade      = h_blade;
key_values.D_mean       = inputs.D_mean_LOx;
key_values.D_tip        = D_tip;
key_values.D_hub        = D_hub;
key_values.HTR          = HTR;
key_values.A_flow       = A_flow;

end

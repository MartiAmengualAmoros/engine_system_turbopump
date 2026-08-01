function [properties, key_values] = sub_Turbine_LOx(inputs, properties)
% SUB_TURBINE_LOx  Multi-stage turbine — LOx shaft (working fluid: methane).
% Stage count auto-computed to satisfy Ns >= Ns_min AND D_mean <= D_mean_max_LOx
% Efficiency mode controlled by inputs.eta_turbine_mode:
%   'fixed'  -> uses inputs.eta_turbine_LOx directly
%   'balje'  -> computes eta from Balje Ns correlation with design margin

p_in = properties.p;
T_in = properties.T;
p_out = inputs.p_turbine_LOx_out;
p_out_safe = max(min(p_out, 1000e5), 0.1e5);

%% 1. TOTAL THERMODYNAMICS
h_in   = py.CoolProp.CoolProp.PropsSI('H','T',T_in,'P',p_in,'Methane');
s_in   = py.CoolProp.CoolProp.PropsSI('S','T',T_in,'P',p_in,'Methane');
try
    h_out_s = py.CoolProp.CoolProp.PropsSI('H','P',p_out_safe,'S',s_in,'Methane');
catch
    cp_in = py.CoolProp.CoolProp.PropsSI('CPMASS','T',T_in,'P',p_in,'Methane');
    cv_in = py.CoolProp.CoolProp.PropsSI('CVMASS','T',T_in,'P',p_in,'Methane');
    gamma = cp_in/cv_in;
    T_out_s_est = max(T_in*(p_out_safe/p_in)^((gamma-1)/gamma), 93.3);
    h_out_s = py.CoolProp.CoolProp.PropsSI('H','T',T_out_s_est,'P',p_out_safe,'Methane');
end

delta_h_is_total = h_in - h_out_s;

if delta_h_is_total <= 0
    key_values = zero_kv_LOx(inputs, T_in);
    properties.T   = T_in; properties.p = p_out;
    properties.rho = py.CoolProp.CoolProp.PropsSI('D','T',T_in,'P',p_out_safe,'Methane');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS','T',T_in,'P',p_out_safe,'Methane');
    return
end

%% 2. EXIT DENSITY (isentropic — used for Ns and HTR)
try
    rho_out_est = py.CoolProp.CoolProp.PropsSI('D','P',p_out_safe,'S',s_in,'Methane');
catch
    cp_in = py.CoolProp.CoolProp.PropsSI('CPMASS','T',T_in,'P',p_in,'Methane');
    cv_in = py.CoolProp.CoolProp.PropsSI('CVMASS','T',T_in,'P',p_in,'Methane');
    gamma = cp_in/cv_in;
    T_out_s_est = max(T_in*(p_out_safe/p_in)^((gamma-1)/gamma), 93.3);
    rho_out_est = py.CoolProp.CoolProp.PropsSI('D','T',T_out_s_est,'P',p_out_safe,'Methane');
end
Q_vol = inputs.m_dot_fuel / rho_out_est;

%% 3. AUTO STAGE COUNT
omega         = inputs.N_shaft_LOx * pi / 30;
C0_full       = sqrt(2 * delta_h_is_total);
Ns_1stage     = (inputs.N_shaft_LOx/60) * sqrt(Q_vol) / delta_h_is_total^0.75;
D_mean_1stage = 2 * inputs.nu_target_LOx * C0_full / omega;

Ns_min  = 0.055;
n_Ns    = ceil((Ns_min / max(Ns_1stage, 1e-12))^(1/0.75));
n_D     = ceil((D_mean_1stage / inputs.D_mean_max_LOx)^2);
n_stages = max([n_Ns, n_D, 1]);

%% 4. PER-STAGE GEOMETRY
delta_h_stage = delta_h_is_total / n_stages;
C0_stage = sqrt(2 * delta_h_stage);
u1       = inputs.nu_target_LOx * C0_stage;
D_rotor  = 2 * u1 / omega;
r1       = D_rotor / 2;

alpha1_rad = deg2rad(70);
W_euler    = delta_h_stage;
C1u = W_euler / u1;
C1r = C1u / tan(alpha1_rad);
C1  = sqrt(C1u^2 + C1r^2);
W1u = C1u - u1;
W1  = sqrt(W1u^2 + C1r^2);
beta1_deg = rad2deg(atan2(C1r, abs(W1u)));

r2        = 0.60 * r1;
D_exducer = 2 * r2;
u2        = omega * r2;
nu_actual = u1 / C0_stage;

%% 5. HTR from exit axial velocity constraint (C2ax/u1 = 0.25)
% Note: dense supercritical methane at inter-turbine pressure often makes
% HTR^2 < 0 (required area > full disc). Fallback HTR=0.40 is physically justified.
C2ax_target  = 0.25 * u1;
A_ex_target  = inputs.m_dot_fuel / (rho_out_est * C2ax_target);
K            = pi/2 * D_exducer^2;
HTR_sq       = (K - A_ex_target) / (K + A_ex_target);
if HTR_sq <= 0 || HTR_sq >= 1
    HTR_ex = 0.40;
else
    HTR_ex = min(max(sqrt(HTR_sq), 0.20), 0.75);
end

D_shroud  = D_exducer * sqrt(2/(1+HTR_ex^2));
D_hub_ex  = HTR_ex * D_shroud;
A_ex      = pi/4 * (D_shroud^2 - D_hub_ex^2);
C2ax      = inputs.m_dot_fuel / (rho_out_est * A_ex);
rho_in    = py.CoolProp.CoolProp.PropsSI('D','T',T_in,'P',p_in,'Methane');
b1        = inputs.m_dot_fuel / (rho_in * C1r * pi * D_rotor);

%% 6. Ns PER STAGE
Ns_stage   = (inputs.N_shaft_LOx/60) * sqrt(Q_vol) / delta_h_stage^0.75;
eta_Ns_max = min(max(0.87 - 6.0*(Ns_stage - 0.18)^2, 0.50), 0.90);

%% 6b. EFFICIENCY — fixed or Balje-based
if isfield(inputs,'eta_turbine_mode') && strcmp(inputs.eta_turbine_mode,'balje')
    eta_margin  = 0.92;
    eta_turbine = eta_margin * eta_Ns_max;
    fprintf('[T-LOx]  Balje mode: Ns=%.4f  eta_Ns_max=%.3f  eta_used=%.3f\n', ...
        Ns_stage, eta_Ns_max, eta_turbine);
else
    eta_turbine = inputs.eta_turbine_LOx;
    if eta_turbine > eta_Ns_max + 0.03
        warning('sub_Turbine_LOx: eta_fixed=%.3f may be optimistic for Ns=%.4f (Balje max~%.3f).', ...
            eta_turbine, Ns_stage, eta_Ns_max);
    end
end

%% 6c. TURBINE TYPE SELECTION (Balje diagram)
%% 6c. TURBINE TYPE SELECTION (Balje diagram + Stage Override)
if n_stages >= 3
    turbine_type = 'Axial impulse turbine (Pressure-compounded)';
elseif Ns_stage < 0.05
    turbine_type = 'Radial inflow turbine (RIT) - partial admission';
elseif Ns_stage <= 0.40
    turbine_type = 'Radial inflow turbine (RIT)';
elseif Ns_stage <= 1.0
    turbine_type = 'Mixed-flow / Axial impulse turbine';
else
    turbine_type = 'Axial reaction turbine';
end

%% 7. THERMODYNAMICS with resolved eta
h_out       = h_in - eta_turbine * delta_h_is_total;
W_available = inputs.m_dot_fuel * (h_in - h_out);
W_euler     = eta_turbine * delta_h_stage;
C1u         = W_euler / u1;
psi_stage   = W_euler / u1^2;
eta_euler   = min(max((u1*C1u) / delta_h_stage, 0), 1);

%% 8. EXIT STATE
T_out   = py.CoolProp.CoolProp.PropsSI('T','P',p_out_safe,'H',h_out,'Methane');
rho_out = py.CoolProp.CoolProp.PropsSI('D','P',p_out_safe,'T',T_out,'Methane');
cp_out  = py.CoolProp.CoolProp.PropsSI('CPMASS','T',T_out,'P',p_out_safe,'Methane');
properties.T   = T_out; properties.p = p_out;
properties.rho = rho_out; properties.c_p = cp_out;

surplus = W_available - inputs.P_turbine_LOx_needed;

%% 9. KEY VALUES
key_values.W_available  = W_available;
key_values.surplus      = surplus;
key_values.delta_p      = p_in - p_out;
key_values.delta_h_is   = delta_h_is_total;
key_values.delta_h_stage = delta_h_stage;
key_values.PR           = p_in / p_out;
key_values.T_out        = T_out;
key_values.architecture = turbine_type;
key_values.n_stages     = n_stages;
key_values.nu_actual    = nu_actual;
key_values.C0           = C0_stage;
key_values.C0_full      = C0_full;
key_values.C0_stage     = C0_stage;
key_values.u1           = u1;   key_values.u_blade = u1;
key_values.C1           = C1;   key_values.C1u = C1u;
key_values.C1r          = C1r;  key_values.C1ax = C1r;
key_values.W1           = W1;
key_values.alpha1_deg   = 70;
key_values.beta1_deg    = beta1_deg;
key_values.u2           = u2;   key_values.W2u = u2;
key_values.W2           = u2;   key_values.beta2_deg = 0;
key_values.C2ax         = C2ax; key_values.C2ax_over_u1 = C2ax/u1;
key_values.psi          = psi_stage;
key_values.phi_flow     = C1r/u1;
key_values.R_deg        = 0.5;
key_values.D_rotor      = D_rotor;  key_values.D_mean = D_rotor;
key_values.D_exducer    = D_exducer;
key_values.D_shroud     = D_shroud; key_values.D_tip = D_shroud;
key_values.D_hub_exit   = D_hub_ex; key_values.D_hub = D_hub_ex;
key_values.b1           = b1;       key_values.h_blade = b1;
key_values.A_exducer    = A_ex;     key_values.A_flow = A_ex;
key_values.HTR_exducer  = HTR_ex;   key_values.HTR = HTR_ex;
key_values.Ns           = Ns_stage;
key_values.Ns_1stage    = Ns_1stage;
key_values.eta_turbine  = eta_turbine;
key_values.eta_euler    = eta_euler;
key_values.eta_Ns_max   = eta_Ns_max;
key_values.n_blade_min  = 11;
end

function kv = zero_kv_LOx(inputs, T_in)
fields = {'W_available','surplus','delta_p','delta_h_is','delta_h_stage','PR','T_out', ...
    'nu_actual','C0','C0_full','C0_stage','u1','u_blade','C1','C1u','C1r','C1ax', ...
    'W1','alpha1_deg','beta1_deg','u2','W2u','W2','beta2_deg','C2ax','C2ax_over_u1', ...
    'psi','phi_flow','R_deg','D_rotor','D_mean','D_exducer','D_shroud','D_tip', ...
    'D_hub_exit','D_hub','b1','h_blade','A_exducer','A_flow','HTR_exducer','HTR', ...
    'Ns','Ns_1stage','eta_turbine','eta_euler','eta_Ns_max','n_blade_min'};
for i = 1:numel(fields); kv.(fields{i}) = 0; end
kv.surplus      = -inputs.P_turbine_LOx_needed;
kv.T_out        = T_in;
kv.architecture = 'Infeasible';
kv.n_stages     = 0;
end
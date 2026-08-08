function [key_values, properties_flow] = cycle_solver()

% Solves for [delta_p_pump_LCH4, p_mid] using Newton-Raphson.
% Two shaft-balance equations:
%   (1) W_turbine_LCH4 = P_pump_LCH4_needed
%   (2) W_turbine_LOx  = P_pump_LOx_needed

[inputs, pf0, po0] = engine_inputs();

p_mid_lb      = inputs.p_fuel_injector_inlet * 1.05;
p_T1_in_lb    = p_mid_lb * 1.10;
p_pump_out_lb = p_T1_in_lb + inputs.delta_p_cooling_channels + inputs.delta_p_partial;
dp_pump_lb    = p_pump_out_lb - pf0.p;
lb = [dp_pump_lb; p_mid_lb];

x0 = [inputs.delta_p_pump_LCH4; inputs.p_turbine_LCH4_out];
F  = @(x) system_residual(x, inputs, pf0, po0);

ws = warning('off', 'all');
[x_sol, res_norm] = newton2d(F, x0, lb, 1.0, 1e2, 60);
warning(ws);

inputs.delta_p_pump_LCH4  = x_sol(1);
inputs.p_turbine_LCH4_out = x_sol(2);
inputs.p_turbine_LOx_out  = inputs.p_fuel_injector_inlet;

r = system_residual(x_sol, inputs, pf0, po0);
fprintf('Cycle solver exit: delta_p_pump_LCH4 = %.4f MPa | p_mid = %.4f MPa\n', ...
    x_sol(1)/1e6, x_sol(2)/1e6);
fprintf('  Residuals: LCH4 = %.2f W | LOx = %.2f W | |r| = %.3f W\n', r(1), r(2), res_norm);

if res_norm > 1e3
    error('Cycle infeasible: |r| = %.1f W. Reduce p_CC_req or increase Q_dot.', res_norm);
end
fprintf('Cycle CONVERGED (|r| < 1 kW).\n');

clear key_values properties_flow
try
    [key_values, properties_flow] = run_cycle(inputs, pf0, po0);
catch ME
    fprintf('Final run_cycle FAILED: %s\n', ME.message);
    rethrow(ME);
end
plot_cycle(key_values, properties_flow)

%% Pump report
fprintf('\n--- PUMP REPORT ---\n');
fprintf('[P-LCH4] P_pump          : %.2f kW\n',  key_values.P_Pump_LCH4 / 1e3);
fprintf('[P-LCH4] N_shaft         : %.0f rpm\n', key_values.N_shaft_LCH4);
fprintf('[P-LCH4] D_impeller      : %.1f mm\n',  key_values.D_impeller_LCH4 * 1e3);
%fprintf('[P-LCH4] u2              : %.1f m/s\n', key_values.u2_LCH4);
fprintf('[P-LCH4] psi             : %.3f\n',     key_values.psi_LCH4);
fprintf('[P-LCH4] phi             : %.3f\n',     key_values.phi_LCH4);
fprintf('[P-LCH4] n_stages        : %d\n',       key_values.n_stages_pump_LCH4);
fprintf('[P-LCH4] Head_total      : %.1f m\n',   key_values.Head_total_LCH4);
fprintf('[P-LCH4] NPSH            : %.2f m\n',   key_values.NPSH_LCH4);
fprintf('[P-LCH4] N_max_cavitation: %.0f rpm\n', key_values.N_max_cavitation_LCH4);
fprintf('[P-LCH4] Ns_pump         : %.4f\n',     key_values.Ns_pump_LCH4);
fprintf('[P-LOx]  P_pump          : %.2f kW\n',  key_values.P_Pump_LOx / 1e3);
fprintf('[P-LOx]  N_shaft         : %.0f rpm\n', key_values.N_shaft_LOx);
fprintf('[P-LOx]  D_impeller      : %.1f mm\n',  key_values.D_impeller_LOx * 1e3);
%fprintf('[P-LOx]  u2              : %.1f m/s\n', key_values.u2_LOx);
fprintf('[P-LOx]  psi             : %.3f\n',     key_values.psi_LOx);
fprintf('[P-LOx]  phi             : %.3f\n',     key_values.phi_LOx);
fprintf('[P-LOx]  n_stages        : %d\n',       key_values.n_stages_pump_LOx);
fprintf('[P-LOx]  Head            : %.1f m\n',   key_values.Head_LOx);
fprintf('[P-LOx]  NPSH           : %.2f m\n',   key_values.NPSH_LOx);
fprintf('[P-LOx]  N_max_cavitation: %.0f rpm\n', key_values.N_max_cavitation_LOx);
fprintf('[P-LOx]  Ns_pump         : %.4f\n',     key_values.Ns_pump_LOx);

%% Turbine report
fprintf('\n--- TURBINE REPORT ---\n');
fprintf('[T-LCH4] Architecture : %s\n',      key_values.architecture_LCH4);
fprintf('[T-LCH4] n_stages     : %d\n',      key_values.n_stages_LCH4);
fprintf('[T-LCH4] nu_actual    : %.3f\n',    key_values.nu_actual_LCH4);
fprintf('[T-LCH4] D_mean       : %.1f mm\n', key_values.D_mean_LCH4 * 1e3);
fprintf('[T-LCH4] h_blade      : %.2f mm\n', key_values.h_blade_LCH4 * 1e3);
fprintf('[T-LCH4] HTR          : %.3f\n',    key_values.HTR_LCH4);
fprintf('[T-LCH4] Ns           : %.4f\n',    key_values.Ns_LCH4);
fprintf('[T-LCH4] eta_turbine  : %.3f\n',    key_values.eta_turbine_LCH4);
fprintf('[T-LOx]  Architecture : %s\n',      key_values.architecture_LOx);
fprintf('[T-LOx]  n_stages     : %d\n',      key_values.n_stages_LOx);
fprintf('[T-LOx]  nu_actual    : %.3f\n',    key_values.nu_actual_LOx);
fprintf('[T-LOx]  D_mean       : %.1f mm\n', key_values.D_mean_LOx * 1e3);
fprintf('[T-LOx]  h_blade      : %.2f mm\n', key_values.h_blade_LOx * 1e3);
fprintf('[T-LOx]  HTR          : %.3f\n',    key_values.HTR_LOx);
fprintf('[T-LOx]  Ns           : %.4f\n',    key_values.Ns_LOx);
fprintf('[T-LOx]  eta_turbine  : %.3f\n',    key_values.eta_turbine_LOx);

end

% -------------------------------------------------------------------------

function [x, res_norm] = newton2d(F, x0, lb, ftol, xtol, max_iter)
x        = max(x0(:), lb);
res_norm = inf;
for k = 1:max_iter
    f        = F(x);
    res_norm = norm(f);
    if res_norm < ftol; break; end
    J = zeros(2,2);
    for j = 1:2
        h      = max(abs(x(j)) * 1e-4, 1e4);
        xp     = x; xp(j) = xp(j) + h;
        try; fp = F(xp); catch; fp = f * 2; end
        J(:,j) = (fp - f) / h;
    end
    if rcond(J) < 1e-12
        warning('newton2d: near-singular Jacobian at iter %d — stopping', k);
        break;
    end
    step  = J \ (-f);
    alpha = 1.0;
    for ls = 1:12
        x_try = max(x + alpha * step, lb);
        try; f_try = F(x_try); if norm(f_try) < res_norm; break; end; catch; end
        alpha = alpha * 0.5;
    end
    x = max(x + alpha * step, lb);
    if max(abs(alpha * step)) < xtol; break; end
end
end

% -------------------------------------------------------------------------

function residual = system_residual(x, inputs, pf0, po0)
inputs.delta_p_pump_LCH4  = x(1);
inputs.p_turbine_LCH4_out = x(2);
inputs.p_turbine_LOx_out  = inputs.p_fuel_injector_inlet;
ws = warning('off', 'all');
[~, ~, W_LCH4, W_LOx, P_LCH4, P_LOx] = run_cycle(inputs, pf0, po0);
warning(ws);
residual = [W_LCH4 - P_LCH4; W_LOx - P_LOx];
end

% -------------------------------------------------------------------------

function [key_values, properties_flow, W_LCH4, W_LOx, P_LCH4, P_LOx] = run_cycle(inputs, properties_fuel, properties_oxidizer)

properties_flow.fuel.Tanks     = properties_fuel;
properties_flow.oxidizer.Tanks = properties_oxidizer;

%% LCH4 Pump
[properties_fuel, kv_pump_LCH4] = sub_Pump_LCH4(inputs, properties_fuel);
properties_flow.fuel.Pump_LCH4  = properties_fuel;

key_values.P_Pump_LCH4           = kv_pump_LCH4.P_pump;
key_values.N_shaft_LCH4          = kv_pump_LCH4.N_rpm;
key_values.D_impeller_LCH4       = kv_pump_LCH4.D_impeller;
key_values.psi_LCH4              = kv_pump_LCH4.psi;
key_values.phi_LCH4              = kv_pump_LCH4.phi;
key_values.n_stages_pump_LCH4    = kv_pump_LCH4.n_stages;
key_values.Head_total_LCH4       = kv_pump_LCH4.Head_total;
key_values.Head_per_stage_LCH4   = kv_pump_LCH4.Head_per_stage;
key_values.NPSH_LCH4             = kv_pump_LCH4.NPSH_system;
key_values.N_max_cavitation_LCH4 = kv_pump_LCH4.N_max_cavitation;
key_values.Ns_pump_LCH4          = kv_pump_LCH4.Ns_pump;

inputs.P_turbine_LCH4_needed = kv_pump_LCH4.P_pump;
inputs.N_shaft_LCH4          = kv_pump_LCH4.N_rpm;

%% LOx Pump
[properties_oxidizer, kv_pump_LOx] = sub_Pump_LOx(inputs, properties_oxidizer);
properties_flow.oxidizer.Pump_LOx  = properties_oxidizer;

key_values.P_Pump_LOx            = kv_pump_LOx.P_pump;
key_values.N_shaft_LOx           = kv_pump_LOx.N_rpm;
key_values.D_impeller_LOx        = kv_pump_LOx.D_impeller;
key_values.psi_LOx               = kv_pump_LOx.psi;
key_values.phi_LOx               = kv_pump_LOx.phi;
key_values.n_stages_pump_LOx     = kv_pump_LOx.n_stages;
key_values.Head_LOx              = kv_pump_LOx.Head;
key_values.NPSH_LOx              = kv_pump_LOx.NPSH;
key_values.N_max_cavitation_LOx  = kv_pump_LOx.N_max_cavitation;
key_values.Ns_pump_LOx           = kv_pump_LOx.Ns_pump;

inputs.P_turbine_LOx_needed = kv_pump_LOx.P_pump;
inputs.N_shaft_LOx          = kv_pump_LOx.N_rpm;

%% Cooling Channels
[properties_fuel, kv_cooling]         = sub_Cooling_channels(inputs, properties_fuel);
key_values.delta_T_Cooling_Channels   = kv_cooling;
key_values.Q_dot_Cooling              = inputs.Q_dot;
properties_flow.fuel.Cooling_Channels = properties_fuel;

%% LCH4 Turbine
[properties_fuel, kv_turb_LCH4]   = sub_Turbine_LCH4(inputs, properties_fuel);
properties_flow.fuel.Turbine_LCH4  = properties_fuel;

key_values.delta_p_Turbine_LCH4 = kv_turb_LCH4.delta_p;
key_values.W_Turbine_LCH4       = kv_turb_LCH4.W_available;
key_values.architecture_LCH4    = kv_turb_LCH4.architecture;
key_values.n_stages_LCH4        = kv_turb_LCH4.n_stages;
key_values.nu_actual_LCH4       = kv_turb_LCH4.nu_actual;
key_values.h_blade_LCH4         = kv_turb_LCH4.h_blade;
key_values.HTR_LCH4             = kv_turb_LCH4.HTR;
key_values.Ns_LCH4              = kv_turb_LCH4.Ns;
key_values.D_mean_LCH4          = kv_turb_LCH4.D_mean;
key_values.eta_turbine_LCH4     = kv_turb_LCH4.eta_turbine;

%% LOx Turbine
[properties_fuel, kv_turb_LOx]   = sub_Turbine_LOx(inputs, properties_fuel);
properties_flow.fuel.Turbine_LOx  = properties_fuel;

key_values.delta_p_Turbine_LOx  = kv_turb_LOx.delta_p;
key_values.W_Turbine_LOx        = kv_turb_LOx.W_available;
key_values.architecture_LOx     = kv_turb_LOx.architecture;
key_values.n_stages_LOx         = kv_turb_LOx.n_stages;
key_values.nu_actual_LOx        = kv_turb_LOx.nu_actual;
key_values.h_blade_LOx          = kv_turb_LOx.h_blade;
key_values.HTR_LOx              = kv_turb_LOx.HTR;
key_values.Ns_LOx               = kv_turb_LOx.Ns;
key_values.D_mean_LOx           = kv_turb_LOx.D_mean;
key_values.eta_turbine_LOx      = kv_turb_LOx.eta_turbine;

W_LCH4 = kv_turb_LCH4.W_available;
W_LOx  = kv_turb_LOx.W_available;
P_LCH4 = inputs.P_turbine_LCH4_needed;
P_LOx  = inputs.P_turbine_LOx_needed;

%% Injectors
[properties_fuel,     key_values.TC_Ingoing_fuel]     = sub_Injector_CH4(inputs, properties_fuel);
properties_flow.fuel.Injector = properties_fuel;

[properties_oxidizer, key_values.TC_Ingoing_oxidizer] = sub_Injector_LOx(inputs, properties_oxidizer);
properties_flow.oxidizer.Injector = properties_oxidizer;

%% Thrust Chamber
key_values.Thrust_Chamber = sub_Thrust_Chamber(inputs, properties_oxidizer, properties_fuel);
end

function [key_values, properties_flow] = cycle_solver()

    % Solves for [delta_p_pump_LCH4, p_mid] simultaneously using Newton-Raphson
    % with backtracking line search (no Optimization Toolbox required).
    % Two equations enforced at once:
    %   (1) W_LCH4(delta_p, p_mid) = P_LCH4_needed   — LCH4 shaft balance
    %   (2) W_LOx(delta_p, p_mid)  = P_LOx_needed    — LOx shaft balance

    % TODO: Add *more* visualization tools for the flow values all over the cycle

    [inputs, pf0, po0] = engine_inputs();

    % Initial guess: current pump delta_p + midpoint pressure for p_mid
    p_in_est   = pf0.p + inputs.delta_p_pump_LCH4 - inputs.delta_p_cooling_channels;
   p_mid_init = (inputs.p_fuel_injector_inlet + p_in_est) / 2;
    x0 = [inputs.delta_p_pump_LCH4; p_mid_init];

    % Physical lower bounds: pump must add pressure; p_mid must exceed exit
    lb = [1e6; inputs.p_fuel_injector_inlet + 1e5];

    F = @(x) system_residual(x, inputs, pf0, po0);

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
        error('Cycle infeasible: |r| = %.1f W — turbines cannot produce enough shaft work to drive the pumps. Reduce p_CC_req or increase Q_dot.', res_norm);
    end
    fprintf('Cycle CONVERGED (|r| < 1 kW).\n');

    [key_values, properties_flow] = run_cycle(inputs, pf0, po0);
    plot_cycle(key_values, properties_flow)

    % ---- FINAL TURBINE REPORT (reads from key_values forwarded by run_cycle) ----
fprintf('\n===== TURBINE RESULTS =====\n');
fprintf('[T-LCH4] Architecture : %s\n',   key_values.architecture_LCH4);
fprintf('[T-LCH4] n_stages     : %d\n',   key_values.n_stages_LCH4);
fprintf('[T-LCH4] nu_actual    : %.3f\n', key_values.nu_actual_LCH4);
fprintf('[T-LCH4] h_blade      : %.2f mm\n', key_values.h_blade_LCH4 * 1e3);
fprintf('[T-LCH4] HTR          : %.3f\n', key_values.HTR_LCH4);
fprintf('[T-LOx]  Architecture : %s\n',   key_values.architecture_LOx);
fprintf('[T-LOx]  n_stages     : %d\n',   key_values.n_stages_LOx);
fprintf('[T-LOx]  nu_actual    : %.3f\n', key_values.nu_actual_LOx);
fprintf('[T-LOx]  h_blade      : %.2f mm\n', key_values.h_blade_LOx * 1e3);
fprintf('[T-LOx]  HTR          : %.3f\n', key_values.HTR_LOx);
fprintf('===========================\n');
end

% -------------------------------------------------------------------------

function [x, res_norm] = newton2d(F, x0, lb, ftol, xtol, max_iter)
    % Newton-Raphson with backtracking line search for 2D system F(x) = 0.
    % lb: component-wise lower bounds enforced after every step.
    x = max(x0(:), lb);
    res_norm = inf;
    for k = 1:max_iter
        f = F(x);
        res_norm = norm(f);
        if res_norm < ftol; break; end

        % Finite-difference Jacobian (relative step, minimum 10 kPa)
        J = zeros(2,2);
        for j = 1:2
            h = max(abs(x(j)) * 1e-4, 1e4);
            xp = x; xp(j) = xp(j) + h;
            try
                fp = F(xp);
            catch
                fp = f * 2;
            end
            J(:,j) = (fp - f) / h;
        end

        if rcond(J) < 1e-12
            warning('newton2d: near-singular Jacobian at iter %d — stopping', k);
            break;
        end
        step = J \ (-f);

        % Backtracking: halve alpha until ||F|| strictly decreases
        alpha = 1.0;
        for ls = 1:12
            x_try = max(x + alpha * step, lb);
            try
                f_try = F(x_try);
                if norm(f_try) < res_norm; break; end
            catch
                % infeasible — shrink step
            end
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

    %% Tanks
    properties_flow.fuel.Tanks = properties_fuel;
    properties_flow.oxidizer.Tanks = properties_oxidizer;

    %% Pumps
  [properties_fuel,     key_values.P_Pump_LCH4] = sub_Pump_LCH4(inputs, properties_fuel);
properties_flow.fuel.Pump_LCH4                = properties_fuel;
inputs.P_turbine_LCH4_needed                  = key_values.P_Pump_LCH4;
inputs.N_shaft_LCH4                           = properties_fuel.RPM;
%inputs.D_mean_LCH4  = properties_fuel.D_impeller; % forward from pump

[properties_oxidizer, key_values.P_Pump_LOx]  = sub_Pump_LOx(inputs, properties_oxidizer);
properties_flow.oxidizer.Pump_LOx             = properties_oxidizer;
inputs.P_turbine_LOx_needed                   = key_values.P_Pump_LOx;
inputs.N_shaft_LOx                            = properties_oxidizer.RPM;
%inputs.D_mean_LOx    = properties_oxidizer.D_impeller; % forward from pump
    %% Cooling Channels (only fuel)
    [properties_fuel, key_values.delta_T_Cooling_Channels] = sub_Cooling_channels(inputs, properties_fuel);
    key_values.Q_dot_Cooling = inputs.Q_dot;
    properties_flow.fuel.Cooling_Channels = properties_fuel;

    %% Turbines — p_turbine_LCH4_out and p_turbine_LOx_out set by solver
    %% Turbines
[properties_fuel, kv_LCH4] = sub_Turbine_LCH4(inputs, properties_fuel);
key_values.delta_p_Turbine_LCH4 = kv_LCH4.delta_p;
key_values.W_Turbine_LCH4       = kv_LCH4.W_available;
% BUG FIX: forward geometry fields with _LCH4 suffix so cycle_solver fprintf can read them
key_values.architecture_LCH4    = kv_LCH4.architecture;
key_values.n_stages_LCH4        = kv_LCH4.n_stages;
key_values.nu_actual_LCH4       = kv_LCH4.nu_actual;
key_values.h_blade_LCH4         = kv_LCH4.h_blade;
key_values.HTR_LCH4             = kv_LCH4.HTR;
properties_flow.fuel.Turbine_LCH4 = properties_fuel;

[properties_fuel, kv_LOx] = sub_Turbine_LOx(inputs, properties_fuel);
key_values.delta_p_Turbine_LOx  = kv_LOx.delta_p;
key_values.W_Turbine_LOx        = kv_LOx.W_available;
% BUG FIX: forward geometry fields with _LOx suffix
key_values.architecture_LOx     = kv_LOx.architecture;
key_values.n_stages_LOx         = kv_LOx.n_stages;
key_values.nu_actual_LOx        = kv_LOx.nu_actual;
key_values.h_blade_LOx          = kv_LOx.h_blade;
key_values.HTR_LOx              = kv_LOx.HTR;
properties_flow.fuel.Turbine_LOx = properties_fuel;

W_LCH4 = kv_LCH4.W_available;
W_LOx  = kv_LOx.W_available;
P_LCH4 = inputs.P_turbine_LCH4_needed;
P_LOx  = inputs.P_turbine_LOx_needed;


    %% Injectors
    [properties_fuel, key_values.TC_Ingoing_fuel] = sub_Injector_CH4(inputs, properties_fuel);
    properties_flow.fuel.Injector = properties_fuel;

    [properties_oxidizer, key_values.TC_Ingoing_oxidizer] = sub_Injector_LOx(inputs, properties_oxidizer);
    properties_flow.oxidizer.Injector = properties_oxidizer;

    %% Thrust Chamber
    key_values.Thrust_Chamber = sub_Thrust_Chamber(inputs, properties_oxidizer, properties_fuel);
end

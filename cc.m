[inputs, properties_fuel, properties_oxidizer] = engine_inputs();

% --- CH4 pump ---
pump_CH4 = pump_efficiency_iteration( ...
    inputs.m_dot_fuel, ...
    inputs.delta_p_pump_LCH4, ...
    properties_fuel.rho, 'Methane', ...
    properties_fuel.p, properties_fuel.T);

% --- LOx pump ---
pump_LOx = pump_efficiency_iteration( ...
    inputs.m_dot_oxidizer, ...
    inputs.delta_p_pump_LOx, ...
    properties_oxidizer.rho, 'Oxygen', ...
    properties_oxidizer.p, properties_oxidizer.T);

% --- Update efficiency back into inputs for downstream use ---
inputs.eta_pump_LCH4 = pump_CH4.eta;
inputs.eta_pump_LOx  = pump_LOx.eta;

% Plot: max efficiency vs Omega_s for each pump
Omega_s_CH4 = [pump_CH4.Omega_s];
eta_CH4 = [pump_CH4.eta];
[uniq_Omega_s_CH4, ~, ic_CH4] = unique(Omega_s_CH4);
eta_max_CH4 = accumarray(ic_CH4(:), eta_CH4(:), [], @max);

Omega_s_LOx = [pump_LOx.Omega_s];
eta_LOx = [pump_LOx.eta];
[uniq_Omega_s_LOx, ~, ic_LOx] = unique(Omega_s_LOx);
eta_max_LOx = accumarray(ic_LOx(:), eta_LOx(:), [], @max);

figure;
plot(uniq_Omega_s_CH4, eta_max_CH4, 'b-', 'LineWidth', 2); hold on;
plot(uniq_Omega_s_LOx, eta_max_LOx, 'r-', 'LineWidth', 2);
xlabel('\Omega_s [-]'); ylabel('\eta [-]');
legend('CH4 pump', 'LOx pump');
title('Max Efficiency vs Specific Speed'); grid on;

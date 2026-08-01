function [inputs, properties_fuel, properties_oxidizer] = engine_inputs()

%% Physical constants
inputs.g0    = 9.8067;   % [m/s^2]
inputs.R_gas = 8314.4;   % [J/(kmol*K)]

%% Engine requirements
inputs.F_thrust_req = 40e3;   % [N]
inputs.p_CC_req     = 80e5;   % [Pa] 8 MPa

%% Efficiencies
inputs.eta_pump_LOx     = 0.70;
inputs.eta_pump_LCH4    = 0.68;
inputs.eta_turbine_LCH4 = 0.7;   % used when eta_turbine_mode = 'fixed'
inputs.eta_turbine_LOx  = 0.7;   % used when eta_turbine_mode = 'fixed'
inputs.eta_combustion   = 0.98;
inputs.eta_nozzle       = 0.96;

% Turbine efficiency mode:
%   'fixed' -> turbine files use eta_turbine_LCH4 / eta_turbine_LOx as-is
%   'balje' -> turbine files compute eta from Balje Ns correlation
inputs.eta_turbine_mode = 'balje';

%% Propellant properties
inputs.ROF = 3.09;
[T_CC, M, k] = get_cea_properties(inputs.p_CC_req, inputs.ROF);
inputs.T_CC_ideal          = T_CC;
inputs.Molar_mass_CC_ideal = M;
inputs.kappa_CC_ideal      = k;
inputs.v_e_ideal = sqrt(inputs.eta_nozzle * 2 * (k/(k-1)) * (inputs.R_gas/M) * T_CC * (1 - (101300/inputs.p_CC_req)^((k-1)/k)));

%% Mass flows
inputs.m_dot_tot = 13.154;    % add pressure terms later, when nozzle defined
inputs.m_dot_fuel = 3.233;
inputs.m_dot_oxidizer = 9.921;

%% Pressure losses
inputs.delta_p_inj_percent_LCH4 = 0.2;
inputs.delta_p_inj_percent_LOx  = 0.2;
inputs.delta_p_cooling_channels = 24.3e5;   % [Pa]
inputs.delta_p_feed             = 5e5;    % [Pa]
inputs.delta_p_partial          = inputs.delta_p_feed / 5;

%% Injector inlet pressures
inputs.p_fuel_injector_inlet = 88e5; % [Pa]
inputs.p_lox_injector_inlet  = 91.93e5; % [Pa]

%% Pressure chain (fuel side)
inputs.p_turbine_LOx_out  = 88e5;
inputs.p_turbine_LCH4_out = 103.4e5;

p_pump_LCH4_out_hc     = (inputs.p_turbine_LCH4_out + ...
    inputs.delta_p_cooling_channels + inputs.delta_p_partial) * 1.5;
p_tank_LOx     = 1.951e5;   % [Pa]
p_pump_LOx_out = 91.95e5;   % [Pa]

inputs.delta_p_pump_LOx = p_pump_LOx_out - p_tank_LOx;

p_tank_LCH4     = 1.986e5;    % [Pa]
p_pump_LCH4_out = 155.61e5;   % [Pa]

inputs.delta_p_pump_LCH4 = p_pump_LCH4_out - p_tank_LCH4;

%% Turbine velocity ratio targets
inputs.nu_target_LCH4  = 0.70;
inputs.nu_target_LOx   = 0.55;

%% Max rotor diameter per stage [m]
inputs.D_mean_max_LCH4 = 0.080;
inputs.D_mean_max_LOx  = 0.20;

%% Cooling
inputs.Q_dot = 4e6;   % [W]

%% Initial fuel properties
properties_fuel.T   = 110;
properties_fuel.p   = 1.986e5;
properties_fuel.rho = py.CoolProp.CoolProp.PropsSI('D','P',properties_fuel.p,'T',properties_fuel.T,'Methane');
properties_fuel.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS','T',properties_fuel.T,'P',properties_fuel.p,'Methane');

%% Initial oxidizer properties
properties_oxidizer.T   = 90;
properties_oxidizer.p   = 1.951e5;
properties_oxidizer.rho = py.CoolProp.CoolProp.PropsSI('D','P',properties_oxidizer.p,'T',properties_oxidizer.T,'Oxygen');
properties_oxidizer.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS','T',properties_oxidizer.T,'P',properties_oxidizer.p,'Oxygen');

end
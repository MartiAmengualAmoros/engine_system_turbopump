function [inputs, properties_fuel, properties_oxidizer] = engine_inputs()

    %% Physical constants
    inputs.g0 = 9.8067;                     % [m/s^2]
    inputs.R_gas = 8314.4;                  % [J/(kmol·K)]

    %% Engine requirements
    inputs.F_thrust_req = 40e3;             % [N]
    inputs.p_CC_req = 80e5;                 % [Pa]

    %% Efficiencies
    inputs.eta_pump_LOx = 0.70;             % (Ask teams and revise)
    inputs.eta_pump_LCH4 = 0.60;            % (Ask teams and revise)
    inputs.eta_turbine_LOx = 0.65;          % (Ask teams and revise)
    inputs.eta_turbine_LCH4 = 0.65;         % (Ask teams and revise)
    inputs.eta_combustion = 0.98;           % (Ask teams and revise)
    inputs.eta_nozzle = 0.96;               % (Ask teams and revise)

    %% Propellant properties
    inputs.ROF = 3.09;                       % (Research and ask Thrust chamber)  stoch
    [T_CC, M, k] = get_cea_properties(inputs.p_CC_req, inputs.ROF);
    inputs.T_CC_ideal = T_CC;
    inputs.Molar_mass_CC_ideal = M;
    inputs.kappa_CC_ideal = k;
    inputs.v_e_ideal = sqrt(inputs.eta_nozzle * 2 * (inputs.kappa_CC_ideal/(inputs.kappa_CC_ideal-1)) * (inputs.R_gas/inputs.Molar_mass_CC_ideal) * inputs.T_CC_ideal * (1 - (101300/inputs.p_CC_req)^((inputs.kappa_CC_ideal-1)/inputs.kappa_CC_ideal)));

    %% Get necessary massflow (at given ROF):
    inputs.m_dot_tot = inputs.F_thrust_req / inputs.v_e_ideal;    % add pressure terms later, when nozzle defined
    inputs.m_dot_fuel = inputs.m_dot_tot / (inputs.ROF + 1);
    inputs.m_dot_oxidizer = inputs.m_dot_tot - inputs.m_dot_fuel;

    %% Pressures
    % The team's math translates to a 7% fuel loss and an 8% oxidizer loss
    inputs.delta_p_inj_percent_LOx = 0.10;      
    inputs.delta_p_inj_percent_LCH4 = 0.07;     
    
    inputs.delta_p_cooling_channels = 13e5;   % Team's 16.9 bar drop
    inputs.delta_p_feed = 5e5;              
    inputs.delta_p_partial = inputs.delta_p_feed / 5;   
    inputs.delta_p_pump_LCH4 = 6e6;         

    % SEPARATE Target Inlet Pressures for the Injectors (Working backward from 80 bar)
    inputs.p_fuel_injector_inlet = inputs.p_CC_req / (1 - inputs.delta_p_inj_percent_LCH4);
    inputs.p_lox_injector_inlet = inputs.p_CC_req / (1 - inputs.delta_p_inj_percent_LOx);

    % LOx Line End: The LOx pump feeds the LOx injector directly
    p_tank_LOx = 2e5;                       
    inputs.delta_p_pump_LOx = inputs.p_lox_injector_inlet - p_tank_LOx + inputs.delta_p_partial;
 %% Turbine geometry inputs
    inputs.D_mean_LCH4    = 0.05;    % [m]   Initial estimate — update from impeller sizing
    inputs.nu_target_LCH4 = 0.45;    % [-]   Blade speed ratio target (impulse: 0.35–0.50)
    inputs.D_mean_LOx     = 0.06;    % [m]   Initial estimate
    inputs.nu_target_LOx  = 0.45;    % [-]   Blade speed ratio target
    % Turbine exit pressure — forced by injector inlet requirement.
    % Injector model drops p by (1 - delta_p_inj_percent), so to land exactly at p_CC: p_in = p_CC / (1 - pct)
    inputs.p_turbine_exit = inputs.p_CC_req / (1 - inputs.delta_p_inj_percent_LCH4);

    % LOx pump rise derived from turbine exit pressure so the chain always closes:
    % p_exit_LOx = p_tank + delta_p_pump_LOx - delta_p_partial = p_turbine_exit
    p_tank_LOx = 2e5;                       % [Pa]  oxidizer tank pressure
    inputs.delta_p_pump_LOx = inputs.p_turbine_exit - p_tank_LOx + inputs.delta_p_partial;

    %% Cooling assumptions
    % Bumped up slightly to hit the team's 540 K target
    inputs.Q_dot = 4e6;                   % [J/s]

    %% Fuel properties (input here initial properties
    properties_fuel.T = 110;                % [K]       % (Design choice, fairly easy, just a quick research)
    properties_fuel.p = 2e5;                % [Pa]
    properties_fuel.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties_fuel.p, 'T', properties_fuel.T, 'Methane');
    properties_fuel.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties_fuel.T, 'P', properties_fuel.p, 'Methane');

    %% Fuel properties (input here initial properties);
    properties_oxidizer.T = 90;             % [K]
    properties_oxidizer.p = 2e5;            % [Pa]
    properties_oxidizer.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties_oxidizer.p, 'T', properties_oxidizer.T, 'Oxygen');
    properties_oxidizer.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties_oxidizer.T, 'P', properties_oxidizer.p, 'Oxygen');

end
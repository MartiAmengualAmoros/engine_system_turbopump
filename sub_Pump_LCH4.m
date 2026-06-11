function [properties, key_values] = sub_Pump_LCH4(inputs, properties)

    % 0. Capture Inlet States (Needed for Hydraulics & NPSH)
    p_in = properties.p;
    T_in = properties.T;
    rho_in = py.CoolProp.CoolProp.PropsSI('D', 'T', T_in, 'P', p_in, 'Methane');

    % 1. Define Pressures
    p_intermediate = p_in + inputs.delta_p_pump_LCH4;
    p_final = p_intermediate - inputs.delta_p_partial;
    
    % 2. Calculate Inlet Thermodynamic State
    h_in = py.CoolProp.CoolProp.PropsSI('H', 'T', T_in, 'P', p_in, 'Methane');
    s_in = py.CoolProp.CoolProp.PropsSI('S', 'T', T_in, 'P', p_in, 'Methane');
    
    % 3. Calculate Ideal (Isentropic) Compression
    h_out_s = py.CoolProp.CoolProp.PropsSI('H', 'P', p_intermediate, 'S', s_in, 'Methane');
    
    % 4. Calculate Real Compression
    h_out_intermediate = h_in + (h_out_s - h_in) / inputs.eta_pump_LCH4;
    
    % 5. Calculate Required Shaft Power [W]
    P_pump = inputs.m_dot_fuel * (h_out_intermediate - h_in);

   %% --- HYDRAULICS, CAVITATION (NPSH), AND SIZING ---
    
    p_vapor = py.CoolProp.CoolProp.PropsSI('P', 'T', T_in, 'Q', 0, 'Methane');
    NPSH = (p_in - p_vapor) / (rho_in * inputs.g0);
    Head_Total = inputs.delta_p_pump_LCH4 / (rho_in * inputs.g0);
    Q = inputs.m_dot_fuel / rho_in;
    
    % --- 2-STAGE PUMP WITH INDUCER ARCHITECTURE ---
    n_stages = 4;
    Head_per_stage = Head_Total / n_stages;
    
    % Inducer Suction Specific Speed
    N_ss = 300; 
    
    % Calculate RPM based on the inducer limits
    N_max_rpm = N_ss * (NPSH^0.75) / sqrt(Q);
    Omega = N_max_rpm * (pi / 30);
    
    % Calculate Universal Specific Speed PER STAGE
    N_s_universal = (Omega * sqrt(Q)) / (inputs.g0 * Head_per_stage)^0.75;

    %% --- OUTPUTS ---
    % Update standard thermodynamic properties for downstream components
    properties.p = p_final;
    properties.T = py.CoolProp.CoolProp.PropsSI('T', 'P', properties.p, 'H', h_out_intermediate, 'Methane');
    properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Methane');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Methane');

    % Attach sizing parameters SILENTLY to the properties struct
    properties.NPSH = NPSH;
    properties.Head_Total = Head_Total;
    properties.Head_Per_Stage = Head_per_stage;
    properties.Volumetric_Q = Q;
    properties.RPM = N_max_rpm;
    properties.Specific_Speed_Ns_Universal = N_s_universal;

    % Output MUST remain a scalar to not break shared scripts
    key_values = P_pump;
end
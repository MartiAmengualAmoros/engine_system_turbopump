function [properties, key_values] = sub_Pump_LOx(inputs, properties)

    % 0. Capture Inlet States (Needed for Hydraulics & NPSH)
    p_in = properties.p;
    T_in = properties.T;
    rho_in = py.CoolProp.CoolProp.PropsSI('D', 'T', T_in, 'P', p_in, 'Oxygen');

    % 1. Define Pressures
    p_intermediate = p_in + inputs.delta_p_pump_LOx;
    p_final = p_intermediate - inputs.delta_p_partial;
    
    % 2. Calculate Inlet Thermodynamic State
    h_in = py.CoolProp.CoolProp.PropsSI('H', 'T', T_in, 'P', p_in, 'Oxygen');
    s_in = py.CoolProp.CoolProp.PropsSI('S', 'T', T_in, 'P', p_in, 'Oxygen');
    
    % 3. Calculate Ideal (Isentropic) Compression
    h_out_s = py.CoolProp.CoolProp.PropsSI('H', 'P', p_intermediate, 'S', s_in, 'Oxygen');
    
    % 4. Calculate Real Compression
    h_out_intermediate = h_in + (h_out_s - h_in) / inputs.eta_pump_LOx;
    
    % 5. Calculate Required Shaft Power [W]
    P_pump = inputs.m_dot_oxidizer * (h_out_intermediate - h_in);

   %% --- HYDRAULICS, CAVITATION (NPSH), AND SIZING ---
    
    p_vapor = py.CoolProp.CoolProp.PropsSI('P', 'T', T_in, 'Q', 0, 'Oxygen');
    NPSH = (p_in - p_vapor) / (rho_in * inputs.g0);
    Head = inputs.delta_p_pump_LOx / (rho_in * inputs.g0);
    Q = inputs.m_dot_oxidizer / rho_in;
    
    N_ss = 150; 
    N_max_rpm = N_ss * (NPSH^0.75) / sqrt(Q);
    
    % --- UNIVERSAL SPECIFIC SPEED CALCULATION ---
    % 1. Convert RPM to angular velocity (rad/s)
    Omega = N_max_rpm * (pi / 30);
    
    % 2. Calculate true dimensionless Specific Speed
    N_s_universal = (Omega * sqrt(Q)) / (inputs.g0 * Head)^0.75;

    %% --- OUTPUTS ---
    % Update standard thermodynamic properties for downstream components
    properties.p = p_final;
    properties.T = py.CoolProp.CoolProp.PropsSI('T', 'P', properties.p, 'H', h_out_intermediate, 'Oxygen');
    properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Oxygen');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Oxygen');

    % Attach sizing parameters SILENTLY to the properties struct
    properties.NPSH = NPSH;
    properties.Head = Head;
    properties.Volumetric_Q = Q;
    properties.RPM = N_max_rpm;
    properties.Specific_Speed_Ns_Universal = N_s_universal; % Updated name

    % Output MUST remain a scalar to not break shared scripts
    key_values = P_pump;
end
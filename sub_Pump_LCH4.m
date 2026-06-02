function [properties, key_values] = sub_Pump_LCH4(inputs, properties)

    % REAL-FLUID CRYOGENIC PUMP MODEL (METHANE)
    % Replaces the incompressible assumption (W = m*dp/rho) with a rigorous 
    % isentropic enthalpy method suitable for supercritical/high-pressure LCH4.

    % 1. Define Pressures
    % The pump rotor does physical work to raise the pressure to p_intermediate.
    % The delta_p_partial represents friction/valve losses *after* the pump discharge.
    p_intermediate = properties.p + inputs.delta_p_pump_LCH4;
    p_final = p_intermediate - inputs.delta_p_partial;
    
    % 2. Calculate Inlet Thermodynamic State
    % Extract the exact enthalpy and entropy of the LCH4 entering the pump.
    h_in = py.CoolProp.CoolProp.PropsSI('H', 'T', properties.T, 'P', properties.p, 'Methane');
    s_in = py.CoolProp.CoolProp.PropsSI('S', 'T', properties.T, 'P', properties.p, 'Methane');
    
    % 3. Calculate Ideal (Isentropic) Compression
    % What the enthalpy would be if the pump were 100% efficient.
    h_out_s = py.CoolProp.CoolProp.PropsSI('H', 'P', p_intermediate, 'S', s_in, 'Methane');
    
    % 4. Calculate Real Compression
    % Apply the hydraulic efficiency to find the actual enthalpy at pump discharge.
    % Because pumps *consume* work, efficiency divides the enthalpy rise.
    h_out_intermediate = h_in + (h_out_s - h_in) / inputs.eta_pump_LCH4;
    
    % 5. Calculate Required Shaft Power [W]
    P_pump = inputs.m_dot_fuel * (h_out_intermediate - h_in);
    
    % 6. Apply Feed Line Pressure Drop
    % Pressure drop through pipes and valves (delta_p_partial) is typically 
    % modeled as an isenthalpic process (enthalpy remains constant).
    h_out_final = h_out_intermediate;
    
    % 7. Update Properties Struct for Downstream Components
    properties.p = p_final;
    
    % Use CoolProp to find the true temperature based on the final pressure and enthalpy.
    % This replaces the old (T + P*(1-eta)/(m*c_p)) approximation.
    properties.T = py.CoolProp.CoolProp.PropsSI('T', 'P', properties.p, 'H', h_out_final, 'Methane');
    
    % Update density and heat capacity at the new high-pressure state.
    properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Methane');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Methane');

    % Return the power required so the cycle_solver can balance the LCH4 turbine
    key_values = P_pump;
end
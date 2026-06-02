function [properties, key_values] = sub_Pump_LOx(inputs, properties)

    % TODO: Research validity of Formulary equations in real life
    %       Polytropic exponent?
    %       Add 1-eta energy to temperature?
    %       Research more realistic equations/transformations
    
    p_intermediate = properties.p + inputs.delta_p_pump_LOx;
    properties.p = p_intermediate - inputs.delta_p_partial;

    % Isentropic specific work for the pump
    delta_h_isen = inputs.delta_p_pump_LOx / properties.rho;
    % Actual specific work with pump efficiency
    delta_h_actual = delta_h_isen / inputs.eta_pump_LOx;
    P_pump = inputs.m_dot_oxidizer * delta_h_actual;

    % Excess enthalpy due to inefficiency heats the fluid
    delta_h_loss = delta_h_actual - delta_h_isen;
    properties.T = properties.T + delta_h_loss / properties.c_p;

    properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Oxygen');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Oxygen');

    key_values = P_pump;
end
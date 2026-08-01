function [properties, key_values] = sub_Cooling_channels(inputs, properties)

% TODO: Research/ask for a more accurate estimate of Q_dot
% Have Q_dot dependent on other quantities? percentage of energy in thrust chamber? Research
% --> Check correct use of CoolProp over state transition into supercritical
% Research validity of isobaric assumption
% --> Research Heat going other places

% Guard: if inlet pressure is non-physical, return passthrough state
p_min_valid = 1e3; % [Pa] absolute minimum physical pressure (1 kPa)
if properties.p <= p_min_valid
    % Return unchanged properties and zero delta_T — solver will self-correct
    key_values = 0;
    return
end

% m_dot * h_in + Q_dot = m_dot * h_out
try
    h_in  = py.CoolProp.CoolProp.PropsSI('H', 'P', properties.p, 'T', properties.T, 'Methane');
    h_out = h_in + inputs.Q_dot / inputs.m_dot_fuel;

    T_in         = properties.T;
    properties.T = py.CoolProp.CoolProp.PropsSI('T', 'P', properties.p, 'H', h_out, 'Methane');
    delta_T      = properties.T - T_in;

    properties.p = properties.p - inputs.delta_p_cooling_channels - inputs.delta_p_partial;

    % Guard against negative pressure after drop
    if properties.p <= p_min_valid
        properties.p = p_min_valid;
    end

    properties.rho = py.CoolProp.CoolProp.PropsSI('D',      'P', properties.p, 'T', properties.T, 'Methane');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Methane');

    key_values = delta_T;
catch
    % CoolProp failed on an unphysical Newton probe — return passthrough
    key_values = 0;
end
end
function [properties, key_values] = sub_Injector_CH4(inputs, properties)

    % TODO: Check ALL assumptions, ask porous injector team
    %       wtf is a porous injector?

    properties.T = properties.T; % ?

    properties.p = properties.p * (1 - inputs.delta_p_inj_percent_LCH4);

    properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Methane');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Methane');

    key_values.T_ingoing = properties.T;
    key_values.p_ingoing = properties.p;
end
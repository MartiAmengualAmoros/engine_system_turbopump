% sub_Pump_LCH4.m
function [properties, key_values] = sub_Pump_LCH4(inputs, properties)

    % 1. Thermodynamics & Hydraulics (Merged)
    p_in = properties.p;
    T_in = properties.T;
    
    rho_in = double(py.CoolProp.CoolProp.PropsSI('D', 'T', T_in, 'P', p_in, 'Methane'));
    mu_in  = double(py.CoolProp.CoolProp.PropsSI('VISCOSITY', 'T', T_in, 'P', p_in, 'Methane'));
    h_in   = double(py.CoolProp.CoolProp.PropsSI('H', 'T', T_in, 'P', p_in, 'Methane'));
    s_in   = double(py.CoolProp.CoolProp.PropsSI('S', 'T', T_in, 'P', p_in, 'Methane'));

    p_intermediate = p_in + inputs.delta_p_pump_LCH4;
    p_final = p_intermediate - inputs.delta_p_partial;
    
    h_out_s = double(py.CoolProp.CoolProp.PropsSI('H', 'P', p_intermediate, 'S', s_in, 'Methane'));
    h_out_intermediate = h_in + (h_out_s - h_in) / inputs.eta_pump_LCH4;
    
    P_pump = inputs.m_dot_fuel * (h_out_intermediate - h_in);

    p_vapor = double(py.CoolProp.CoolProp.PropsSI('P', 'T', T_in, 'Q', 0, 'Methane'));
    NPSH = (p_in - p_vapor) / (rho_in * inputs.g0);
    
    Head_Total = inputs.delta_p_pump_LCH4 / (rho_in * inputs.g0);
    Q = inputs.m_dot_fuel / rho_in;
    
    n_stages = 4;
    Head_per_stage = Head_Total / n_stages;
    
    N_ss = 300; 
    N_max_rpm = N_ss * (NPSH^0.75) / sqrt(Q);
    Omega = N_max_rpm * (pi / 30);

    %% 2. Geometric Sizing  

    % ------------- INPUTS --------------

    use_inducer = true;

    eta_inducer = 0.60;
    delta_p_inducer = 1.0e5;        % [Pa] initial guess
    alpha1_inducer = 0.0;         % [rad] effective prewhirl at impeller inlet
    c1m_factor_inducer = 0.95;    % modifies meridional inlet velocity
    w1u_factor_inducer = 0.50;    % reduces relative tangential inlet component

    eta_hyd = inputs.eta_pump_LCH4; 
    alpha1 = 0.0;     
    phi = 0.13;       
    psi = 0.54;       
    R = 0.60;         
    delta_h = 0.30;   


    % ------------- IMPELLER --------------


    W_stage = inputs.g0 * Head_per_stage / eta_hyd; 
    u2 = sqrt(W_stage / psi); 
    D2 = 60 * u2 / (pi * N_max_rpm); 
    c1m = u2 * phi; 
    A1 = Q / c1m; 

    delta_t = sqrt(delta_h^2 + 4 * A1 / (pi * D2^2)); 
    D1h = delta_h * D2; 
    D1t = delta_t * D2; 
    b1 = (D1t - D1h) / 2; 

    term_xi = 1 - (2 * psi / phi^2) * (R + psi / 2 - 1) + tan(alpha1)^2 * (1 - delta_t^2) - (2 * psi * delta_t / phi) * tan(alpha1); 
    xi = sqrt(term_xi); 

    alpha2 = atan( psi / (xi * phi) + (delta_t / xi) * tan(alpha1) ); 
    beta1 = atan( delta_t / phi - tan(alpha1) ); 
    beta2 = atan( (1 / (phi * xi)) * (1 - psi) - (delta_t / xi) * tan(alpha1) ); 
    
    c2m = u2 * xi * phi;
    A2 = Q / c2m; 
    b2 = A2 / (pi * D2); 

    w1 = u2 * phi * sqrt(1 + tan(alpha1)^2);
    w1u = u2 * phi * tan(beta1);


    N_s_universal = (Omega * sqrt(Q)) / (inputs.g0 * Head_per_stage)^0.75; 


        % --- Effective inlet conditions after inducer ---
    alpha1_eff = alpha1;
    c1m_eff = c1m;
    w1u_eff = w1u;
    p_in_eff = p_in;

    if use_inducer
        alpha1_eff = alpha1_inducer;
        c1m_eff = c1m * c1m_factor_inducer;
        w1u_eff = w1u * w1u_factor_inducer;
        p_in_eff = p_in + delta_p_inducer;
    end

    % Recompute available NPSH at impeller inlet
    NPSH = (p_in_eff - p_vapor) / (rho_in * inputs.g0);

    m = 0.08  * (1 / Q)^0.15 * (0.85 / N_s_universal)^0.06; 

    eta_hyd = 1 - 0.065 * (1 / Q)^m - 0.23 * (0.3 - log10(2.3 * N_s_universal))^2 * (1 / Q)^0.05; 


    lambda_c = -0.103 * log10(N_s_universal) + 1.1855;  

    lambda_w = 0.2144 * exp(0.1745 * N_s_universal); 

    NPSH_R = lambda_c * (c1m_eff^2 / (2 * inputs.g0)) + lambda_w * (w1u_eff^2 / (2 * inputs.g0)); 





    % 3. Outputs
    properties.p = p_final;
    properties.T = py.CoolProp.CoolProp.PropsSI('T', 'P', properties.p, 'H', h_out_intermediate, 'Methane');
    properties.rho = py.CoolProp.CoolProp.PropsSI('D', 'P', properties.p, 'T', properties.T, 'Methane');
    properties.c_p = py.CoolProp.CoolProp.PropsSI('CPMASS', 'T', properties.T, 'P', properties.p, 'Methane');

    properties.NPSH = NPSH;
    properties.Head_Total = Head_Total;
    properties.Head_Per_Stage = Head_per_stage;
    properties.Volumetric_Q = Q;
    properties.RPM = N_max_rpm;
    properties.Specific_Speed_Ns_Universal = N_s_universal;
    
    properties.n_stages = n_stages;
    properties.u2 = u2;
    properties.D2 = D2;
    properties.D1t = D1t;
    properties.b2 = b2;
    properties.delta_t = delta_t;


    properties.eta_hyd = eta_hyd;
    properties.NPSH_R = NPSH_R;


    % Attach angle geometry for the solver to print
    properties.alpha1 = alpha1;
    properties.beta1  = beta1;
    properties.alpha2 = alpha2;
    properties.beta2  = beta2;
    properties.beta2B = beta2B; 

    key_values = P_pump;
end
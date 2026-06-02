function pumps = pump_efficiency_iteration(m_dot, delta_p, rho, propellant_name, p_in, T_in, varargin)
% Iterates pump efficiency as a function of rotational specific speed
% and inlet diameter until convergence.
%
% INPUTS:
%   m_dot           - mass flow rate [kg/s]
%   delta_p         - pump pressure rise [Pa] or vector of values for sweep
%   rho             - inlet fluid density [kg/m^3]
%   propellant_name - string: 'Methane' or 'Oxygen' for CoolProp
%   p_in            - inlet pressure [Pa]
%   T_in            - inlet temperature [K]
%   Optional name/value pairs:
%     'Omega_s_range' - vector of dimensionless specific speeds
%     'psi_range'     - vector of head coefficients
%     'phi_range'     - vector of flow coefficients
%     'Plot'          - true/false to plot efficiency results
%     'PlotX'         - x-axis variable for plot: 'Omega_s', 'psi_target', 'phi_target', 'delta_p', 'n_rpm', or 'D_in'
%
% OUTPUT:
%   pumps struct array with converged parameters per sweep case

%% ---- Constants -------------------------------------------------------
g0 = 9.8067;                % [m/s^2]

%% ---- Derived hydraulic quantities ------------------------------------
V_dot = m_dot / rho;                    % volumetric flow [m^3/s]

% Allow delta_p to be a vector for sweep calculations
% delta_p is the pump pressure rise [Pa] for each design case.
delta_p_list = delta_p(:)';

%% ---- Fluid viscosity via CoolProp ------------------------------------
mu = py.CoolProp.CoolProp.PropsSI('VISCOSITY', 'P', p_in, 'T', T_in, propellant_name);
nu = mu / rho;                              % kinematic viscosity [m^2/s]

%% ---- Iteration setup -------------------------------------------------
p = inputParser;
addParameter(p,'Omega_s_range',0.50);
addParameter(p,'psi_range',0.47);
addParameter(p,'phi_range',0.15);
addParameter(p,'Plot',false);
addParameter(p,'PlotX','Omega_s');
addParameter(p,'tol',1e-5);
addParameter(p,'max_iter',200);
parse(p,varargin{:});

Omega_s_range = p.Results.Omega_s_range(:)';
psi_range     = p.Results.psi_range(:)';
phi_range     = p.Results.phi_range(:)';
plot_on       = p.Results.Plot;
plot_x        = validatestring(p.Results.PlotX,{'Omega_s','psi_target','phi_target','delta_p','n_rpm','D_in'});
tol           = p.Results.tol;
max_iter      = p.Results.max_iter;

fprintf('\n--- Pump Efficiency Sweep: %s ---\n', propellant_name);
fprintf('%4s | %6s | %6s | %8s | %8s | %8s | %8s | %10s | %8s\n', ...
    'run','Omega','psi','phi','n [rpm]','D_in','D2','Re','eta');

% Pre-allocate result container
pumps = struct([]);
run_idx = 0;

% Allow delta_p to be a vector (sweep over heads).
delta_p_list = delta_p(:)';

for ip_dp = 1:numel(delta_p_list)
    delta_p_val = delta_p_list(ip_dp);
    delta_H = delta_p_val / (rho * g0);

    for iOmega = 1:numel(Omega_s_range)
        Omega_s_target = Omega_s_range(iOmega);
        for ipsi = 1:numel(psi_range)
            psi_target = psi_range(ipsi);
            for iphi = 1:numel(phi_range)
                phi_target = phi_range(iphi);

                % Each combination gets its own iterative convergence
                eta_assumed  = 0.70;          % initial assumed efficiency
                converged    = false;

                for iter = 1:max_iter

                    %% Step 1: Shaft power and actual head at current efficiency
                    W_pump    = (m_dot * delta_p_val) / (rho * eta_assumed);   % shaft power [W]

                    %% Step 2: Compute shaft speed from specific speed target
                    % Omega_s = n * sqrt(Q) / (g H)^(3/4) in SI units (rad/s)
                    omega = Omega_s_target * (g0 * delta_H)^(3/4) / sqrt(V_dot);  % [rad/s]
                    n_rpm = omega * 60 / (2 * pi);                                  % [rpm]

                    %% Step 3: Impeller tip speed and diameter from head coefficient
                    % psi = g H / u2^2  =>  u2 = sqrt(g H / psi)
                    u2 = sqrt(g0 * delta_H / psi_target);
                    D2 = 2 * u2 / omega;                  % impeller tip diameter [m]

                    %% Step 4: Inlet eye area from flow coefficient
                    % phi = Q / (A_in * u2)
                    A_in = V_dot / (phi_target * u2);
                    D_in = sqrt(4 * A_in / pi);            % inlet eye diameter [m]

                    %% Step 5: Reynolds number at inlet
                    c_in = V_dot / A_in;                   % mean axial velocity at inlet [m/s]
                    Re = c_in * D_in / nu;                 % Reynolds number [-]

                    %% Step 6: Empirical pump efficiency correlation
                    % Peak efficiency depends on specific speed and flow coefficient.
                    eta_peak = 0.88;
                    eta_Omega = 1 - 1.2 * ((Omega_s_target - 0.55) / 0.35).^2;
                    eta_phi = 1 - 1.0 * ((phi_target - 0.15) / 0.12).^2;
                    eta_Omega = max(0.55, min(1.00, eta_Omega));
                    eta_phi = max(0.50, min(1.00, eta_phi));

                    if Re > 0
                        f_Re = 1 - 0.18 * min(1, (1e6 / Re)^0.2);
                    else
                        f_Re = 0.80;
                    end

                    f_disk = 1 - 0.02 * (D2 / D_in)^2 / (Re^0.2 + 1e-12);
                    f_disk = max(0.82, min(1.00, f_disk));

                    eta_computed = eta_peak * eta_Omega * eta_phi * f_Re * f_disk;
                    eta_computed = max(0.45, min(0.92, eta_computed));

                    %% Step 7: Convergence check and relaxation update
                    if abs(eta_computed - eta_assumed) < tol
                        converged = true;
                        break;
                    end
                    relax       = 0.4;
                    eta_assumed = eta_assumed + relax * (eta_computed - eta_assumed);

                end % iter

                % Store results for this run
                run_idx = run_idx + 1;
                pumps(run_idx).eta          = eta_computed;       % converged efficiency [-]
                pumps(run_idx).n_rpm        = n_rpm;              % shaft speed [rpm]
                pumps(run_idx).omega        = omega;              % shaft speed [rad/s]
                pumps(run_idx).Omega_s      = Omega_s_target;     % dimensionless specific speed [-]
                pumps(run_idx).D_in         = D_in;               % inlet eye diameter [m]
                pumps(run_idx).D2           = D2;                 % impeller tip diameter [m]
                pumps(run_idx).u2           = u2;                 % tip speed [m/s]
                pumps(run_idx).Re           = Re;                 % Reynolds number [-]
                pumps(run_idx).W_shaft      = W_pump;             % shaft power [W]
                pumps(run_idx).V_dot        = V_dot;              % volumetric flow [m^3/s]
                pumps(run_idx).delta_H      = delta_H;            % developed head [m]
                pumps(run_idx).c_in         = c_in;               % inlet axial velocity [m/s]
                pumps(run_idx).psi_target   = psi_target;
                pumps(run_idx).phi_target   = phi_target;
                pumps(run_idx).delta_p      = delta_p_val;
                pumps(run_idx).converged    = converged;
                pumps(run_idx).iterations   = iter;

                fprintf('%4d | %6.3f | %6.3f | %8.3f | %8.1f | %8.4f | %8.4f | %10.2e | %8.4f\n', ...
                    run_idx, Omega_s_target, psi_target, phi_target, n_rpm, D_in, D2, Re, eta_computed);

            end % phi
        end % psi
    end % Omega
end % delta_p

% Plot efficiency if requested
if plot_on
    numRuns = numel(pumps);
    xvals = zeros(numRuns,1);
    yvals = zeros(numRuns,1);
    labels = cell(numRuns,1);
    for ii = 1:numRuns
        xvals(ii) = pumps(ii).(plot_x);
        yvals(ii) = pumps(ii).eta;
        labels{ii} = sprintf('Omega_s=%.3f psi=%.3f phi=%.3f dp=%.0f', ...
            pumps(ii).Omega_s, pumps(ii).psi_target, pumps(ii).phi_target, pumps(ii).delta_p);
    end
    uniqueLabels = unique(labels, 'stable');

    figure;
    hold on;
    for jj = 1:numel(uniqueLabels)
        idx = strcmp(labels, uniqueLabels{jj});
        [x_sorted, sort_idx] = sort(xvals(idx));
        y_sorted = yvals(idx);
        y_sorted = y_sorted(sort_idx);
        plot(x_sorted, y_sorted, '-o', 'DisplayName', uniqueLabels{jj}, 'LineWidth', 1.4);
    end
    grid on;
    xlabel(plot_x, 'Interpreter', 'none');
    ylabel('Efficiency [-]');
    title(sprintf('Pump Efficiency vs %s', plot_x), 'Interpreter', 'none');
    legend('Location','best');
    hold off;
end

% If only one result, return scalar struct for compatibility
if numel(pumps) == 1
    pumps = pumps(1);
end

end
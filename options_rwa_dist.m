function [struct_out, srp_info] = options_rwa_dist(sim_time, fs, rpm_start)
%% =====================================================================
%  DYNAMIC SIMULATION OPTIONS: SRP-DRIVEN RWA MOMENTUM STORAGE
%  =====================================================================
%  Inputs (all optional):
%  - sim_time  : [Double] Total simulation duration in seconds (default: 3600)
%  - fs        : [Double] Sampling frequency in Hz (default: 1000)
%  - rpm_start : [Double] Initial RWA spin rate in RPM (default: 0)
%
%  Outputs:
%  - struct_out : Simulation options struct (fields: t, forces)
%  - srp_info   : Diagnostic struct with SRP force, torque, and RWA scaling

    if nargin < 1 || isempty(sim_time),  sim_time  = 3600; end
    if nargin < 2 || isempty(fs),        fs        = 1000; end
    if nargin < 3 || isempty(rpm_start), rpm_start = 0;    end

    %% 1. Merge FEM Models and Extract Global Mass & Geometry
    SV = fem.merge_fem_models(config_bus(), config_Iso_Obench(), config_Payload());

    M_total = SV.M_total;   % Total spacecraft mass [kg]
    I_total = SV.I_total;   % Total spacecraft rotational inertia [kg*m^2]
    x_cg    = SV.CG(1);     % Spacecraft CG x-coordinate [m]

    % Extract nodal coordinates to determine bottom-face width & centerline
    all_x = [SV.nodes.x];
    x_centerline = 0.5 * (min(all_x) + max(all_x));

    % Calculate the width of the bottom of the vehicle
    width_bottom = max(all_x) - min(all_x);


    %% 2. Solar Radiation Pressure (SRP) Force & Disturbance Torque
    S_0   = 1361;             % Solar constant at 1 AU [W/m^2]
    c     = 2.99792458e8;     % Speed of light [m/s]
    P_srp = S_0 / c;          % Nominal solar radiation pressure [N/m^2]

    C_r      = 1.4;           % Surface reflectivity coeff
    sv_depth = 2;  % Assumed out-of-plane depth for 2D FEM [m]
    A_bottom = width_bottom * sv_depth; % Projected SRP area [m^2]

    % Net SRP acts 1.0 m to the right of the vehicle centerline
    x_cp      = x_centerline + 1.0; % Center of pressure x-coordinate [m]
    lever_arm = x_cp - x_cg;        % Moment arm about CG [m]

    F_srp = P_srp * A_bottom * C_r; % Net SRP force [N]
    T_srp = F_srp * lever_arm;      % Net SRP disturbance torque [N*m]

    %% 3. RWA Sizing & Time-Scaled Spin Rate Calculation
    % Detect RWA nodes from the merged model
    node_ids = string({SV.nodes.id});
    rwa_mask = contains(node_ids, "RWA", "IgnoreCase", true);
    rwa_ids  = node_ids(rwa_mask);
    num_rwas = numel(rwa_ids);

    if num_rwas == 0
        error('No RWA nodes found in merged SV model.');
    end

    % RWA hardware parameters
    unbalance     = 0.05;                % Static unbalance [kg*m]
    rotor_inertia = 0.5 * 10 * 0.2^2;    % Polar inertia [kg*m^2]

    % Divide required counter-torque equally across all RWAs
    T_wheel     = -T_srp / num_rwas;              % Torque per wheel [N*m]
    alpha_wheel = T_wheel / rotor_inertia;        % Angular accel [rad/s^2]

    % Scale total speed change linearly with simulation time
    delta_omega = alpha_wheel * sim_time;         % [rad/s]
    delta_rpm   = delta_omega * (60 / (2 * pi));  % [RPM]
    rpm_end     = rpm_start + delta_rpm;          % Final wheel speed [RPM]

    %% 4. Generate Time Vector & RWA Disturbance Profiles
    t = 0 : (1/fs) : sim_time;

    [Fx, Fy, M_theta] = fem.generate_rwa_disturbance(...
        t, rpm_start, rpm_end, unbalance, rotor_inertia);

    %% 5. Assemble Simulation Options Struct
    struct_out   = struct();
    struct_out.t = t;

    idx = 1;
    for k = 1:num_rwas
        struct_out.forces(idx) = struct(...
            'node_id', rwa_ids(k), 'dof', "x", 'profile', Fx);
        struct_out.forces(idx + 1) = struct(...
            'node_id', rwa_ids(k), 'dof', "y", 'profile', Fy);
        % Uncomment below to also apply the reaction torque to the model:
        % struct_out.forces(idx + 2) = struct(...
        %     'node_id', rwa_ids(k), 'dof', "theta", 'profile', M_theta);
        idx = idx + 2;
    end

    %% 6. Package Diagnostics (Optional Second Output)
    srp_info = struct(...
        'M_total',      M_total, ...
        'I_total',      I_total, ...
        'x_cg',         x_cg, ...
        'x_centerline', x_centerline, ...
        'width_bottom', width_bottom, ...
        'A_bottom',     A_bottom, ...
        'lever_arm',    lever_arm, ...
        'F_srp',        F_srp, ...
        'T_srp',        T_srp, ...
        'alpha_sv_uncontrolled', T_srp / I_total, ...
        'num_rwas',     num_rwas, ...
        'T_wheel',      T_wheel, ...
        'alpha_wheel',  alpha_wheel, ...
        'rpm_start',    rpm_start, ...
        'rpm_end',      rpm_end, ...
        'delta_rpm',    delta_rpm);
end
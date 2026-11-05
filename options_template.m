function struct_out = options_template()
    %% =====================================================================
    %  DYNAMIC SIMULATION OPTIONS DEFINITIONS
    %  =====================================================================
    %  The options struct configures the state-space simulation parameters
    %  and external forcing functions.
    %
    %  Fields:
    %  - t      : [1xN Double] Time vector in seconds (s). Determines simulation
    %             length and resolution.
    %  - z      : [Double] Proportional damping ratio (zeta). e.g., 0.01 = 1%.
    %  - forces : [1xM Struct Array] Defines external perturbing forces.
    %      Sub-fields for each force:
    %      - node_id : [String] Must perfectly match an existing node's 'id'.
    %      - dof     : [String] Direction of application.
    %                  Available Options: "x", "y", "theta".
    %      - profile : [1xN Double] Time-history vector of force/moment magnitude.
    %                  Length MUST equal the length of 't'.
    %                  Units: Newtons (N) for x/y, Newton-meters (N*m) for theta.
    
    struct_out = struct();
    
    % 1. Define time parameters and generate the time vector
    fs = 1000; % Sampling frequency in Hz
    t_span = [0, 20]; % Start and end time in seconds
    t = t_span(1) : 1/fs : t_span(2);
    
    % Assign base simulation parameters
    struct_out.t = t;
    struct_out.z = 0.01; % 1% structural damping
    
    % 2. Define external forcing profiles
    % Create a 2 Hz harmonic profile with 0.1 amplitude (e.g., Newtons)
    harmonic_profile = 0.1 * cos(2 * pi * 2 * t);
    
    % 3. Apply forces to specific structural Degrees of Freedom
    % Apply the harmonic profile to the X translation of the "Top" node
    struct_out.forces(1) = struct(...
        'node_id', "Top", ...
        'dof', "x", ...
        'profile', harmonic_profile);
    
    % Apply the same harmonic profile to the Y translation of the "Top" node
    % (Note: Multiple forces can be applied to the same node, or even the
    % same node & DOF to superimpose their effects).
    struct_out.forces(2) = struct(...
        'node_id', "Top", ...
        'dof', "y", ...
        'profile', harmonic_profile);
    
    end
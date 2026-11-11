function struct_out = options_sm_motor()
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

% Define time parameters and generate the time vector
fs = 1000; % Sampling frequency in Hz
t_span = [0, 60]; % Start and end time in seconds
t = t_span(1) : 1/fs : t_span(2);

% Assign base simulation parameters
struct_out.t = t;
struct_out.z = 0.01; % 1% structural damping

% Define external forcing profile
dist = 0.1 * sin(2 * pi * t);

% Disturbances: Left RWA
struct_out.forces(1) = struct(...
    'node_id', "SecondaryMirror", ...
    'dof', "theta", ...
    'profile', dist);
end
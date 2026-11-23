function [struct_out, options_out] = config_fem_template()
    % CONFIG_FEM_TEMPLATE Generates a well-documented FEM structural template
    % along with a simulation options template for dynamic response analysis.
    %
    % This function outputs:
    %   1. struct_out  : Formatted for the 'fem' class matrix generation methods.
    %   2. options_out : Formatted for the 'fem.calculate_rms' and 'fem.plot_rms' methods.
    
    struct_out = struct();
    
    %% =====================================================================
    %  GLOBAL STRUCTURAL PROPERTIES
    %  =====================================================================
    %  - z : [Double] Global damping ratio as a decimal percentage
    %        (e.g., 0.01 = 1% damping, 0.25 = 25% damping).
    %        Applied to all elements unless overridden by an element's c_damping field.
    struct_out.z = 0.01;
    
    %% =====================================================================
    %  NODE DEFINITIONS
    %  =====================================================================
    %  Each node requires the following fields:
    %  - id          : [String] Unique identifier for the node.
    %  - x, y        : [Double] Spatial coordinates in meters (m).
    %  - mass        : [Double] Point translational mass in kg.
    %  - inertia     : [Double] Point rotational mass moment of inertia in kgm^2.
    %  - constraints : [1x3 Array] Boundary conditions for [x, y, theta].
    %                  0 = Free to move.
    %                  1 = Fixed / Grounded to the "universe".
    %  - k_rot       : [Double] Rotational stiffness in Nm/rad. Used to ground
    %                  rotational DOFs (e.g., modeling an ACS or preventing
    %                  singular matrices on pin joints).
    
    % Node 1: Fully Fixed Base
    struct_out.nodes(1) = struct(...
    'id', "Bottom Left", ...
    'x', -1, 'y', 0, ...
    'mass', 10, ...
    'inertia', 20, ...
    'constraints', [1, 1, 1], ...
    'k_rot', 0);
    
    % Node 2: Free space node with an attached rotational spring
    struct_out.nodes(2) = struct(...
    'id', "Top", ...
    'x', 0,  'y', 2, ...
    'mass', 100, ...
    'inertia', 200, ...
    'constraints', [0, 0, 0], ...
    'k_rot', 10);
    
    % Node 3: Pinned connection
    struct_out.nodes(3) = struct(...
    'id', "Bottom Right", ...
    'x', 1,  'y', 0, ...
    'mass', 10, ...
    'inertia', 20, ...
    'constraints', [1, 1, 0], ...
    'k_rot', 0);
    
    % Node 4: Spring attached payload mass
    struct_out.nodes(4) = struct(...
    'id', "Payload", ...
    'x', 0,  'y', 3, ...
    'mass', 50, ...
    'inertia', 50, ...
    'constraints', [0, 0, 0], ...
    'k_rot', 0);
    
    %% =====================================================================
    %  ELEMENT DEFINITIONS (CELL ARRAY OF STRUCTS)
    %  =====================================================================
    %  Use curly braces {1}, {2}, ... so each element only includes the fields
    %  relevant to its type (no unused fields or dissimilar-struct errors).
    %
    %  Required for ALL elements:
    %  - id          : [String] Unique identifier for the element.
    %  - type        : [String] Element formulation type: "beam" or "spring".
    %  - nodes       : [1x2 String Array] Exact IDs of the connecting nodes.
    %
    %  Required for "beam" elements ONLY:
    %  - E           : [Double] Young's Modulus in Pascals (Pa).
    %  - I           : [Double] Area Moment of Inertia in m^4.
    %  - A           : [Double] Cross-sectional area in m^2.
    %  - rho         : [Double] Mass density per unit length in kg/m.
    %
    %  Required for "spring" elements ONLY:
    %  - k_spring    : [1x3 Array] Local stiffness [k_axial, k_shear, k_torsion]
    %                  in [N/m, N/m, N*m/rad].
    %
    %  Optional for ANY element:
    %  - c_damping   : [1x3 Array] Local DOF damping ratios [x_damp, y_damp, theta_damp]
    %                  expressed as decimal percentages (e.g., 0.25 = 25% damping).
    %                  Overrides global struct_out.z for this element.
    %                  Omit this field entirely to use global struct_out.z;
    %                  individual NaN entries also fall back to struct_out.z.
    
    % Element 1: Left Beam (only beam fields; uses global damping z)
    struct_out.elements{1} = struct(...
    'id', "Beam Left", ...
    'type', "beam", ...
    'nodes', ["Bottom Left", "Top"], ...
    'E', 70e9, 'I', 1e-4, 'A', 0.01, 'rho', 1);
    
    % Element 2: Right Beam (only beam fields; uses global damping z)
    struct_out.elements{2} = struct(...
    'id', "Beam Right", ...
    'type', "beam", ...
    'nodes', ["Top", "Bottom Right"], ...
    'E', 70e9, 'I', 1e-4, 'A', 0.01, 'rho', 2);
    
    % Element 3: Base Beam (with optional per-DOF damping override)
    struct_out.elements{3} = struct(...
    'id', "Beam Base", ...
    'type', "beam", ...
    'nodes', ["Bottom Right", "Bottom Left"], ...
    'E', 70e9, 'I', 1e-4, 'A', 0.01, 'rho', 3, ...
    'c_damping', [0.05, NaN, NaN]); % 5% axial damping; y & theta fall back to global z
    
    % Element 4: Discrete Isolator Spring (only spring fields + damping override)
    struct_out.elements{4} = struct(...
    'id', "Payload Isolator", ...
    'type', "spring", ...
    'nodes', ["Top", "Payload"], ...
    'k_spring', [1e5, 1e4, 5e3], ...      % [Axial, Shear, Torsion] stiffness
    'c_damping', [0.25, 0.20, 0.10]);     % [x_damp, y_damp, theta_damp] override
    
    %% =====================================================================
    %  DYNAMIC SIMULATION OPTIONS DEFINITIONS
    %  =====================================================================
    %  Note: Damping (z) is defined on struct_out rather than options_out.
    
    options_out = struct();
    
    fs = 1000;
    t_span = [0, 20];
    t = t_span(1) : 1/fs : t_span(2);
    
    options_out.t = t;
    
    harmonic_profile = 0.1 * cos(2 * pi * 2 * t);
    
    options_out.forces(1) = struct('node_id', "Top", 'dof', "x", 'profile', harmonic_profile);
    options_out.forces(2) = struct('node_id', "Top", 'dof', "y", 'profile', harmonic_profile);
    
    end
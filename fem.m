classdef fem
    methods (Static)
        function el = get_element(elements, idx)
            % GET_ELEMENT Helper to retrieve an element struct from either a
            % cell array of structs or a legacy struct array.
            if iscell(elements)
                el = elements{idx};
            else
                el = elements(idx);
            end
        end

        % ===============================================================================
        
        function model = generate_fem_matrices(model)
        % GENERATE_FEM_MATRICES Assembles global M, K, and C matrices, mass properties,
        % natural frequencies, mode shapes, and modal participation factors.
        %
        % SYNTAX:
        %   model = fem.generate_fem_matrices(model)
        %
        % INPUTS:
        %   model - Struct containing structural definitions:
        %     .nodes (1xN struct array) with fields:
        %       - id          : String unique node identifier
        %       - x, y        : Node coordinates [m]
        %       - mass        : Point translational mass [kg]
        %       - inertia     : Point rotational mass moment of inertia [kg*m^2]
        %       - constraints : 1x3 binary array [x, y, theta] (1 = fixed, 0 = free)
        %       - k_rot       : Point rotational stiffness [N*m/rad]
        %     .elements (1xM struct array) with fields:
        %       - id          : String unique element identifier
        %       - type        : String ("beam", etc.)
        %       - nodes       : 1x2 string array of connecting node IDs
        %       - E           : Young's Modulus [Pa]
        %       - I           : Area Moment of Inertia [m^4]
        %       - A           : Cross-sectional Area [m^2]
        %       - rho         : Mass density per unit length [kg/m]
        %       - k_discrete  : Discrete spring stiffness (if non-beam) [N/m]
        %       - c_damping   : Optional 1x3 array [x_damp, y_damp, theta_damp]
        %                       in local element coordinates that overrides model.z.
        %                       Use [] or omit to use global z; NaN entries also
        %                       fall back to global z.
        %     .z - Global damping ratio (e.g., 0.01 for 1%, 0.25 for 25%).
        %     
        %
        % OUTPUTS:
        %   model - Updated input struct with appended fields:
        %     .K       : Global Stiffness Matrix (3N x 3N) [N/m, N*m/rad]
        %     .M       : Global Consistent Mass Matrix (3N x 3N) [kg, kg*m^2]
        %     .C       : Global Damping Matrix (3N x 3N) [N/(m/s), N*m/(rad/s)]
        %     .mapping : 1xN struct array linking node IDs to global DOF indices
        %                Fields: .node_id, .global_dofs (1x3 array: [x, y, theta])
        %     .fn      : Natural Frequencies of all system modes (3Nx1)
        %     .modes   : All system mode shapes (3Nx3N)
        %
        % DEPENDENCIES:
        %   None (Built-in MATLAB functions only)

        % Extract arrays for easier referencing
        nodes = model.nodes;
        elements = model.elements;
        
        num_nodes = length(nodes);
        num_elements = length(elements);
        num_dofs = num_nodes * 3; % 3 DOFs per node in 2D: [x, y, theta]
        
        if isfield(model, 'z') && ~isempty(model.z)
            z_global = model.z;
        else
            z_global = 0;
        end
        
        % 1. Create DOF Mapping and Identify Constrained / Free DOFs
        node_ids = [nodes.id];
        dof_map = struct('node_id', {}, 'global_dofs', {});
        fixed_mask = false(num_dofs, 1);
        
        for i = 1:num_nodes
            dofs = [3*i-2, 3*i-1, 3*i];
            dof_map(i).node_id = nodes(i).id;
            dof_map(i).global_dofs = dofs;
            
            if isfield(nodes(i), 'constraints') && ~isempty(nodes(i).constraints)
                if length(nodes(i).constraints) ~= 3
                    error('Node "%s" constraints must be a 1x3 vector [x, y, theta].', nodes(i).id);
                end
                fixed_mask(dofs) = logical(nodes(i).constraints(:));
            end
        end
        
        free_dofs = find(~fixed_mask);
        fixed_dofs = find(fixed_mask);
        
        if isempty(free_dofs)
            error('All degrees of freedom are constrained. At least one free DOF is required.');
        end
        
        model.mapping = dof_map;
        
        % 2. Initialize Global Matrices and Element Caches (for damping assembly)
        K_global = zeros(num_dofs, num_dofs);
        M_global = zeros(num_dofs, num_dofs);
        
        el_k_local = cell(1, num_elements);
        el_T = cell(1, num_elements);
        el_dofs = cell(1, num_elements);
        
        % 3. Assemble Element Matrices
        for i = 1:num_elements
            el = fem.get_element(elements, i);
            
            % Find indices for connecting nodes
            n1_idx = find(node_ids == el.nodes(1));
            n2_idx = find(node_ids == el.nodes(2));
    
            n1 = nodes(n1_idx);
            n2 = nodes(n2_idx);
    
            % Element Geometry & Transformation
            dx = n2.x - n1.x;
            dy = n2.y - n1.y;
            L = sqrt(dx^2 + dy^2);
            
            if L > 1e-12
                c = dx / L;
                s = dy / L;
            else
                c = 1; s = 0; % Fallback for coincident nodes if they exist
            end
            
            T = [ c,  s,  0,  0,  0,  0;
                 -s,  c,  0,  0,  0,  0;
                  0,  0,  1,  0,  0,  0;
                  0,  0,  0,  c,  s,  0;
                  0,  0,  0, -s,  c,  0;
                  0,  0,  0,  0,  0,  1 ];
        
            if el.type == "beam"
                % Local Stiffness Matrix
                E = el.E; I = el.I; A = el.A;
                k_local = [ E*A/L, 0,          0,         -E*A/L, 0,          0;
                    0,     12*E*I/L^3, 6*E*I/L^2, 0,     -12*E*I/L^3, 6*E*I/L^2;
                    0,     6*E*I/L^2,  4*E*I/L,   0,     -6*E*I/L^2,  2*E*I/L;
                    -E*A/L, 0,          0,          E*A/L, 0,          0;
                    0,    -12*E*I/L^3,-6*E*I/L^2, 0,      12*E*I/L^3,-6*E*I/L^2;
                    0,     6*E*I/L^2,  2*E*I/L,   0,     -6*E*I/L^2,  4*E*I/L ];
        
                % Local Consistent Mass Matrix
                rho = el.rho;
                m_local = (rho*L / 420) * ...
                    [ 140, 0,      0,       70,  0,      0;
                    0,   156,    22*L,    0,   54,    -13*L;
                    0,   22*L,   4*L^2,   0,   13*L,  -3*L^2;
                    70,  0,      0,       140, 0,      0;
                    0,   54,     13*L,    0,   156,   -22*L;
                    0,  -13*L,  -3*L^2,   0,  -22*L,   4*L^2 ];
                    
            elseif el.type == "spring"
                % Local Stiffness Matrix (Uncoupled Axial, Shear, Torsion)
                ka = el.k_spring(1);
                ks = el.k_spring(2);
                kt = el.k_spring(3);
                
                k_local = [ ka,  0,  0, -ka,  0,  0;
                             0, ks,  0,   0,-ks,  0;
                             0,  0, kt,   0,  0,-kt;
                           -ka,  0,  0,  ka,  0,  0;
                             0,-ks,  0,   0, ks,  0;
                             0,  0,-kt,   0,  0, kt ];
                
                % Springs are considered massless elements
                m_local = zeros(6,6);
            else
                error('Unknown element type "%s" for element "%s". Use "beam" or "spring".', el.type, el.id);
            end
    
            % Transform to Global Frame
            K_e = T' * k_local * T;
            M_e = T' * m_local * T;
    
            % Identify Global DOFs for Assembly
            dofs = [dof_map(n1_idx).global_dofs, dof_map(n2_idx).global_dofs];
            
            % Cache element stiffness, transform, and DOFs for damping matrix assembly
            el_k_local{i} = k_local;
            el_T{i} = T;
            el_dofs{i} = dofs;
    
            % Add to Global Matrices
            K_global(dofs, dofs) = K_global(dofs, dofs) + K_e;
            M_global(dofs, dofs) = M_global(dofs, dofs) + M_e;
        end
        
        % 4. Add Nodal Point Masses, Inertias, and Stiffness
        for i = 1:num_nodes
            dofs = dof_map(i).global_dofs;
        
            % Translational Mass
            M_global(dofs(1), dofs(1)) = M_global(dofs(1), dofs(1)) + nodes(i).mass;
            M_global(dofs(2), dofs(2)) = M_global(dofs(2), dofs(2)) + nodes(i).mass;
        
            % Rotational Inertia and Stiffness
            M_global(dofs(3), dofs(3)) = M_global(dofs(3), dofs(3)) + nodes(i).inertia;
            K_global(dofs(3), dofs(3)) = K_global(dofs(3), dofs(3)) + nodes(i).k_rot;
        end
        
        % Enforce exact numerical symmetry on full matrices
        K_global = 0.5 * (K_global + K_global');
        M_global = 0.5 * (M_global + M_global');
        
        % 5. Build Rigid Body Transformation Matrix and Mass Properties
        phi_matrix = zeros(num_dofs, 3);
        for i = 1:num_nodes
            dofs = dof_map(i).global_dofs;
            x = nodes(i).x;
            y = nodes(i).y;
            
            % Rigid body displacement mappings for [X_trans, Y_trans, Z_rot]
            phi_matrix(dofs(1), :) = [1, 0, -y];
            phi_matrix(dofs(2), :) = [0, 1,  x];
            phi_matrix(dofs(3), :) = [0, 0,  1];
        end
        
        % Compute Rigid Body Mass Matrix to extract properties
        M_rb = phi_matrix' * M_global * phi_matrix; 
        
        M_total = M_rb(1, 1); 
        if M_total > 1e-12
            X_cg = M_rb(2, 3) / M_total;
            Y_cg = -M_rb(1, 3) / M_total;
        else
            X_cg = 0;
            Y_cg = 0;
        end
        
        % Total Inertia about the CG using the Parallel Axis Theorem
        I_total = M_rb(3, 3) - M_total * (X_cg^2 + Y_cg^2);
        
        % 6. Calculate Natural Frequencies and Mode Shapes (Reduced to Free DOFs)
        K_ff = K_global(free_dofs, free_dofs);
        M_ff = M_global(free_dofs, free_dofs);
        
        [V_f, D_f] = eig(K_ff, M_ff);
        
        % Explicitly Mass-Normalize the eigenvectors on the free DOFs
        m_modal = real(diag(V_f' * M_ff * V_f));
        V_f_norm = V_f * diag(1 ./ sqrt(m_modal));
    
        evals = real(diag(D_f));
        evals(abs(evals) < 1e-10) = 0; 
    
        omega = sqrt(abs(evals));
        fn_all = omega / (2 * pi);
    
        [fn_sorted, sort_idx] = sort(fn_all);
        V_f_sorted = V_f_norm(:, sort_idx);
        
        % Expand mode shapes back to full DOF space (0 at fixed DOFs)
        modes_sorted = zeros(num_dofs, length(free_dofs));
        modes_sorted(free_dofs, :) = V_f_sorted;
        
        % 7. Assemble Global Damping Matrix (C)
        %    Baseline: Exact modal damping (z_global across all flexible modes on free DOFs)
        %    Overrides: Incremental local stiffness-proportional damping for elements
        %               that explicitly define c_damping
        omega_sorted = fn_sorted * 2 * pi;
        omega_modal_damp = omega_sorted;
        omega_modal_damp(fn_sorted <= 1e-4) = 0; % Keep rigid-body modes undamped
        
        % C_ff_base = M_ff * V_f * diag(2 * z_global * omega_i) * V_f' * M_ff
        M_f_Phi = M_ff * V_f_sorted;
        C_global = zeros(num_dofs, num_dofs);
        C_global(free_dofs, free_dofs) = M_f_Phi * diag(2 * z_global * omega_modal_damp) * M_f_Phi';
        
        % Reference frequency for incremental element damping overrides
        if isfield(model, 'f_ref') && ~isempty(model.f_ref) && model.f_ref > 0
            omega_ref = model.f_ref * 2 * pi;
        else
            flex_fn = fn_sorted(fn_sorted > 1e-4);
            if ~isempty(flex_fn)
                omega_ref = flex_fn(1) * 2 * pi;
            else
                omega_ref = max(omega_sorted);
            end
        end
        
        % Add incremental damping only for elements with a c_damping override
        if omega_ref > 1e-12
            for i = 1:num_elements
                el = fem.get_element(elements, i);
                
                if isfield(el, 'c_damping') && ~isempty(el.c_damping)
                    z_el = el.c_damping;
                    if length(z_el) ~= 3
                        error('Element "%s" c_damping must be a 1x3 vector [x_damp, y_damp, theta_damp].', el.id);
                    end
                    % Fall back to global z for any NaN entries
                    z_el(isnan(z_el)) = z_global;
                    
                    % Incremental damping relative to the modal baseline
                    delta_z = z_el - z_global;
                    if any(abs(delta_z) > 1e-12)
                        beta_delta = (2 .* delta_z) ./ omega_ref;
                        B_mat = diag([beta_delta, beta_delta]);
                        
                        % Symmetric real local incremental damping (handles +/- delta_z)
                        delta_c_local = 0.5 * (B_mat * el_k_local{i} + el_k_local{i} * B_mat);
                        delta_C_e = el_T{i}' * delta_c_local * el_T{i};
                        
                        dofs = el_dofs{i};
                        C_global(dofs, dofs) = C_global(dofs, dofs) + delta_C_e;
                    end
                end
            end
        end
        
        % Enforce exact numerical symmetry
        C_global = 0.5 * (C_global + C_global');
        
        % 8. Calculate Modal Participation and Effective Mass
        % L = Phi_norm^T * M * Phi_rb (fixed DOFs have 0 in modes_sorted)
        participation_factors = modes_sorted' * M_global * phi_matrix; 
        effective_mass = participation_factors.^2;
        
        % Percentages (Handling zero-mass directions safely)
        M_rb_diag = diag(M_rb)';
        M_rb_diag(M_rb_diag < 1e-12) = 1; 
        effective_mass_pct = (effective_mass ./ M_rb_diag) * 100;
    
        % 9. Pack outputs back into the struct
        model.z = z_global;
        model.free_dofs = free_dofs;
        model.fixed_dofs = fixed_dofs;
        model.K = K_global;
        model.M = M_global;
        model.C = C_global;
        model.phi_matrix = phi_matrix;
        model.M_total = M_total;
        model.I_total = I_total;
        model.CG = [X_cg, Y_cg];
        model.fn = fn_sorted;
        model.modes = modes_sorted;
        model.participation_factors = participation_factors;
        model.effective_mass = effective_mass;
        model.effective_mass_pct = effective_mass_pct;
        end
    
    % ===============================================================================

        function [X_rms, X] = calculate_rms(model, options)
            % CALCULATE_PLOT_RMS Simulates transient dynamic response under external
            % perturbing forces using state-space representation and plots RMS metrics.
            %
            % SYNTAX:
            %   [X_rms, X] = fem.calculate_plot_rms(model, options)
            %
            % INPUTS:
            %   model   - Struct containing assembled system matrices and mapping:
            %             .M       : Mass matrix (n_dof x n_dof)
            %             .K       : Stiffness matrix (n_dof x n_dof)
            %             .C       : Damping matrix (n_dof x n_dof)
            %             .mapping : DOF mapping struct array
            %   options - Dynamic simulation setup struct:
            %             .t      : Time vector [s]
            %             .forces : Struct array of forces with fields:
            %                       - node_id : String node ID target
            %                       - dof     : Direction string ("x", "y", "theta")
            %                       - profile : Time history force vector (1 x length(t))
            %
            % OUTPUTS:
            %   X_rms - RMS displacement vector for all global DOFs (n_dof x 1)
            %   X     - Time-history displacement matrix (n_dof x length(t))
            %
            % DEPENDENCIES:
            %   Control System Toolbox (ss, lsim)
            
            % Extract matrices and parameters
            M = model.M;
            K = model.K;
            C = model.C;
            mapping = model.mapping;

            t = options.t;
            n_dof = size(M, 1);

            if isfield(model, 'free_dofs') && ~isempty(model.free_dofs)
                free_dofs = model.free_dofs;
            else
                free_dofs = (1:n_dof)';
            end

            % Build Full Force Matrix
            F_matrix = zeros(n_dof, length(t));
            if isfield(options, 'forces')
                for i = 1:length(options.forces)
                    f = options.forces(i);

                    node_idx = find([mapping.node_id] == f.node_id);
                    if isempty(node_idx)
                        error('Node ID %s not found in mapping.', f.node_id);
                    end

                    node_dofs = mapping(node_idx).global_dofs;

                    if strcmpi(f.dof, 'x')
                        g_dof = node_dofs(1);
                    elseif strcmpi(f.dof, 'y')
                        g_dof = node_dofs(2);
                    elseif strcmpi(f.dof, 'theta')
                        g_dof = node_dofs(3);
                    else
                        error('Invalid DOF type. Use ''x'', ''y'', or ''theta''.');
                    end

                    F_matrix(g_dof, :) = F_matrix(g_dof, :) + f.profile;
                end
            end

            % Reduce system to free DOFs (enforcing nodal constraints)
            M_f = M(free_dofs, free_dofs);
            K_f = K(free_dofs, free_dofs);
            C_f = C(free_dofs, free_dofs);
            F_f = F_matrix(free_dofs, :);
            n_f = length(free_dofs);

            % Convert to State-Space Formulation using pre-assembled model.C
            A = [zeros(n_f), eye(n_f);
                -M_f \ K_f, -M_f \ C_f];
            B = [zeros(n_f);
                M_f \ eye(n_f)];
            C_ss = [eye(n_f), zeros(n_f)]; 
            D_ss = zeros(n_f);

            sys = ss(A, B, C_ss, D_ss);

            % Simulate response on free DOFs
            X_f = lsim(sys, F_f', t)'; 

            % Expand response back to full DOF space (0 on fixed DOFs)
            X = zeros(n_dof, length(t));
            X(free_dofs, :) = X_f;

            % Calculate RMS response
            X_rms = sqrt(mean(X.^2, 2));
        end

        % ===============================================================================

        function combined_model = merge_fem_models(varargin)
            % MERGE_FEM_MODELS Combines multiple sub-structure models based on spatial
            % node coordinates within a 1e-6 tolerance.
            %
            % SYNTAX:
            %   combined_model = fem.merge_fem_models(model_1, model_2, ...)
            %
            % INPUTS:
            %   varargin - Variable number of model structs to combine.
            %
            % BEHAVIOR & CONSTRAINTS:
            %   - Coincident nodes within 1e-6 m are merged.
            %   - Node IDs are joined with "-" if differing; duplicates retained once.
            %   - Node properties default to primary; warnings issued on mismatch.
            %   - Constraints combine using logical OR (most restrictive state).
            %   - Element IDs prepend parent variable name prefix formatted as "Name - ID".
            %
            % OUTPUTS:
            %   combined_model - Merged struct containing consolidated nodes, elements,
            %                    and re-generated global K, M, and mapping matrices.
            %
            % DEPENDENCIES:
            %   fem.generate_fem_matrices
            
            tol = 1e-6;
            merged_nodes = [];
            merged_elements = {}; % Cell array allows heterogeneous element structs
            merged_z = [];
            
            map_old_id = cell(1, nargin);
            map_new_idx = cell(1, nargin);
            
            for i = 1:nargin
                mdl = varargin{i};
                
                % Reconcile global damping factor z across models
                if isfield(mdl, 'z') && ~isempty(mdl.z)
                    if isempty(merged_z)
                        merged_z = mdl.z;
                    elseif abs(merged_z - mdl.z) > tol
                        warning('Global damping z mismatch between models. Keeping primary z: %f', merged_z);
                    end
                end
                
                map_old_id{i} = strings(1, length(mdl.nodes));
                map_new_idx{i} = zeros(1, length(mdl.nodes));
                
                %% 1. Process and Merge Nodes
                for j = 1:length(mdl.nodes)
                    cn = mdl.nodes(j);
                    match_found = false;
                    
                    for k = 1:length(merged_nodes)
                        dx = merged_nodes(k).x - cn.x;
                        dy = merged_nodes(k).y - cn.y;
                        
                        if sqrt(dx^2 + dy^2) <= tol
                            match_found = true;
                            
                            if abs(merged_nodes(k).mass - cn.mass) > tol
                                warning('Mass mismatch at merged node %s. Keeping primary mass: %f', merged_nodes(k).id, merged_nodes(k).mass);
                            end
                            if abs(merged_nodes(k).inertia - cn.inertia) > tol
                                warning('Inertia mismatch at merged node %s. Keeping primary inertia: %f', merged_nodes(k).id, merged_nodes(k).inertia);
                            end
                            
                            merged_nodes(k).constraints = merged_nodes(k).constraints | cn.constraints;
                            
                            if merged_nodes(k).id ~= cn.id
                                if ~contains(merged_nodes(k).id, cn.id)
                                    merged_nodes(k).id = merged_nodes(k).id + "-" + cn.id;
                                end
                            end
                            
                            map_old_id{i}(j) = cn.id;
                            map_new_idx{i}(j) = k;
                            break;
                        end
                    end
                    
                    if ~match_found
                        if isempty(merged_nodes)
                            merged_nodes = cn;
                        else
                            merged_nodes(end+1) = cn;
                        end
                        map_old_id{i}(j) = cn.id;
                        map_new_idx{i}(j) = length(merged_nodes);
                    end
                end
            end
            
            %% 2. Process Elements and Remap to Merged Nodes
            for i = 1:nargin
                mdl = varargin{i};
                
                parent_name = inputname(i);
                if isempty(parent_name)
                    parent_name = "Model " + i; 
                end
                parent_name = strrep(string(parent_name), "_", " ");
                
                for j = 1:length(mdl.elements)
                    ce = fem.get_element(mdl.elements, j);
                    
                    ce.id = parent_name + " - " + ce.id;
                    
                    n1_old_idx = find(map_old_id{i} == ce.nodes(1));
                    new_n1_id = merged_nodes(map_new_idx{i}(n1_old_idx)).id;
                    
                    n2_old_idx = find(map_old_id{i} == ce.nodes(2));
                    new_n2_id = merged_nodes(map_new_idx{i}(n2_old_idx)).id;
                    
                    ce.nodes = [new_n1_id, new_n2_id];
                    
                    merged_elements{end+1} = ce; %#ok<AGROW>
                end
            end
            
            %% 3. Package and Generate Global Matrices
            combined_model.nodes = merged_nodes;
            combined_model.elements = merged_elements;
            if ~isempty(merged_z)
                combined_model.z = merged_z;
            end
            
            combined_model = fem.generate_fem_matrices(combined_model);
        end
        
        %% =============================================================================== 
        function [Fx, Fy, M_theta] = generate_rwa_disturbance(t, rpm_start, rpm_end, unbalance, rotor_inertia)
            % GENERATE_RWA_DISTURBANCE Generates the time-domain force and torque profiles
            % for the Reaction Wheel Assembly (RWA) sweeps across all three DOFs.
            %
            % SYNTAX:
            %   [Fx, Fy, M_theta] = generate_rwa_disturbance(t, rpm_start, rpm_end, unbalance, rotor_inertia)
            %
            % INPUTS:
            %   t             - Linearly spaced time vector [s] (e.g., 0:0.01:3600)
            %   rpm_start     - Starting angular velocity [RPM] (e.g., -5000)
            %   rpm_end       - Ending angular velocity [RPM] (e.g., 5000)
            %   unbalance     - Mass unbalance property of the RWA [kgm]
            %   rotor_inertia - Rotational inertia of the wheel's rotor [kgm^2]
            %
            % OUTPUTS:
            %   Fx      - Force profile in the x-direction [N] (1 x length(t))
            %   Fy      - Force profile in the y-direction [N] (1 x length(t))
            %   M_theta - Torque profile in the theta-direction [N*m] (1 x length(t))
            %
            % DEPENDENCIES:
            %   None

            % Normalize time vector to start at 0 for calculations
            t_sweep = t - t(1);
            t_total = t_sweep(end);

            % Convert RPM bounds to angular velocity (rad/s)
            omega_start = rpm_start * (pi / 30);
            omega_end = rpm_end * (pi / 30);

            % Calculate constant angular acceleration (alpha = dw/dt)
            alpha = (omega_end - omega_start) / t_total;

            % Instantaneous angular velocity over time
            omega = omega_start + alpha .* t_sweep;

            % Cumulative rotational phase angle via integration of omega
            theta_rot = omega_start .* t_sweep + 0.5 * alpha .* (t_sweep.^2);

            % Compute magnitude of the centrifugal force (F = m*r*omega^2)
            force_magnitude = unbalance .* (omega.^2);

            % Resolve force into x and y components
            Fx = force_magnitude .* cos(theta_rot);
            Fy = force_magnitude .* sin(theta_rot);

            % Compute the reaction torque on the structure (T_reaction = - I * alpha)
            M_theta = (rotor_inertia * alpha) .* ones(size(t));


        end

        %% =============================================================================== 
        function [error_total, components] = calculate_optical_pointing_error(model, X, optical_nodes)
            % CALCULATE_OPTICAL_POINTING_ERROR Converts raw global DOF time-history
            % response into the optical pointing error metric (radians).
            %
            % SYNTAX:
            %   [pointing_error, components] = calculate_optical_pointing_error(model, X, optical_nodes)
            %
            % INPUTS:
            %   model         - Struct containing assembled system matrices and mapping:
            %                   .mapping : DOF mapping struct array linking node IDs
            %                              to global matrix indices [x, y, theta].
            %   X             - Time-history displacement matrix (n_dof x length(t))
            %                   outputted by the dynamic simulation.
            %   optical_nodes - Struct defining the relevant optical node IDs:
            %                   .sm_id  : String ID of the Secondary Mirror node
            %                   .pm_ids : 1D Cell Array of String IDs for Primary Mirror nodes
            %
            % OUTPUTS:
            %   pointing_error - Maximum total pointing error over time across all PM paths (1 x length(t))
            %   components     - Struct containing the individual pointing errors:
            %                    .sm    : Error contribution from the SM (1 x length(t))
            %                    .pms   : Error contributions from each PM (num_pms x length(t))
            %                    .total : Total error for each PM/SM path (num_pms x length(t))
            %
            % DEPENDENCIES:
            %   None

            % Extract the mapping dictionary
            mapping = model.mapping;
            nodes = model.nodes;
            all_map_ids = string({mapping.node_id});
            all_node_ids = string({nodes.id});

            function dofs = get_node_dofs(target_id)
                target_str = string(target_id);
                is_match = strcmp(all_map_ids, target_str);

                if ~any(is_match)
                    error('OpticalNode:NotFound', ...
                        'The optical node "%s" was not found in the model mapping. Please check for typos.', target_str);
                end

                matched_node = mapping(is_match);
                dofs = matched_node(1).global_dofs;
            end

            function x_pos = get_node_x(target_id)
                target_str = string(target_id);
                is_match = strcmp(all_node_ids, target_str);

                if ~any(is_match)
                    error('OpticalNode:NotFound', ...
                        'The optical node "%s" was not found in model.nodes. Please check for typos.', target_str);
                end

                matched_node = nodes(is_match);
                x_pos = matched_node(1).x;
            end

            % Secondary mirror pointing error (from rotational DOF)
            sm_dofs = get_node_dofs(optical_nodes.sm_id);
            theta_sm = X(sm_dofs(3), :);
            error_sm = 17 * theta_sm;

            % Primary mirror pointing error (linear least-squares slope of y-displacements)
            num_pms = length(optical_nodes.pm_ids);
            num_time_steps = size(X, 2);

            pm_x = zeros(num_pms, 1);
            y_pms = zeros(num_pms, num_time_steps);
            theta_pms = zeros(num_pms, num_time_steps);

            for i = 1:num_pms
                pm_id = optical_nodes.pm_ids{i};
                pm_dofs = get_node_dofs(pm_id);
                pm_x(i) = get_node_x(pm_id);
                y_pms(i, :) = X(pm_dofs(2), :);
                theta_pms(i, :) = X(pm_dofs(3), :);
            end

            x_centered = pm_x - mean(pm_x);
            denom = sum(x_centered.^2);

            if num_pms >= 2 && denom > 1e-12
                y_centered = y_pms - mean(y_pms, 1);
                pm_slope = (x_centered' * y_centered) / denom;
                theta_pm = atan(pm_slope);
            else
                warning('OpticalNode:SingleOrCollinearPM', ...
                    'Primary mirror requires at least 2 nodes with distinct x-coordinates for linear slope fitting. Falling back to nodal theta.');
                theta_pm = mean(theta_pms, 1);
            end

            error_pm = 2 * theta_pm;
            error_total_rad = abs(error_pm) + abs(error_sm);

            % Convert output to arc-seconds (1 rad = (180/pi)*3600 arc-sec)
            rad_to_arcsec = (180 / pi) * 3600;

            components.sm = error_sm * rad_to_arcsec;
            components.pm = error_pm * rad_to_arcsec;
            components.pms = components.pm; % Retained for backwards compatibility
            components.total = error_total_rad * rad_to_arcsec;

            error_total = components.total;
        end

        %% =============================================================================== 
        function [A, B, C, D] = build_state_space(model, options, optical_nodes, zeta)
        % BUILD_STATE_SPACE Constructs the continuous-time state-space matrices
        % (A, B, C, D) tailored for calculating optical pointing error.
        %
        % SYNTAX:
        %   [A, B, C, D] = build_state_space(model, options, optical_nodes, zeta)
        %
        % INPUTS:
        %   model         - System model struct (must contain .M, .K, .C, .mapping, .fn)
        %   options       - Dynamic simulation setup struct (must contain .forces)
        %                   Forces fields: .node_id (String), .dof ("x", "y", or "theta")
        %   optical_nodes - Struct defining mirrors (for C matrix)
        %                   Fields: .sm_id (String), .pm_ids (Cell array of Strings)
        %   zeta          - Damping ratio (e.g., 0.01 for 1%)
        %
        % OUTPUTS:
        %   A, B, C, D    - State-space matrices for the continuous-time system
        
        % Total number of degrees of freedom
        n_full = size(model.M, 1);

        if isfield(model, 'free_dofs') && ~isempty(model.free_dofs)
            free_dofs = model.free_dofs;
        else
            free_dofs = (1:n_full)';
        end

        M_f = model.M(free_dofs, free_dofs);
        K_f = model.K(free_dofs, free_dofs);
        n = length(free_dofs);

        if nargin < 4 || isempty(zeta)
            C_damp_f = model.C(free_dofs, free_dofs);
        else
            flex_fn = model.fn(model.fn > 1e-4);
            w1 = flex_fn(1) * 2 * pi;
            w2 = flex_fn(min(2, length(flex_fn))) * 2 * pi;

            if abs(w2 - w1) > 1e-6
                W = [1/(2*w1), w1/2; 
                    1/(2*w2), w2/2];
                ab = W \ [zeta; zeta];
                C_damp_f = ab(1) * M_f + ab(2) * K_f;
            else
                C_damp_f = (2 * zeta / w1) * K_f;
            end
        end

        A = [zeros(n, n), eye(n);
            -M_f \ K_f, -M_f \ C_damp_f];

        if isfield(options, 'forces')
            forces_array = options.forces;
        else
            forces_array = options; 
        end

        num_inputs = length(forces_array);
        B_top = zeros(n, num_inputs); 
        B_bottom = zeros(n, num_inputs);

        map_ids = strtrim(string({model.mapping.node_id}));

        for i = 1:num_inputs
            target_id = strtrim(string(forces_array(i).node_id));
            dof_type = lower(strtrim(string(forces_array(i).dof)));

            idx = find(strcmp(map_ids, target_id));

            if isempty(idx)
                error('Force node "%s" not found in model mapping.', target_id);
            end

            if strcmp(dof_type, "x")
                dof_idx = model.mapping(idx).global_dofs(1);
            elseif strcmp(dof_type, "y")
                dof_idx = model.mapping(idx).global_dofs(2);
            elseif strcmp(dof_type, "theta")
                dof_idx = model.mapping(idx).global_dofs(3);
            else
                error('Invalid DOF "%s" for force node "%s". Use "x", "y", or "theta".', dof_type, target_id);
            end

            F_vec = zeros(n_full, 1);
            F_vec(dof_idx) = 1;
            B_bottom(:, i) = M_f \ F_vec(free_dofs);
        end

        B = [B_top; B_bottom];

        C_disp_full = zeros(1, n_full);

            function t_idx = get_theta(node_string)
                target = strtrim(string(node_string));
                node_idx = find(strcmp(map_ids, target));
                if isempty(node_idx)
                    error('Optical node "%s" not found in model mapping.', target);
                end
                t_idx = model.mapping(node_idx).global_dofs(3);
            end

            sm_dof = get_theta(optical_nodes.sm_id);
            C_disp_full(sm_dof) = 2 * 8.5;

            for i = 1:length(optical_nodes.pm_ids)
                pm_dof = get_theta(optical_nodes.pm_ids{i});
                C_disp_full(pm_dof) = -2;
            end

            C = [C_disp_full(free_dofs), zeros(1, n)];
            D = zeros(1, num_inputs);
        end

        %% =============================================================================== 
        function rms_error = compute_psd_rms(A, B, C, D, f_range, input_idx)
            % COMPUTE_PSD_RMS Integrates a noise Power Spectral Density (PSD) profile
            % over a frequency range using the system's transfer function to find
            % the 1-sigma RMS pointing error.
            %
            % SYNTAX:
            %   rms_error = compute_psd_rms(A, B, C, D, f_range)
            %   rms_error = compute_psd_rms(A, B, C, D, f_range, input_idx)
            %
            % INPUTS:
            %   A, B, C, D - State-space matrices for the continuous-time system.
            %   f_range    - Frequency vector in Hz (e.g., linspace(10, 100, 1000)).
            %   input_idx  - (Optional) Scalar index specifying which input column of B and D
            %                corresponds to the SM motor disturbance. Defaults to 1.
            %
            % OUTPUTS:
            %   rms_error  - Scalar double representing the 1-sigma RMS pointing error.
            %
            % DEPENDENCIES:
            %   Control System Toolbox (ss, freqresp)

            % Default to the first input if no index is provided
            if nargin < 6
                input_idx = 1;
            end

            if size(B, 2) > 1
                warning('Multiple inputs detected in B matrix. Automatically slicing to use input index %d.', input_idx);
                B = B(:, input_idx);
                D = D(:, input_idx);
            end

            sys = ss(A, B, C, D);
            w_range = 2 * pi * f_range;

            % freqresp returns a 3D array: [outputs x inputs x frequencies]
            H_complex = squeeze(freqresp(sys, w_range));

            % If there's only 1 output, squeeze makes it a column vector, we want it to be 1 x N
            if size(C, 1) == 1
                H_complex = H_complex(:)';
            end

            f_range_row = f_range(:)';
            S_xx = 0.01 ./ f_range_row;

            % Implicit expansion will multiply S_xx across all rows (outputs)
            G_yy = (abs(H_complex).^2) .* S_xx;

            % Integrate across the frequency dimension (dimension 2)
            variance = trapz(f_range_row, G_yy, 2);
            rms_error = sqrt(variance);
        end

    end
    
end

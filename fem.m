classdef fem
    methods (Static)
        function model = generate_fem_matrices(model)
        % GENERATE_FEM_MATRICES Assembles global stiffness (K) and mass (M)
        % matrices for a 2D beam/element model and maps degrees of freedom.
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
        %
        % OUTPUTS:
        %   model - Updated input struct with appended fields:
        %     .K       : Global Stiffness Matrix (3N x 3N) [N/m, N*m/rad]
        %     .M       : Global Consistent Mass Matrix (3N x 3N) [kg, kg*m^2]
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
        num_dofs = num_nodes * 3; % 3 DOFs per node in 2D: [x, y, theta]
        
        % 1. Create DOF Mapping
        node_ids = [nodes.id];
        dof_map = struct('node_id', {}, 'global_dofs', {});
        
        for i = 1:num_nodes
            dof_map(i).node_id = nodes(i).id;
            dof_map(i).global_dofs = [3*i-2, 3*i-1, 3*i];
        end
        model.mapping = dof_map;
        
        % 2. Initialize Global Matrices
        K_global = zeros(num_dofs, num_dofs);
        M_global = zeros(num_dofs, num_dofs);
        
        % 3. Assemble Element Matrices
        for i = 1:length(elements)
            el = elements(i);
            
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
            end
    
            % Transform to Global Frame
            K_e = T' * k_local * T;
            M_e = T' * m_local * T;
    
            % Identify Global DOFs for Assembly
            dofs = [dof_map(n1_idx).global_dofs, dof_map(n2_idx).global_dofs];
    
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
        
        % 5. Calculate Natural Frequencies and Mode Shapes
        % Solve the generalized eigenvalue problem K*V = w^2*M*V
        [V, D] = eig(K_global, M_global);

        % Extract eigenvalues (omega^2)
        evals = diag(D);

        % Clean up numerical noise for rigid body modes
        evals(abs(evals) < 1e-10) = 0; 

        % Calculate natural frequencies in Hz
        omega = sqrt(abs(evals));
        fn_all = omega / (2 * pi);

        % Sort frequencies and corresponding modes
        [fn_sorted, sort_idx] = sort(fn_all);
        modes_sorted = V(:, sort_idx);

        % Pack outputs back into the struct
        model.K = K_global;
        model.M = M_global;
        model.fn = fn_sorted;
        model.modes = modes_sorted;
    end
    
    % ===============================================================================
        
        function plot_fem_model(model)
        % PLOT_FEM_MODEL Visualizes the geometry, relative element stiffness,
        % and relative nodal mass/inertia distribution across 3 subplots.
        %
        % SYNTAX:
        %   fem.plot_fem_model(model)
        %
        % INPUTS:
        %   model - Struct containing:
        %     .nodes    : Array of node structs (x, y, mass, inertia)
        %     .elements : Array of element structs (nodes, E, A, k_discrete)
        %     .mapping  : Mapping struct containing global_dofs per node
        %
        % OUTPUTS:
        %   Generates a figure window containing 3 subplots:
        %     1. Geometry & Mapping (Node IDs, DOFs, Element IDs)
        %     2. Stiffness Mapping (Line width scaled by element axial/discrete stiffness)
        %     3. Mass Mapping (Bubble size scaled by node mass + rotational inertia)
        %
        % DEPENDENCIES:
        %   None (Built-in MATLAB graphics routines)

            % Extract arrays for cleaner referencing
            nodes = model.nodes;
            elements = model.elements;
            node_ids = [nodes.id];
            
        % Define specific colors
        grey_color = [0.5, 0.5, 0.5];
        stiffness_color = '#F9BE08';
        mass_color = '#8A9A5B';
        royal_blue = '#4169E1'; 
        
        % Setup Figure
        figure('Name', 'FEM Struct Visualization', 'Position', [100, 100, 1500, 500]);
        
        %% Plot 1: Pure Geometry with Labels
        subplot(1, 3, 1);
        hold on; grid on; axis equal;
        title('Geometry and Mapping');
        xlabel('X Position (m)'); ylabel('Y Position (m)');
        
        % Plot Elements (Grey)
        for i = 1:length(elements)
            el = elements(i);
            n1_idx = find(node_ids == el.nodes(1));
            n2_idx = find(node_ids == el.nodes(2));
            
            x = [nodes(n1_idx).x, nodes(n2_idx).x];
            y = [nodes(n1_idx).y, nodes(n2_idx).y];
            
            plot(x, y, '-', 'Color', grey_color, 'LineWidth', 1.5);
            
            % Label Element at midpoint
            text(mean(x), mean(y), "  " + el.id, 'Color', grey_color, 'FontSize', 9);
        end
        
        % Plot Nodes (Royal Blue) and Label DOFs
        for i = 1:length(nodes)
            scatter(nodes(i).x, nodes(i).y, 50, 'filled', 'MarkerFaceColor', royal_blue);
            
            dofs = model.mapping(i).global_dofs;
            label = sprintf(' %s\n DOFs: [%d, %d, %d]', nodes(i).id, dofs(1), dofs(2), dofs(3));
            text(nodes(i).x, nodes(i).y, label, 'VerticalAlignment', 'bottom', 'FontSize', 9);
        end
        hold off;
        
        %% Plot 2: Stiffness Mapping
        subplot(1, 3, 2);
        hold on; grid on; axis equal;
        title('Relative Element Stiffness (Axial)');
        xlabel('X Position (m)'); ylabel('Y Position (m)');
        
        % Calculate stiffness magnitudes (Using Axial Stiffness for visualization)
        k_mags = zeros(1, length(elements));
        for i = 1:length(elements)
            el = elements(i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            L = sqrt((n2.x - n1.x)^2 + (n2.y - n1.y)^2);
            
            if el.type == "beam"
                k_mags(i) = (el.E * el.A) / L;
            elseif el.type == "spring"
                k_mags(i) = el.k_spring(1); % Using Axial stiffness for plot weight
            end
        end
        
        % Normalize stiffness for line width plotting (range 1 to 5)
        k_min = min(k_mags); k_max = max(k_mags);
        if k_max == k_min
            widths = ones(1, length(elements)) * 3;
        else
            widths = 1 + 4 * (k_mags - k_min) / (k_max - k_min);
        end
        
        % Plot elements with varying line widths
        for i = 1:length(elements)
            el = elements(i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            
            if el.type == "spring"
                % Make springs dashed to differentiate them visually
                plot([n1.x, n2.x], [n1.y, n2.y], '--', 'Color', stiffness_color, 'LineWidth', widths(i));
            else
                plot([n1.x, n2.x], [n1.y, n2.y], '-', 'Color', stiffness_color, 'LineWidth', widths(i));
            end
        end
        
        for i = 1:length(nodes)
            scatter(nodes(i).x, nodes(i).y, 50, 'filled', 'MarkerFaceColor', royal_blue);
        end
        hold off;
        
        %% Plot 3: Mass Mapping
        subplot(1, 3, 3);
        hold on; grid on; axis equal;
        title('Relative Nodal Mass & Inertia');
        xlabel('X Position (m)'); ylabel('Y Position (m)');
        
        % Plot faint element lines for context
        for i = 1:length(elements)
            el = elements(i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            plot([n1.x, n2.x], [n1.y, n2.y], '-', 'Color', mass_color);
        end
        
        % Calculate mass + inertia magnitudes for scaling
        m_mags = arrayfun(@(n) n.mass + n.inertia, nodes);
        m_min = min(m_mags); m_max = max(m_mags);
        
        if m_max == m_min
            sizes = ones(1, length(nodes)) * 100;
        else
            sizes = 30 + 270 * (m_mags - m_min) / (m_max - m_min);
        end
        
        for i = 1:length(nodes)
            scatter(nodes(i).x, nodes(i).y, sizes(i), 'filled', 'MarkerFaceColor', mass_color, 'MarkerFaceAlpha', 0.6);
        end
        
        for i = 1:length(nodes)
            scatter(nodes(i).x, nodes(i).y, 30, 'filled', 'MarkerFaceColor', royal_blue);
            text(nodes(i).x, nodes(i).y, sprintf('  M:%.1f, I:%.1f', nodes(i).mass, nodes(i).inertia), 'FontSize', 8);
        end
        hold off;
    end
    
    % ===============================================================================

        function [X_rms, X] = calculate_plot_rms(model, options)
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
            %             .mapping : DOF mapping struct array
            %   options - Dynamic simulation setup struct:
            %             .t      : Time vector [s]
            %             .z      : Proportional damping ratio (e.g., 0.01 for 1%)
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
            mapping = model.mapping;
            
            t = options.t;
            z = options.z;
            n_dof = size(M, 1);
            
            % Build Force Matrix from options
            F_matrix = zeros(n_dof, length(t));
            if isfield(options, 'forces')
                for i = 1:length(options.forces)
                    f = options.forces(i);
                    
                    % Find the global dofs for this node
                    node_idx = find([mapping.node_id] == f.node_id);
                    if isempty(node_idx)
                        error('Node ID %s not found in mapping.', f.node_id);
                    end
                    
                    node_dofs = mapping(node_idx).global_dofs;
                    
                    % Map to the specific directional DOF
                    if strcmpi(f.dof, 'x')
                        g_dof = node_dofs(1);
                    elseif strcmpi(f.dof, 'y')
                        g_dof = node_dofs(2);
                    elseif strcmpi(f.dof, 'theta')
                        g_dof = node_dofs(3);
                    else
                        error('Invalid DOF type. Use ''x'', ''y'', or ''theta''.');
                    end
                    
                    % Add force profile (allows superposition if multiple forces hit one DOF)
                    F_matrix(g_dof, :) = F_matrix(g_dof, :) + f.profile;
                end
            end
            
            % Generate labels for plotting using the mapping
            dof_labels = strings(n_dof, 1);
            for i = 1:length(mapping)
                dof_labels(mapping(i).global_dofs(1)) = mapping(i).node_id + " (X)";
                dof_labels(mapping(i).global_dofs(2)) = mapping(i).node_id + " (Y)";
                dof_labels(mapping(i).global_dofs(3)) = mapping(i).node_id + " (\theta)";
            end
            
            % Calculate proportional Rayleigh Damping matrix (C = alpha*M + beta*K)
            [~, D] = eig(K, M);
            omega = sqrt(max(diag(D))); % approximate fundamental omega
            beta = 2 * z / omega;
            C = beta * K; 
            
            % Convert to State-Space Formulation
            A = [zeros(n_dof), eye(n_dof);
                 -M\K, -M\C];
            B = [zeros(n_dof);
                 M\eye(n_dof)];
            C_ss = [eye(n_dof), zeros(n_dof)]; 
            D_ss = zeros(n_dof);
            
            % Create State-Space model
            sys = ss(A, B, C_ss, D_ss);
            
            % Simulate system response using lsim
            X = lsim(sys, F_matrix', t)'; 
            
            % Calculate RMS response and maximum absolute displacement
            X_rms = sqrt(mean(X.^2, 2));
            max_abs_X = max(abs(X), [], 2);
            
            % Sort indices so we plot the largest magnitude first (background)
            [~, sort_idx] = sort(max_abs_X, 'descend');
            
            % Generate consistent colors for each DOF
            colors = lines(n_dof); 
            
            % Plot 1: Time History Response
            figure('Position', [100, 100, 800, 600]);
            subplot(2,1,1);
            hold on;
            % Plot in sorted order: largest to smallest
            for i = 1:n_dof
                idx = sort_idx(i);
                % Only plot if there is actual movement to keep legend clean
                if max_abs_X(idx) > 1e-12 
                    plot(t, X(idx, :), 'LineWidth', 1.5, 'Color', colors(idx, :), ...
                        'DisplayName', dof_labels(idx));
                end
            end
            hold off;
            title('Time History of Displacements (Layered Largest to Smallest)');
            xlabel('Time (s)');
            ylabel('Displacement (m, rad)');
            grid on;
            legend('Location', 'eastoutside');
            
            % Plot 2: RMS Bar Chart
            subplot(2,1,2);
            bar(1:n_dof, X_rms, 'FaceColor', [0.8500 0.3250 0.0980]);
            title('RMS Displacement by Node & DOF');
            xlabel('Degree of Freedom');
            ylabel('X_{RMS} (m, rad)');
            xticks(1:n_dof);
            xticklabels(dof_labels);
            xtickangle(45);
            grid on;
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
            merged_elements = [];
            
            % Store mapping from old string IDs to the merged node indices per model
            map_old_id = cell(1, nargin);
            map_new_idx = cell(1, nargin);
            
            for i = 1:nargin
                mdl = varargin{i};
                map_old_id{i} = strings(1, length(mdl.nodes));
                map_new_idx{i} = zeros(1, length(mdl.nodes));
                
                %% 1. Process and Merge Nodes
                for j = 1:length(mdl.nodes)
                    cn = mdl.nodes(j);
                    match_found = false;
                    
                    % Check against already merged nodes
                    for k = 1:length(merged_nodes)
                        dx = merged_nodes(k).x - cn.x;
                        dy = merged_nodes(k).y - cn.y;
                        
                        if sqrt(dx^2 + dy^2) <= tol
                            match_found = true;
                            
                            % Check for property mismatches
                            if abs(merged_nodes(k).mass - cn.mass) > tol
                                warning('Mass mismatch at merged node %s. Keeping primary mass: %f', merged_nodes(k).id, merged_nodes(k).mass);
                            end
                            if abs(merged_nodes(k).inertia - cn.inertia) > tol
                                warning('Inertia mismatch at merged node %s. Keeping primary inertia: %f', merged_nodes(k).id, merged_nodes(k).inertia);
                            end
                            
                            % Combine constraints (most restrictive defaults to 1 via logical OR)
                            merged_nodes(k).constraints = merged_nodes(k).constraints | cn.constraints;
                            
                            % Combine string IDs if they differ, without duplicating identical IDs
                            % Updated to use "-" instead of "_"
                            if merged_nodes(k).id ~= cn.id
                                if ~contains(merged_nodes(k).id, cn.id)
                                    merged_nodes(k).id = merged_nodes(k).id + "-" + cn.id;
                                end
                            end
                            
                            % Store index mapping for element connectivity update
                            map_old_id{i}(j) = cn.id;
                            map_new_idx{i}(j) = k;
                            break;
                        end
                    end
                    
                    % If no match within tolerance, add as a new node
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
                
                % Extract parent variable name, replace underscores, and handle fallback
                parent_name = inputname(i);
                if isempty(parent_name)
                    parent_name = "Model " + i; 
                end
                parent_name = strrep(string(parent_name), "_", " ");
                
                for j = 1:length(mdl.elements)
                    ce = mdl.elements(j);
                    
                    % Prepend prefix based on parent model name using " - "
                    ce.id = parent_name + " - " + ce.id;
                    
                    % Find updated string IDs for the connected nodes
                    n1_old_idx = find(map_old_id{i} == ce.nodes(1));
                    new_n1_id = merged_nodes(map_new_idx{i}(n1_old_idx)).id;
                    
                    n2_old_idx = find(map_old_id{i} == ce.nodes(2));
                    new_n2_id = merged_nodes(map_new_idx{i}(n2_old_idx)).id;
                    
                    % Update element's node array with the unified node IDs
                    ce.nodes = [new_n1_id, new_n2_id];
                    
                    if isempty(merged_elements)
                        merged_elements = ce;
                    else
                        merged_elements(end+1) = ce;
                    end
                end
            end
            
            %% 3. Package and Generate Global Matrices
            combined_model.nodes = merged_nodes;
            combined_model.elements = merged_elements;
            
            % Pass the consolidated struct to the generator to map DOFs and output M, K matrices
            combined_model = fem.generate_fem_matrices(combined_model);
        end
        
            % ===============================================================================
    
        function plot_mode_shapes(model, freq_range)
            % PLOT_MODE_SHAPES Visualizes structural eigenmodes.
            % Inputs:
            %   model      : FEM struct with evaluated fn and modes.
            %   freq_range : (Optional) Single frequency to find the closest match, 
            %                or [min, max] bound to plot a subset. 
            
            if nargin < 2
                freq_range = [];
            end
            
            fn = model.fn;
            modes = model.modes;
            nodes = model.nodes;
            elements = model.elements;
            mapping = model.mapping;
            node_ids = [nodes.id];
            
            % Determine which modes to plot based on freq_range
            if isempty(freq_range)
                modes_to_plot = 1:length(fn); % Plot all
            elseif length(freq_range) == 1
                % Find closest mode
                [~, min_idx] = min(abs(fn - freq_range));
                modes_to_plot = min_idx;
            elseif length(freq_range) == 2
                % Find modes within bounds
                modes_to_plot = find(fn >= freq_range(1) & fn <= freq_range(2));
                if isempty(modes_to_plot)
                    warning('No modes found in the specified frequency range.');
                    return;
                end
            else
                error('freq_range must be empty, a scalar, or a 2-element vector.');
            end
            
            % Set colors
            bg_grey = [0.7 0.7 0.7];
            moss_green = '#8A9A5B';
            
            % Calculate automatic spatial scaling for the mode shape visualizations.
            % Extract neutral node positions
            x_vals = [nodes.x];
            y_vals = [nodes.y];
            max_dim = max(max(x_vals) - min(x_vals), max(y_vals) - min(y_vals));
            if max_dim == 0; max_dim = 1; end % Fallback for single node
            
            % Limit maximum visual displacement to 15% of the model's max dimension
            target_disp = 0.15 * max_dim; 
            
            % Length of the rotation indicating line (10% of model dim)
            rot_line_len = 0.10 * max_dim; 
            
            % Create figure with a flowing tiled layout
            figure('Name', 'Mode Shapes', 'Position', [100, 100, 1000, 800]);
            t = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
            
            for i = 1:length(modes_to_plot)
                m_idx = modes_to_plot(i);
                mode_vec = modes(:, m_idx);
                
                % Calculate scaling factor for this specific mode
                % Extract only translational DOFs (x and y) to find max translation
                trans_dofs = [];
                for k = 1:length(mapping)
                    trans_dofs = [trans_dofs, mapping(k).global_dofs(1:2)];
                end
                max_trans_val = max(abs(mode_vec(trans_dofs)));
                
                if max_trans_val > 1e-12
                    scale_factor = target_disp / max_trans_val;
                else
                    scale_factor = 0; % Prevent division by zero for purely rotational modes
                end
                
                nexttile;
                hold on; grid on; axis equal;
                title(sprintf('Mode %d: %.3f Hz', m_idx, fn(m_idx)));
                
                % 1. Plot Background (Undeformed) Structure
                for e = 1:length(elements)
                    el = elements(e);
                    n1_idx = find(node_ids == el.nodes(1));
                    n2_idx = find(node_ids == el.nodes(2));
                    
                    plot([nodes(n1_idx).x, nodes(n2_idx).x], ...
                         [nodes(n1_idx).y, nodes(n2_idx).y], ...
                         'Color', bg_grey, 'LineWidth', 1);
                end
                for n = 1:length(nodes)
                    scatter(nodes(n).x, nodes(n).y, 20, 'MarkerEdgeColor', bg_grey, 'MarkerFaceColor', bg_grey);
                end
                
                % 2. Plot Deformed Mode Shape
                % Extract deformed node positions
                def_x = zeros(1, length(nodes));
                def_y = zeros(1, length(nodes));
                def_theta = zeros(1, length(nodes));
                
                for n = 1:length(nodes)
                    dofs = mapping(n).global_dofs;
                    def_x(n) = nodes(n).x + mode_vec(dofs(1)) * scale_factor;
                    def_y(n) = nodes(n).y + mode_vec(dofs(2)) * scale_factor;
                    def_theta(n) = mode_vec(dofs(3)) * scale_factor; 
                end
                
                % Draw deformed elements
                for e = 1:length(elements)
                    el = elements(e);
                    n1_idx = find(node_ids == el.nodes(1));
                    n2_idx = find(node_ids == el.nodes(2));
                    
                    if el.type == "spring"
                        line_style = '--';
                    else
                        line_style = '-';
                    end
                    
                    plot([def_x(n1_idx), def_x(n2_idx)], ...
                         [def_y(n1_idx), def_y(n2_idx)], ...
                         line_style, 'Color', moss_green, 'LineWidth', 2);
                end
                
                % Draw deformed nodes (Squares with rotation lines)
                for n = 1:length(nodes)
                    % Plot square node
                    scatter(def_x(n), def_y(n), 40, 's', 'MarkerEdgeColor', moss_green, 'MarkerFaceColor', moss_green);
                    
                    % Plot rotation indicator line
                    theta = def_theta(n);
                    dx = (rot_line_len/2) * cos(theta);
                    dy = (rot_line_len/2) * sin(theta);
                    
                    plot([def_x(n) - dx, def_x(n) + dx], ...
                         [def_y(n) - dy, def_y(n) + dy], ...
                         'k-', 'LineWidth', 1.5);
                end
                
                hold off;
            end
        end
        
    end
    
end

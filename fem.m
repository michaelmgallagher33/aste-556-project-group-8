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
        [V, D] = eig(K_global, M_global);
    
        evals = diag(D);
        evals(abs(evals) < 1e-10) = 0; 
    
        omega = sqrt(abs(evals));
        fn_all = omega / (2 * pi);
    
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
        
        grey_color = [0.5, 0.5, 0.5];
        stiffness_color = '#F9BE08';
        mass_color = '#8A9A5B';
        royal_blue = '#4169E1'; 
        
        figure('Name', 'FEM Struct Visualization', 'Position', [100, 100, 1500, 500]);
        
        %% Plot 1: Pure Geometry with Labels
        subplot(1, 3, 1);
        hold on; grid on; axis equal;
        title('Geometry and Mapping');
        xlabel('X Position (m)'); ylabel('Y Position (m)');
        
        for i = 1:length(elements)
            el = elements(i);
            n1_idx = find(node_ids == el.nodes(1));
            n2_idx = find(node_ids == el.nodes(2));
            
            x = [nodes(n1_idx).x, nodes(n2_idx).x];
            y = [nodes(n1_idx).y, nodes(n2_idx).y];
            
            plot(x, y, '-', 'Color', grey_color, 'LineWidth', 1.5);
            text(mean(x), mean(y), "  " + el.id, 'Color', grey_color, 'FontSize', 9);
        end
        
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
        
        k_mags = zeros(1, length(elements));
        for i = 1:length(elements)
            el = elements(i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            L = sqrt((n2.x - n1.x)^2 + (n2.y - n1.y)^2);
            
            if el.type == "beam"
                k_mags(i) = (el.E * el.A) / L;
            elseif el.type == "spring"
                k_mags(i) = el.k_spring(1);
            end
        end
        
        k_min = min(k_mags); k_max = max(k_mags);
        if k_max == k_min
            widths = ones(1, length(elements)) * 3;
        else
            widths = 1 + 4 * (k_mags - k_min) / (k_max - k_min);
        end
        
        for i = 1:length(elements)
            el = elements(i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            
            if el.type == "spring"
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
        
        for i = 1:length(elements)
            el = elements(i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            plot([n1.x, n2.x], [n1.y, n2.y], '-', 'Color', mass_color);
        end
        
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

            % Build Force Matrix
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

            % Proportional Rayleigh Damping matrix
            % Filter out rigid body modes (< 1e-4 Hz) to avoid singularity
            flex_modes = model.fn(model.fn > 1e-4);
            if length(flex_modes) < 2
                % Fallback if strictly rigid
                beta = 0; alpha = 0; 
            else
                w1 = flex_modes(1) * 2 * pi;
                w2 = flex_modes(2) * 2 * pi;
                W = [1/(2*w1), w1/2; 1/(2*w2), w2/2];
                ab = W \ [z; z];
                alpha = ab(1);
                beta = ab(2);
            end
            C = alpha * M + beta * K; 

            % Convert to State-Space Formulation
            A = [zeros(n_dof), eye(n_dof);
                -M\K, -M\C];
            B = [zeros(n_dof);
                M\eye(n_dof)];
            C_ss = [eye(n_dof), zeros(n_dof)]; 
            D_ss = zeros(n_dof);

            sys = ss(A, B, C_ss, D_ss);

            % Simulate response
            X = lsim(sys, F_matrix', t)'; 

            % Calculate RMS response
            X_rms = sqrt(mean(X.^2, 2));
        end

        % ===============================================================================

        function plot_rms(model, options, X, X_rms, plot_ids)
            % PLOT_RMS Plots dynamic time-history displacement profiles categorized
            % separately by degree of freedom (X, Y, Theta), along with an RMS bar chart.
            %
            % SYNTAX:
            %   fem.plot_rms(model, options, X, X_rms)
            %   fem.plot_rms(model, options, X, X_rms, plot_ids)
            %
            % INPUTS:
            %   model    - Struct containing .mapping
            %   options  - Struct containing .t (time vector)
            %   X        - Time-history displacement matrix (n_dof x length(t)) from calculate_rms
            %   X_rms    - RMS response vector (n_dof x 1) from calculate_rms
            %   plot_ids - (Optional) Cell array or string array of node IDs to isolate in the plots.
            %
            % OUTPUTS:
            %   Generates a figure window containing 4 subplots:
            %     1. X-Translation Time History
            %     2. Y-Translation Time History
            %     3. Theta-Rotation Time History
            %     4. Combined RMS Displacement Bar Chart

            if nargin < 5 || isempty(plot_ids)
                filter_nodes = false;
            else
                filter_nodes = true;
                plot_ids = string(plot_ids);
            end
    
            mapping = model.mapping;
            t = options.t;
            n_dof = size(X, 1);
            
            % Extract DOFs categorized by direction
            x_dofs = arrayfun(@(m) m.global_dofs(1), mapping);
            y_dofs = arrayfun(@(m) m.global_dofs(2), mapping);
            theta_dofs = arrayfun(@(m) m.global_dofs(3), mapping);
            
            % Generate labels
            dof_labels = strings(n_dof, 1);
            for i = 1:length(mapping)
                dof_labels(mapping(i).global_dofs(1)) = mapping(i).node_id + " (X)";
                dof_labels(mapping(i).global_dofs(2)) = mapping(i).node_id + " (Y)";
                dof_labels(mapping(i).global_dofs(3)) = mapping(i).node_id + " (\theta)";
            end
            
            colors = lines(length(mapping));
            
            figure('Name', 'Dynamic Response Analysis', 'Position', [100, 100, 1200, 900]);
            
            %% Subplot 1: X Displacements
            subplot(2, 2, 1);
            hold on; grid on;
            title('X-Translation Time History');
            xlabel('Time (s)'); ylabel('Displacement (m)');
            for i = 1:length(x_dofs)
                node_id = mapping(i).node_id;
                if filter_nodes && ~ismember(node_id, plot_ids)
                    continue;
                end
                
                dof = x_dofs(i);
                if max(abs(X(dof, :))) > 1e-12
                    plot(t, X(dof, :), 'LineWidth', 1.5, 'Color', colors(i, :), ...
                        'DisplayName', node_id);
                end
            end
            legend('Location', 'best');
            hold off;
            
            %% Subplot 2: Y Displacements
            subplot(2, 2, 2);
            hold on; grid on;
            title('Y-Translation Time History');
            xlabel('Time (s)'); ylabel('Displacement (m)');
            for i = 1:length(y_dofs)
                node_id = mapping(i).node_id;
                if filter_nodes && ~ismember(node_id, plot_ids)
                    continue;
                end
                
                dof = y_dofs(i);
                if max(abs(X(dof, :))) > 1e-12
                    plot(t, X(dof, :), 'LineWidth', 1.5, 'Color', colors(i, :), ...
                        'DisplayName', node_id);
                end
            end
            legend('Location', 'best');
            hold off;
            
            %% Subplot 3: Theta Rotations
            subplot(2, 2, 3);
            hold on; grid on;
            title('\Theta-Rotation Time History');
            xlabel('Time (s)'); ylabel('Rotation (rad)');
            for i = 1:length(theta_dofs)
                node_id = mapping(i).node_id;
                if filter_nodes && ~ismember(node_id, plot_ids)
                    continue;
                end
                
                dof = theta_dofs(i);
                if max(abs(X(dof, :))) > 1e-12
                    plot(t, X(dof, :), 'LineWidth', 1.5, 'Color', colors(i, :), ...
                        'DisplayName', node_id);
                end
            end
            legend('Location', 'best');
            hold off;
            
            %% Subplot 4: RMS Summary Bar Chart
            subplot(2, 2, 4);
            
            if filter_nodes
                dofs_to_keep = [];
                for i = 1:length(mapping)
                    if ismember(mapping(i).node_id, plot_ids)
                        dofs_to_keep = [dofs_to_keep, mapping(i).global_dofs];
                    end
                end
                bar_data = X_rms(dofs_to_keep);
                bar_labels = dof_labels(dofs_to_keep);
            else
                bar_data = X_rms;
                bar_labels = dof_labels;
            end
            
            bar(1:length(bar_data), bar_data, 'FaceColor', [0.8500 0.3250 0.0980]);
            title('RMS Displacement by Node & DOF');
            xlabel('Degree of Freedom');
            ylabel('X_{RMS} (m, rad)');
            xticks(1:length(bar_data));
            xticklabels(bar_labels);
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
                    ce = mdl.elements(j);

                    ce.id = parent_name + " - " + ce.id;

                    n1_old_idx = find(map_old_id{i} == ce.nodes(1));
                    new_n1_id = merged_nodes(map_new_idx{i}(n1_old_idx)).id;

                    n2_old_idx = find(map_old_id{i} == ce.nodes(2));
                    new_n2_id = merged_nodes(map_new_idx{i}(n2_old_idx)).id;

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

            combined_model = fem.generate_fem_matrices(combined_model);
        end
        
        %% =============================================================================== 
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
        
        if isempty(freq_range)
            modes_to_plot = 1:length(fn); 
        elseif length(freq_range) == 1
            [~, min_idx] = min(abs(fn - freq_range));
            modes_to_plot = min_idx;
        elseif length(freq_range) == 2
            modes_to_plot = find(fn >= freq_range(1) & fn <= freq_range(2));
            if isempty(modes_to_plot)
                warning('No modes found in the specified frequency range.');
                return;
            end
        else
            error('freq_range must be empty, a scalar, or a 2-element vector.');
        end
        
        bg_grey = [0.7 0.7 0.7];
        moss_green = '#8A9A5B';
        
        x_vals = [nodes.x];
        y_vals = [nodes.y];
        max_dim = max(max(x_vals) - min(x_vals), max(y_vals) - min(y_vals));
        if max_dim == 0; max_dim = 1; end 
        
        target_disp = 0.15 * max_dim; 
        rot_line_len = 0.10 * max_dim; 
        
        figure('Name', 'Mode Shapes', 'Position', [100, 100, 1000, 800]);
        t_layout = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
        
        for i = 1:length(modes_to_plot)
            m_idx = modes_to_plot(i);
            mode_vec = modes(:, m_idx);
            
            trans_dofs = [];
            for k = 1:length(mapping)
                trans_dofs = [trans_dofs, mapping(k).global_dofs(1:2)];
            end
            max_trans_val = max(abs(mode_vec(trans_dofs)));
            
            if max_trans_val > 1e-12
                scale_factor = target_disp / max_trans_val;
            else
                scale_factor = 0; 
            end
            
            nexttile;
            hold on; grid on; axis equal;
            title(sprintf('Mode %d: %.3f Hz', m_idx, fn(m_idx)));
            
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
            
            def_x = zeros(1, length(nodes));
            def_y = zeros(1, length(nodes));
            def_theta = zeros(1, length(nodes));
            
            for n = 1:length(nodes)
                dofs = mapping(n).global_dofs;
                def_x(n) = nodes(n).x + mode_vec(dofs(1)) * scale_factor;
                def_y(n) = nodes(n).y + mode_vec(dofs(2)) * scale_factor;
                def_theta(n) = mode_vec(dofs(3)) * scale_factor; 
            end
            
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
            
            for n = 1:length(nodes)
                scatter(def_x(n), def_y(n), 40, 's', 'MarkerEdgeColor', moss_green, 'MarkerFaceColor', moss_green);
                
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

            function idx = get_theta_idx(target_id)
                all_ids = string({mapping.node_id});
                target_str = string(target_id);

                is_match = strcmp(all_ids, target_str);

                if ~any(is_match)
                    error('OpticalNode:NotFound', ...
                        'The optical node "%s" was not found in the model mapping. Please check for typos.', target_str);
                end

                matched_node = mapping(is_match);
                idx = matched_node(1).global_dofs(3); 
            end

            sm_idx = get_theta_idx(optical_nodes.sm_id);
            theta_sm = X(sm_idx, :);
            error_sm = 17 * theta_sm;

            num_pms = length(optical_nodes.pm_ids);
            num_time_steps = size(X, 2);

            theta_pms = zeros(num_pms, num_time_steps);
            error_pms = zeros(num_pms, num_time_steps);
            error_total_rad = zeros(num_pms, num_time_steps);

            for i = 1:num_pms
                pm_idx = get_theta_idx(optical_nodes.pm_ids{i});
                theta_pms(i, :) = X(pm_idx, :);
                error_pms(i, :) = 2 * theta_pms(i, :);

                error_total_rad(i, :) = abs(error_pms(i, :)) + abs(error_sm);
            end

            % Convert output to arc-seconds (1 rad = (180/pi)*3600 arc-sec)
            rad_to_arcsec = (180 / pi) * 3600;

            components.sm = error_sm * rad_to_arcsec;
            components.pms = error_pms * rad_to_arcsec;
            components.total = error_total_rad * rad_to_arcsec;

            error_total = components.total;
        end

        %% =============================================================================== 
        function plot_pointing_error_requirement(t, error_total, threshold, window_size)
            % PLOT_POINTING_ERROR_REQUIREMENT Plots system error vs requirement and evaluates 
            % compliance in specific time windows.
            %
            % SYNTAX:
            %   fem.plot_pointing_error_requirement(t, error_total, threshold, window_size)

            if nargin < 3 || isempty(threshold)
                threshold = 1.88e-3; % arc-seconds requirement
            end
            if nargin < 4 || isempty(window_size)
                window_size = 30; % seconds
            end

            % If there are multiple payload mirrors, take the worst-case error at each step
            if size(error_total, 1) > 1
                max_error = max(abs(error_total), [], 1);
            else
                max_error = abs(error_total);
            end

            % Calculate compliance in time windows
            t_max = t(end);
            num_windows = max(1, floor(t_max / window_size));
            meets_req = zeros(1, num_windows);
            window_starts = 0:window_size:(num_windows-1)*window_size;
            window_ends = window_starts + window_size;

            for w = 1:num_windows
                idx = (t >= window_starts(w)) & (t <= window_ends(w));
                % Requirement: Error must strictly stay below the threshold in the window
                if max(max_error(idx)) <= threshold
                    meets_req(w) = 1;
                end
            end

            compliance_pct = (sum(meets_req) / num_windows) * 100;

            % Plotting
            figure('Name', 'Optical Pointing Error vs Requirement', 'Position', [150, 150, 1000, 700]);

            % Top Plot: Error over time
            subplot(2, 1, 1);
            hold on; grid on;
            plot(t, max_error, 'b', 'LineWidth', 1.5, 'DisplayName', 'Absolute System Error');
            yline(threshold, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('Requirement (%.2e arcsec)', threshold));
            title('Overall System Pointing Error Over Time');
            xlabel('Time (s)');
            ylabel('Absolute Error (arcsec)');
            legend('Location', 'best');
            hold off;

            % Bottom Plot: Window Compliance
            subplot(2, 1, 2);
            hold on; grid on;

            % Plot boolean stairs
            stairs([window_starts, window_ends(end)], [meets_req, meets_req(end)], 'k', 'LineWidth', 1.5);

            % Shade passing/failing windows
            for w = 1:num_windows
                if meets_req(w)
                    patch([window_starts(w) window_ends(w) window_ends(w) window_starts(w)], ...
                        [0 0 1 1], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                else
                    patch([window_starts(w) window_ends(w) window_ends(w) window_starts(w)], ...
                        [0 0 1 1], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
            end

            ylim([-0.2, 1.2]);
            yticks([0, 1]);
            yticklabels({'Fail', 'Pass'});
            title(sprintf('Requirement Compliance (%d-Second Windows)', window_size));
            xlabel('Time (s)');
            ylabel('Compliance Status');

            % Add Percentage Text Overlay
            text_x = (window_starts(1) + window_ends(end)) / 2;
            text(text_x, 0.5, sprintf('Total Compliance: %.1f%%', compliance_pct), ...
                'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                'BackgroundColor', 'w', 'EdgeColor', 'k');

            hold off;
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
        %   model         - System model struct (must contain .M, .K, .mapping, .fn)
        %   options       - Dynamic simulation setup struct (must contain .forces)
        %                   Forces fields: .node_id (String), .dof ("x", "y", or "theta")
        %   optical_nodes - Struct defining mirrors (for C matrix)
        %                   Fields: .sm_id (String), .pm_ids (Cell array of Strings)
        %   zeta          - Damping ratio (e.g., 0.01 for 1%)
        %
        % OUTPUTS:
        %   A, B, C, D    - State-space matrices for the continuous-time system
        
        % Total number of degrees of freedom
        n = size(model.M, 1);

        % Filter out rigid body modes to compute Rayleigh damping
        flex_modes = model.fn(model.fn > 1e-4);
        if length(flex_modes) < 2
            error('Not enough flexible modes to compute Rayleigh damping.');
        end
        w1 = flex_modes(1) * 2 * pi;
        w2 = flex_modes(2) * 2 * pi;

        W = [1/(2*w1), w1/2; 
            1/(2*w2), w2/2];
        ab = W \ [zeta; zeta];
        alpha = ab(1);
        beta = ab(2);

        C_damp = alpha * model.M + beta * model.K;

        A = [zeros(n, n), eye(n);
            -model.M \ model.K, -model.M \ C_damp];

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

            F_vec = zeros(n, 1);
            F_vec(dof_idx) = 1;
            B_bottom(:, i) = model.M \ F_vec;
        end

        B = [B_top; B_bottom];

        % Build C matrix with one row per Primary Mirror
            function t_idx = get_theta(node_string)
                target = strtrim(string(node_string));
                node_idx = find(strcmp(map_ids, target));
                if isempty(node_idx)
                    error('Optical node "%s" not found in model mapping.', target);
                end
                t_idx = model.mapping(node_idx).global_dofs(3);
            end

            num_pms = length(optical_nodes.pm_ids);
            C_disp = zeros(num_pms, n);
            sm_dof = get_theta(optical_nodes.sm_id);

            for i = 1:num_pms
                pm_dof = get_theta(optical_nodes.pm_ids{i});
                C_disp(i, sm_dof) = 17; % 2 * 8.5 for SM
                C_disp(i, pm_dof) = 2;
            end

            C = [C_disp, zeros(num_pms, n)];
            D = zeros(num_pms, num_inputs);
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

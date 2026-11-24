classdef viz
    methods (Static)   
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
        num_elements = length(elements);
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
        
        for i = 1:num_elements
            el = fem.get_element(elements, i);
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
        
        k_mags = zeros(1, num_elements);
        for i = 1:num_elements
            el = fem.get_element(elements, i);
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
            widths = ones(1, num_elements) * 3;
        else
            widths = 1 + 4 * (k_mags - k_min) / (k_max - k_min);
        end
        
        for i = 1:num_elements
            el = fem.get_element(elements, i);
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
        
        for i = 1:num_elements
            el = fem.get_element(elements, i);
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
        
        %% =============================================================================== 
        function plot_mode_shapes(model, freq_range)
        if nargin < 2
            freq_range = [];
        end
        
        fn = model.fn;
        modes = model.modes;
        nodes = model.nodes;
        elements = model.elements;
        num_elements = length(elements);
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
            
            % Title Generation with Effective Mass Percentage
            title_str = sprintf('Mode %d: %.3f Hz', m_idx, fn(m_idx));
            if isfield(model, 'effective_mass_pct') && ~isempty(model.effective_mass_pct)
                emp = model.effective_mass_pct(m_idx, :);
                title_str = [title_str, sprintf('\nEff Mass %% [X:%.1f, Y:%.1f, \\theta:%.1f]', emp(1), emp(2), emp(3))];
            end
            title(title_str);
            
            for e = 1:num_elements
                el = fem.get_element(elements, e);
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
            
            for e = 1:num_elements
                el = fem.get_element(elements, e);
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

            % Ensure error is a 1D absolute time history vector
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

    end
    
end

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
        num_nodes = length(nodes);
        node_ids = [nodes.id];
    
        % Color palette
        grey_color    = [0.5, 0.5, 0.5];
        axial_color   = '#F9BE08';
        shear_color   = '#E05A47';
        rot_color     = '#7B68EE';
        mass_color    = '#8A9A5B';
        el_mass_color = '#2E8B57';
        royal_blue    = '#4169E1';
    
        % Precompute element midpoints and properties
        k_ax  = zeros(1, num_elements);
        k_sh  = zeros(1, num_elements);
        k_rot = zeros(1, num_elements);
        m_el  = zeros(1, num_elements);
        el_xm = zeros(1, num_elements);
        el_ym = zeros(1, num_elements);
    
        for i = 1:num_elements
            el = fem.get_element(elements, i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            L  = hypot(n2.x - n1.x, n2.y - n1.y);
    
            el_xm(i) = 0.5 * (n1.x + n2.x);
            el_ym(i) = 0.5 * (n1.y + n2.y);
    
            if el.type == "beam" && L > 0
                k_ax(i)  = (el.E * el.A) / L;
                k_sh(i)  = (12 * el.E * el.I) / (L^3);
                k_rot(i) = (4 * el.E * el.I) / L;
                m_el(i)  = el.rho * L;
            elseif el.type == "spring"
                k_ax(i)  = el.k_spring(1);
                k_sh(i)  = el.k_spring(2);
                k_rot(i) = el.k_spring(3);
                m_el(i)  = 0;
            end
        end
    
        node_x = [nodes.x];
        node_y = [nodes.y];
        bbox = [min(node_x), max(node_x), min(node_y), max(node_y)];
    
        % Visual scaling
        ax_widths  = scale_visual(k_ax,  1.2, 6.0, true);
        sh_sizes   = scale_visual(k_sh,  30,  220, true);
        rot_sizes  = scale_visual(k_rot, 60,  360, true);
        m_el_width = scale_visual(m_el,  1.2, 6.0, false);
        m_el_sizes = scale_visual(m_el,  30,  200, false);
    
        % Dynamic figure height based on model size (clamped between 600 and 900 px)
        fig_h = min(max(600, 450 + 18 * (num_elements + num_nodes)), 900);
        fig_w = 1050;
    
        %% Figure 1: Pure Geometry & Mapping
        figure('Name', 'FEM - Geometry and Mapping', 'Position', [50, 80, fig_w, fig_h]);
        ax1 = axes;
        hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
        title(ax1, 'Geometry and Mapping');
        xlabel(ax1, 'X Position (m)'); ylabel(ax1, 'Y Position (m)');
    
        geom_labels = strings(1, num_elements + num_nodes);
        geom_colors = cell(1, num_elements + num_nodes);
    
        for i = 1:num_elements
            el = fem.get_element(elements, i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
            plot(ax1, [n1.x, n2.x], [n1.y, n2.y], '-', 'Color', grey_color, 'LineWidth', 1.5);
    
            geom_labels(i) = string(el.id);
            geom_colors{i} = grey_color;
        end
    
        for i = 1:num_nodes
            scatter(ax1, nodes(i).x, nodes(i).y, 50, 'filled', 'MarkerFaceColor', royal_blue);
            dofs = model.mapping(i).global_dofs;
            idx = num_elements + i;
            geom_labels(idx) = sprintf('%s\nDOFs: [%d, %d, %d]', nodes(i).id, dofs(1), dofs(2), dofs(3));
            geom_colors{idx} = royal_blue;
        end
    
        draw_callout_labels(ax1, [el_xm, node_x], [el_ym, node_y], geom_labels, geom_colors, bbox);
        hold(ax1, 'off');
    
        %% Figure 2: Element Stiffness Mapping
        figure('Name', 'FEM - Element Stiffness', 'Position', [100, 60, fig_w, fig_h]);
        ax2 = axes;
        hold(ax2, 'on'); grid(ax2, 'on'); axis(ax2, 'equal');
        title(ax2, 'Element Stiffness (Axial, Shear, Rotational)');
        xlabel(ax2, 'X Position (m)'); ylabel(ax2, 'Y Position (m)');
    
        h_ax  = plot(ax2, nan, nan, '-', 'Color', axial_color, 'LineWidth', 3, 'DisplayName', 'Axial (Line Width)');
        h_sh  = scatter(ax2, nan, nan, 70, 'd', 'filled', 'MarkerFaceColor', shear_color, 'DisplayName', 'Shear (Diamond)');
        h_rot = scatter(ax2, nan, nan, 100, 'o', 'MarkerEdgeColor', rot_color, 'LineWidth', 1.5, 'DisplayName', 'Rotational (Ring)');
    
        stiff_labels = strings(1, num_elements);
        stiff_colors = cell(1, num_elements);
    
        for i = 1:num_elements
            el = fem.get_element(elements, i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
    
            line_style = '-';
            if el.type == "spring", line_style = '--'; end
    
            plot(ax2, [n1.x, n2.x], [n1.y, n2.y], line_style, 'Color', axial_color, 'LineWidth', ax_widths(i), 'HandleVisibility', 'off');
            if k_sh(i) > 0
                scatter(ax2, el_xm(i), el_ym(i), sh_sizes(i), 'd', 'filled', ...
                    'MarkerFaceColor', shear_color, 'MarkerFaceAlpha', 0.75, 'HandleVisibility', 'off');
            end
            if k_rot(i) > 0
                scatter(ax2, el_xm(i), el_ym(i), rot_sizes(i), 'o', ...
                    'MarkerEdgeColor', rot_color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
    
            stiff_labels(i) = sprintf('%s\nAx: %.1e | Sh: %.1e\nRot: %.1e', el.id, k_ax(i), k_sh(i), k_rot(i));
            stiff_colors{i} = grey_color;
        end
    
        for i = 1:num_nodes
            scatter(ax2, nodes(i).x, nodes(i).y, 40, 'filled', 'MarkerFaceColor', royal_blue, 'HandleVisibility', 'off');
        end
    
        draw_callout_labels(ax2, el_xm, el_ym, stiff_labels, stiff_colors, bbox);
        legend(ax2, [h_ax, h_sh, h_rot], 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 9);
        hold(ax2, 'off');
    
        %% Figure 3: Nodal & Element Mass Mapping
        figure('Name', 'FEM - Nodal & Element Mass', 'Position', [150, 40, fig_w, fig_h]);
        ax3 = axes;
        hold(ax3, 'on'); grid(ax3, 'on'); axis(ax3, 'equal');
        title(ax3, 'Nodal & Element Mass Distribution');
        xlabel(ax3, 'X Position (m)'); ylabel(ax3, 'Y Position (m)');
    
        h_nm = scatter(ax3, nan, nan, 90, 'o', 'filled', 'MarkerFaceColor', mass_color, 'MarkerFaceAlpha', 0.6, 'DisplayName', 'Nodal Mass+Inertia');
        h_em = scatter(ax3, nan, nan, 70, 's', 'filled', 'MarkerFaceColor', el_mass_color, 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Element Mass');
    
        mass_labels = strings(1, num_elements + num_nodes);
        mass_colors = cell(1, num_elements + num_nodes);
    
        for i = 1:num_elements
            el = fem.get_element(elements, i);
            n1 = nodes(node_ids == el.nodes(1));
            n2 = nodes(node_ids == el.nodes(2));
    
            if m_el(i) > 0
                plot(ax3, [n1.x, n2.x], [n1.y, n2.y], '-', 'Color', el_mass_color, 'LineWidth', m_el_width(i), 'HandleVisibility', 'off');
                scatter(ax3, el_xm(i), el_ym(i), m_el_sizes(i), 's', 'filled', ...
                    'MarkerFaceColor', el_mass_color, 'MarkerFaceAlpha', 0.7, 'HandleVisibility', 'off');
                mass_colors{i} = el_mass_color;
            else
                plot(ax3, [n1.x, n2.x], [n1.y, n2.y], '--', 'Color', grey_color, 'LineWidth', 1.0, 'HandleVisibility', 'off');
                mass_colors{i} = grey_color;
            end
            mass_labels(i) = sprintf('%s\nm_el: %.2f kg', el.id, m_el(i));
        end
    
        m_mags = arrayfun(@(n) n.mass + n.inertia, nodes);
        node_sizes = scale_visual(m_mags, 30, 300, false);
    
        for i = 1:num_nodes
            scatter(ax3, nodes(i).x, nodes(i).y, node_sizes(i), 'filled', ...
                'MarkerFaceColor', mass_color, 'MarkerFaceAlpha', 0.6, 'HandleVisibility', 'off');
            scatter(ax3, nodes(i).x, nodes(i).y, 30, 'filled', ...
                'MarkerFaceColor', royal_blue, 'HandleVisibility', 'off');
    
            idx = num_elements + i;
            mass_labels(idx) = sprintf('%s\nM: %.1f kg | I: %.1f', nodes(i).id, nodes(i).mass, nodes(i).inertia);
            mass_colors{idx} = royal_blue;
        end
    
        draw_callout_labels(ax3, [el_xm, node_x], [el_ym, node_y], mass_labels, mass_colors, bbox);
        legend(ax3, [h_nm, h_em], 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 9);
        hold(ax3, 'off');
    end
    
    %% Local Helper: External Callout Labels with Leader Lines
    function draw_callout_labels(ax, anchors_x, anchors_y, labels, colors, bbox)
        % Places labels in sorted left/right columns outside the model bounding box
        % with leader lines pointing back to each (anchor_x, anchor_y).
        x_min = bbox(1); x_max = bbox(2);
        y_min = bbox(3); y_max = bbox(4);
    
        dx = max(x_max - x_min, 1e-3);
        dy = max(y_max - y_min, 1e-3);
        x_mid = 0.5 * (x_min + x_max);
    
        % Horizontal margin where text columns sit
        margin_x = 0.25 * dx;
        elbow_offset = 0.08 * dx;
    
        % Split into left and right groups
        left_idx  = find(anchors_x < x_mid);
        right_idx = find(anchors_x >= x_mid);
    
        % Helper to lay out one side
        function place_side(idxs, is_left)
            if isempty(idxs), return; end
    
            % Sort bottom-to-top by anchor Y to prevent leader-line crossings
            [~, order] = sort(anchors_y(idxs));
            idxs = idxs(order);
            n = numel(idxs);
    
            % Spread labels evenly across expanded vertical range
            y_pad = 0.15 * dy;
            y_targets = linspace(y_min - y_pad, y_max + y_pad, max(n, 2));
            if n == 1
                y_targets = anchors_y(idxs(1));
            end
    
            if is_left
                x_text  = x_min - margin_x;
                x_elbow = x_text + elbow_offset;
                h_align = 'right';
            else
                x_text  = x_max + margin_x;
                x_elbow = x_text - elbow_offset;
                h_align = 'left';
            end
    
            for k = 1:n
                idx = idxs(k);
                ax_x = anchors_x(idx);
                ax_y = anchors_y(idx);
                ty   = y_targets(k);
    
                % Leader line: anchor -> elbow -> text
                plot(ax, [ax_x, x_elbow, x_text], [ax_y, ty, ty], '-', ...
                    'Color', [0.6, 0.6, 0.6, 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
                % Small anchor dot
                plot(ax, ax_x, ax_y, '.', 'Color', colors{idx}, 'MarkerSize', 6, 'HandleVisibility', 'off');
    
                % External label
                text(ax, x_text, ty, " " + labels(idx) + " ", ...
                    'FontSize', 7.5, 'Color', colors{idx}, ...
                    'HorizontalAlignment', h_align, 'VerticalAlignment', 'middle', ...
                    'Interpreter', 'none');
            end
        end
    
        place_side(left_idx, true);
        place_side(right_idx, false);
    
        % Expand axis limits so external text columns are never clipped
        xlim(ax, [x_min - 0.85 * dx, x_max + 0.85 * dx]);
        ylim(ax, [y_min - 0.25 * dy, y_max + 0.25 * dy]);
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

% === Local Functions ===
% ===============================================================================
function draw_callout_labels(ax, anchors_x, anchors_y, labels, colors, bbox)
    % Places labels in sorted left/right columns outside the model bounding box
    % with leader lines pointing back to each (anchor_x, anchor_y).
    x_min = bbox(1); x_max = bbox(2);
    y_min = bbox(3); y_max = bbox(4);
    
    dx = max(x_max - x_min, 1e-3);
    dy = max(y_max - y_min, 1e-3);
    x_mid = 0.5 * (x_min + x_max);
    
    % Horizontal margin where text columns sit
    margin_x = 0.25 * dx;
    elbow_offset = 0.08 * dx;
    
    % Split into left and right groups
    left_idx  = find(anchors_x < x_mid);
    right_idx = find(anchors_x >= x_mid);
    
    % Helper to lay out one side
        function place_side(idxs, is_left)
            if isempty(idxs), return; end
    
            % Sort bottom-to-top by anchor Y to prevent leader-line crossings
            [~, order] = sort(anchors_y(idxs));
            idxs = idxs(order);
            n = numel(idxs);
    
            % Spread labels evenly across expanded vertical range
            y_pad = 0.15 * dy;
            y_targets = linspace(y_min - y_pad, y_max + y_pad, max(n, 2));
            if n == 1
                y_targets = anchors_y(idxs(1));
            end
    
            if is_left
                x_text  = x_min - margin_x;
                x_elbow = x_text + elbow_offset;
                h_align = 'right';
            else
                x_text  = x_max + margin_x;
                x_elbow = x_text - elbow_offset;
                h_align = 'left';
            end
    
            for k = 1:n
                idx = idxs(k);
                ax_x = anchors_x(idx);
                ax_y = anchors_y(idx);
                ty   = y_targets(k);
    
                % Leader line: anchor -> elbow -> text
                plot(ax, [ax_x, x_elbow, x_text], [ax_y, ty, ty], '-', ...
                    'Color', [0.6, 0.6, 0.6, 0.7], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
                % Small anchor dot
                plot(ax, ax_x, ax_y, '.', 'Color', colors{idx}, 'MarkerSize', 6, 'HandleVisibility', 'off');
    
                % External label
                text(ax, x_text, ty, " " + labels(idx) + " ", ...
                    'FontSize', 7.5, 'Color', colors{idx}, ...
                    'HorizontalAlignment', h_align, 'VerticalAlignment', 'middle', ...
                    'Interpreter', 'none');
            end
        end
    
    place_side(left_idx, true);
    place_side(right_idx, false);
    
    % Expand axis limits so external text columns are never clipped
    xlim(ax, [x_min - 0.85 * dx, x_max + 0.85 * dx]);
    ylim(ax, [y_min - 0.25 * dy, y_max + 0.25 * dy]);
end

%% Local Helper: Visual Scaling
function scaled = scale_visual(vals, min_out, max_out, use_log)
    scaled = ones(size(vals)) * min_out;
    pos_mask = vals > 0;
    if ~any(pos_mask), return; end
    
    v = vals(pos_mask);
    if use_log, v = log10(v); end
    
    v_min = min(v); v_max = max(v);
    if v_max == v_min
        scaled(pos_mask) = 0.5 * (min_out + max_out);
    else
        scaled(pos_mask) = min_out + (max_out - min_out) * (v - v_min) / (v_max - v_min);
    end
end
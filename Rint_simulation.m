%% Ontic Geometric Entanglement – 3D Visualization
% Cross-curvature R_int(theta_A, theta_B, t)
% Diffusion + localized source (models a CNOT-like gate).

clear; clc; close all;

%% 1. Grid and parameters
N  = 80;                 % Angular resolution in theta_A, theta_B
L  = pi;                 % Domain [0, pi] x [0, pi]
dx = L / N;
x  = linspace(0, L, N);  % theta_A
y  = linspace(0, L, N);  % theta_B
[X, Y] = meshgrid(x, y);

% Diffusion / decay parameters
alpha = 0.02;            % Diffusion strength (stable with dt=0.01)
gamma = 0.01;            % Small decay toward symmetric background

dt       = 0.01;                      % Time step
T_total  = 200;                       % Number of time steps
t_vec    = (0:T_total-1) * dt;        % Physical time array

% Cross-curvature field
R_int = zeros(N, N);

%% 2. Localized interaction source (CNOT region)
interaction_center = [0.75*pi, 0.5*pi];  % (theta_A, theta_B) location
interaction_width  = 0.3;

Source_Profile = exp( -((X - interaction_center(1)).^2 + ...
                        (Y - interaction_center(2)).^2) / interaction_width^2 );

kappa = 2.5;              % Coupling strength (height of injected curvature)

% Gate timing in physical time units
gate_start_time = 0.1;    % When the gate turns on
gate_end_time   = 0.4;    % When the gate turns off

%% 3. Snapshot settings (times at which to save figures)
snapshot_times = [0.05, 0.25, 0.60, 1.80];  % [pre-gate, gate-on, early diff., late diff.]
snapshot_files = { ...
    'Rint_t005_pre_gate.png', ...
    'Rint_t025_gate_on.png', ...
    'Rint_t060_diffusion.png', ...
    'Rint_t180_late.png' ...
};
snapshot_taken = false(size(snapshot_times));  % track which ones are saved
time_tolerance = dt / 2;                       % matching tolerance for t

%% 4. Setup figure
fig = figure('Color', 'k', 'Position', [100, 100, 1100, 800]);

for n = 1:T_total
    t = t_vec(n);

    %--------------------------------------------------------------
    % 4.1 Laplacian with periodic boundary conditions
    %--------------------------------------------------------------
    R_up    = [R_int(2:end, :);  R_int(end, :)];
    R_down  = [R_int(1, :);      R_int(1:end-1, :)];
    R_left  = [R_int(:, 2:end),  R_int(:, end)];
    R_right = [R_int(:, 1),      R_int(:, 1:end-1)];
    Laplacian_R = (R_up + R_down + R_left + R_right - 4*R_int) / dx^2;

    %--------------------------------------------------------------
    % 4.2 Source term: gate ON/OFF
    %--------------------------------------------------------------
    if t >= gate_start_time && t <= gate_end_time
        Current_Source = kappa * Source_Profile;
        status_msg     = 'GATE ON: Injection Phase';
        title_color    = [1 0.2 0.2];   % warm red
    else
        Current_Source = 0;
        status_msg     = 'GATE OFF: Diffusion / Relaxation';
        title_color    = [0.6 0.7 1];   % cool blue
    end

    %--------------------------------------------------------------
    % 4.3 Time update: diffusion + decay + source
    %     dR/dt = alpha * Laplacian(R) - gamma * R + Source
    %--------------------------------------------------------------
    R_int = R_int + dt * (alpha * Laplacian_R - gamma * R_int + Current_Source);

    %--------------------------------------------------------------
    % 4.4 3D Visualization of |R_int|
    %--------------------------------------------------------------
    clf(fig);
    R_plot = abs(R_int);     % We interpret this as |R_int|

    s = surf(X, Y, R_plot);
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    colormap(hot);

    % Lighting for depth
    light('Position',[-1 -1 2],'Style','local');
    lighting gouraud;
    material shiny;

    % View and axis settings
    view(-45, 30);
    axis tight;

    % Dynamical z-limits: show growth without saturating too quickly
    maxR = max(R_plot(:));
    zmax = max(0.5, min(3, 1.1 * maxR));   % clamp between 0.5 and 3
    zlim([0, zmax]);

    xlabel('\theta_A (Control)', 'Color', 'w');
    ylabel('\theta_B (Target)',  'Color', 'w');
    zlabel('|R_{int}|',          'Color', 'w');
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');

    title({ '\bf Ontic Entanglement Flow (3D)', ...
        sprintf('Time: %.2f  |  %s', t, status_msg) }, ...
        'Color', title_color, 'FontSize', 14);

    drawnow;

    %--------------------------------------------------------------
    % 4.5 Save snapshots at selected times
    %--------------------------------------------------------------
    for k = 1:numel(snapshot_times)
        if ~snapshot_taken(k) && abs(t - snapshot_times(k)) < time_tolerance
            exportgraphics(fig, snapshot_files{k}, 'Resolution', 300);
            fprintf('Saved snapshot at t = %.2f s -> %s\n', t, snapshot_files{k});
            snapshot_taken(k) = true;
        end
    end
end

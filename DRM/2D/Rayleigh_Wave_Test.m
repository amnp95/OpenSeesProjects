% Clear all variables and figures
clear all

% Material properties
nu  = 0.35;                % Poisson's ratio
Vs  = 229;                 % Shear wave velocity (m/s)
rho = 1.8;                 % Density (t/m^3)

% Loading parameters
load_duration = 9e-2;      % Duration of loading (s)
initial_time = 0;          % Initial time (s)
initial_position = -10;    % Initial position (m)
final_time = load_duration; % Final time of loading (s)

% Grid of nodes
element_size = 0.5;                        % Spacing between nodes
x = [-10:element_size:-9.5 9.5:element_size:10]; % Abscissa
y = -(0:element_size:10);                   % Ordinate
aux_x = ones(length(y), 1) * x;
aux_y = transpose(y) * ones(1, length(x));
XYnodes = [aux_x(:) aux_y(:)];

x = [-10:element_size:10]; % Abscissa
y = -(9.5:element_size:10);                   % Ordinate
aux_x = ones(length(y), 1) * x;
aux_y = transpose(y) * ones(1, length(x));

XYnodes = union(XYnodes, [aux_x(:) aux_y(:)], 'rows', 'stable');

% Amplitude and wave parameters
amplitude = 0.001;         % Amplitude of the wave
wave_type = 11;            % Type of wave (11 corresponds to a specific type)
omega = [];                % Angular frequency (empty)
dt = 1e-4;                 % Time step (s)
N = 2^12;                   % Number of time steps
t = dt * [0:N-1];          % Time vector

% Generate the time signal
A_t = amplitude * impulse_load(initial_time, load_duration, t, wave_type, omega); % [m]

% Calculate displacements and velocities using Rayleigh wave function
[u, v, ddu, ddv, Vr] = rayleigh_wave(A_t, Vs, nu, rho, dt, XYnodes, t, initial_position);

% Visualization loop
scale = 1000; % Scale factor for visualization
for i = 1:4:length(t)-1
    clf 'reset' % Clear the current figure
    plot(XYnodes(:, 1) + scale * u(:, i), XYnodes(:, 2) + scale * v(:, i), '.')
    axis([-11 11 -11 1]) % Set axis limits
    shg % Show the figure
    hold on
    grid on
    pause(dt / 10000) % Pause for visualization
end

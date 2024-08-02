clc
% define
a = 10000;    % Semi-major axis [km]
e = 0.85;     % Eccentricity
I = 0;        % Inclination [degrees]
Omega = 45;   % Right Ascension of Ascending Node [degrees]
omega = 30;   % Argument of Perigee [degrees]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define constants
G = 6.67430e-20; % Gravitational constant in km^3/kg/s^2
M = 5.972e24;    % Mass of Earth in kg
mu = G * M;      % Standard gravitational parameter in km^3/s^2

% Calculate mean motion (n) based on defined semi-major axis
n = sqrt(mu / a^3);

% Calculate orbital period (T) based on defined semi-major axis
T = 2 * pi / n; % Orbital period in seconds

fprintf('Mean Motion (n): %f rad/s\n', n);
fprintf('Orbital Period (T): %f seconds\n', T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recalculate semi-major axis based on a different T
T_new = 9600; % New orbital period in seconds
a_new = (T_new^2 * mu / (4 * pi^2))^(1/3);

fprintf('Re-calculated Semi-Major Axis (a): %f km\n', a_new);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate and plot the orbit based on Keplerian elements
M_values_calculated = linspace(0, 2*pi, 200);
x_calculated = zeros(1, length(M_values_calculated));
y_calculated = zeros(1, length(M_values_calculated));

for i = 1:length(M_values_calculated)
    M_cl = M_values_calculated(i);
    E_cl = M_cl + e*sin(M_cl) + (e^2/2)*sin(2*M_cl);
    nu_cl = 2 * atan2(sqrt(1+e) * sin(E_cl/2), sqrt(1-e) * cos(E_cl/2));
    r_cl = a * (1 - e^2) / (1 + e * cos(nu_cl));
    x_calculated(i) = r_cl * cos(nu_cl);
    y_calculated(i) = r_cl * sin(nu_cl);
end

plot(x_calculated, y_calculated);
grid on;
title('Calculated Elliptical Orbit Visualization');
xlabel('x [km]');
ylabel('y [km]');
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read orbit data from CSV
orbit_data = readmatrix('orbit_data.csv');
x_csv = orbit_data(:, 1); % X coordinates from CSV
y_csv = orbit_data(:, 2); % Y coordinates from CSV


% Plot the orbit from CSV data
figure; % Open a new figure window
plot(x_csv, y_csv);
grid on;
title('Elliptical Orbit Visualization from CSV Data');
xlabel('x [km]');
ylabel('y [km]');
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define time increment based on the orbital period
num_intervals = 24; % Number of intervals
delta_t = T / num_intervals; % Time increment

% Calculate corresponding delta M
delta_M = (n * delta_t); % in radians

% Generate M values for the entire orbit
M_values = linspace(0, 2*pi, length(orbit_data));

% Choose points corresponding to these intervals
M_intervals = 0:delta_M:(2*pi);
selected_indices = arrayfun(@(M) findNearest(M_values, M), M_intervals);

% Plot the points and lines illustrating Kepler's Second Law
figure; % Open a new figure window
hold on;
plot(x_csv, y_csv);
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    plot(x_csv(idx), y_csv(idx), 'ro'); % Plotting the point
    line([0, x_csv(idx)], [0, y_csv(idx)], 'Color', 'red'); % Line from focus to point
end

hold off;
title('Visualization of Kepler''s Second Law');
xlabel('x [km]');
ylabel('y [km]');
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the points, lines, and enhancements for Kepler's Second Law
figure; % Open a new figure window
hold on;
plot(x_csv, y_csv); % Plot the orbit

% Define a color for filling areas
fill_color = [0.9, 0.9, 0.9]; % Light gray color

for i = 1:length(selected_indices)
    idx = selected_indices(i);
    plot(x_csv(idx), y_csv(idx), 'ro'); % Plotting the point

    % Draw line from focus to point
    line([0, x_csv(idx)], [0, y_csv(idx)], 'Color', 'red'); 

    % Fill every alternate area (optional)
    if mod(i, 2) == 0 && i > 1 % Check for even index (alternate areas)
        prev_idx = selected_indices(i-1);
        x_fill = [0, x_csv(prev_idx), x_csv(idx), 0];
        y_fill = [0, y_csv(prev_idx), y_csv(idx), 0];
        fill(x_fill, y_fill, fill_color, 'LineStyle', 'none');
    end

    % Label the epoch (optional)
    text(x_csv(idx), y_csv(idx), sprintf('Epoch %d', i), 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
end

hold off;
title('Enhanced Visualization of Kepler''s Second Law with Alternate Area Filling');
xlabel('x [km]');
ylabel('y [km]');
axis equal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to find nearest M value in the array
function idx = findNearest(array, value)
    [~, idx] = min(abs(array - value));
end

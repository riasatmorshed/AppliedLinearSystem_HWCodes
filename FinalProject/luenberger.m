% Define system matrices (example, replace with your system)
A = [0 1; -2 -3];  % Example system matrix (2x2)
B = [0; 1];        % Example input matrix (2x1)
C = [1 0];         % Example output matrix (1x2)
D = 0;              % Example direct transmission matrix (scalar)

% Define observer poles (choose poles that are faster than the system poles)
alpha = 1.2;  % Factor for observer poles
sys_poles = eig(A);  % Compute system poles
observer_poles = sys_poles * alpha;  % Shift observer poles to the left

% Compute the observer gain matrix using pole placement
L = place(A', C', observer_poles)';  % Luenberger observer gain (L = place(A', C', observer_poles)')

% Define the observer system (state-space representation)
A_obsv = A - L*C;  % Observer system matrix
B_obsv = B;         % Observer uses the same input matrix
C_obsv = C;         % Observer uses the same output matrix
D_obsv = D;         % Observer uses the same direct transmission matrix

% Create the state-space model for the observer
obsv_sys = ss(A_obsv, B_obsv, C_obsv, D_obsv);

% Time vector for simulation
t = linspace(0, 10, 1000);  % Simulate from t=0 to t=10 with 1000 points

% Define initial state (example, adjust based on your system)
init_state = [1; 0];  % Initial condition for the system (2x1 column vector)

% Simulate the true system (you need the state-space model of the true system)
sys = ss(A, B, C, D);
[y_sys, t_sys] = initial(sys, init_state, t);  % True system output

% Simulate the observer system
[y_obsv, t_obsv] = initial(obsv_sys, init_state, t);  % Observer output

% Compute the observer error (difference between true system and observer output)

% Plot the observer error over time
figure;
plot(t_sys, y_sys(:,1));
hold on
plot(t_obsv, y_obsv(:,1))
title('Observer Error vs Time');
xlabel('Time (s)');
ylabel('Error');

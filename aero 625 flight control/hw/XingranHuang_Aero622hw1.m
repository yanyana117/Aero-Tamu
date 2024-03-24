%% Xingran Huang_Aero 622 Hw 1
%{ 
Purpose use block diagram to find out: 
    (a) The continuous time response to a unit step input. 
    (b) The sampled-data response to a step input for a gain of K = 1.0 and sample period of T = 1.0 seconds.
    Zoh, transfer func, state-space form, characteristic equation, eigenvalues
    (c) You must use the state-space representation of the system in your code.
    (d) You must use the state transition matrix and the discrete control distribution matrix to simulate it.
The flowchart on the next page shows how to code this.
    (e) You will also need to add a zero-order hold to the controls, to keep them constant during the sample
period.
%}
%%
clear all

% inital input:
A = [0 1; 0 -1];
B = [0 ; 1];
C = [1 0];
D = [0];

y_m = 1;
x_k = [0; 0];
C1 = [1 0; 0 1 ; 0 0];% phi
D1 = [0; 0; 1]; % gamma

T = 0.01;   % samlple period
t_initial = 0;
L0 = 0;
t_final = 10;
h = 0.01; % step size. discrete 
k = [1 0];
%% 
[Phi_1, Gamma_1] = c2d(A,B,h) %c2d: Convert model from continuous to discrete time, which discretize A and B
K = [1. 0.];
u_k = y_m - (K * x_k)

% Intial conditions for t = t_0
hold_count = 0;
int_count = 0;
max_count = T / h;
frames = (t_final - t_initial)/h;   % how many discrete n-steps do I have take
-
t = [t_initial];    % time vector: [0:h:t_final]
y = [(C1 * x_k + D1 + u_k)' ];

for i =1: frames; % Run simulation for frames
    % zoh
    if hold_count == max_count; % If statement to enforce the zero order hold
        u_k = y_m - (K * x_k); % Update control input if time T has passed since the last update
        hold_count = 0;
    % else
        % u(k+1) = u(k) % Hold the control if time T has not passed
    end
    
   % Ealucation of closed loop system
    x_k = Phi_1 * x_k + Gamma_1 * u_k; % x(k+1) = Phi*x(k) + Gamma*u(k)
    y = [y; (C1 *x_k + D1 * u_k)' ]; % y(k+1) = C1*x(k+1) + D1*u(k+1)
    t = [t; i * h];
    hold_count = hold_count + 1;
end

figure(1)
subplot(2,1,1)
plot(t, y(:, 1), 'LineWidth',  2)
axis([0. 10. 0. 2.5]);
set(gca, 'XTick', [ 0 2 4 6 8 10 ], 'XTickLabels', [ '  ' ])
set(gca, 'YTick', [ 0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 ], 'FontSize', 18)
ylabel('phi', 'FontSize', 18)
title('k = 1.0 ', 'FontSize', 18)
% grid


subplot(2,1,2)
plot(t, y(:, 3), 'LineWidth',  2)
axis([ 0. 10. -0.4 2.0 ]);
set(gca, 'XTick', [ 0 2 4 6 8 10 ], 'FontSize', 18)
set(gca, 'YTick', [ -0.5 0 0.5 1.0 1.5 2.0 ], 'FontSize', 18)
ylabel('time (sec)', 'FontSize', 18)
title('control', 'FontSize', 18)
% grid

%% 
for i=1:frames
    if hold_count == max_count
        u_k = y_m - (k*x_k);
        hold_count = 0;
    end


    x_k = Phi_1 * x_k + Gamma_1 * u_k;
    y = [y; (C1*x_k + D1*u_k)'];
    t = [t; i*h];
    hold_count = hold_count + 1;

end


figure(2)
subplot(2,1,1)
plot(t, y(:, 1), 'LineWidth',  2)
axis([0. 10. 0. 2.5]);
set(gca, 'XTick', [ 0 2 4 6 8 10 ], 'XTickLabels', [ '  ' ])
set(gca, 'YTick', [ 0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 ], 'FontSize', 18)
ylabel('phi', 'FontSize', 18)
title('k = 1.0 ', 'FontSize', 18)
grid on;

subplot(2,1,2)
plot(t, y(:, 3), 'LineWidth',  2)
axis([ 0. 10. -0.4 2.0 ]);
set(gca, 'XTick', [ 0 2 4 6 8 10 ], 'FontSize', 18)
set(gca, 'YTick', [ -0.5 0 0.5 1.0 1.5 2.0 ], 'FontSize', 18)
ylabel('time (sec)', 'FontSize', 18)
title('control', 'FontSize', 18)
grid on;

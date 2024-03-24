%% Xingran Huang_Aero 622 Hw 1
%{ 
Purpose use block diagram to find out: 
    (a) The continuous time response to a unit step input. 
    (b) The sampled-data response to a step input for a gain of K = 1.0 and sample period of T = 1.0 seconds.
    Zoh, transfer func, state-space form, characteristic equation, eigenvalues
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
C1 = [1 0; 0 1 ; 0 0];
D1 = [0; 0; 1];

T = 0.01;
t_initial = 0;
L0 = 0;
t_final = 10;
h = 0.01;
%% 
[Phi_1, Gamma_1] = c2d(A,B,h) %c2d: Convert model from continuous to discrete time .
K = [1. 0.];
u_k = y_m - (K * x_k)

% Intial conditions for t = t_0
hold_count = 0;
int_count = 0;
max_count = T / h;
frames = (t_final - t_initial)/h;
t = [t_initial];
y = [(C1 * x_k +D1 + u_k)' ];

for i =1: frames;
    % zoh
    if hold_count == max_count;
        u_k = y_m - (K * x_k);
        hold_count = 0;
    end
    
   % Ealucation of closed loop system
    x_k = Phi_1 * x_k + Gamma_1 * u_k;
    y = [y; (C1 *x_k + D1 * u_k)' ];
    t = [t; i * h];
    hold_count = hold_count + 1;
end


%%
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


% figure(2)
% rlocus(A,B,C,D)
% axis([ -2 2 -2 2 ]);
% set(gca, 'XTick', [ -2 -1 0 1 2 ], 'FontSize',28)
% set(gca, 'YTick', [ -2 -1 0 1 2 ], 'FontSize',28)
% % grid
%     

figure(2)
subplot(2,1,1)
plot(t, y(:, 1), 'LineWidth',  2)
axis([0. 10. 0. 2.5]);
set(gca, 'XTick', [ 0 2 4 6 8 10 ], 'XTickLabels', [ '  ' ])
set(gca, 'YTick', [ 0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 ], 'FontSize', 18)
ylabel('phi', 'FontSize', 18)
title('k = 1.0 ', 'FontSize', 18)
% grid

% zoh
subplot(2,1,2)
% tt = 0: 0.5 :10;
plot(stairs(t, y(:, 3)), 'LineWidth',  2)
axis([ 0. 10. -0.4 2.0 ]);
set(gca, 'XTick', [ 0 2 4 6 8 10 ], 'FontSize', 18)
set(gca, 'YTick', [ -0.5 0 0.5 1.0 1.5 2.0 ], 'FontSize', 18)
ylabel('time (sec)', 'FontSize', 18)
title('control', 'FontSize', 18)
    
 










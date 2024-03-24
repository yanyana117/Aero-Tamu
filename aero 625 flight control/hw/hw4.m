%{
- Xingran Huang
- Aero 625 homework 4
- Problem 5
%}
clc;
clear cll;
format short

m = 1;
c = -0.2;
k = 2;
T = 0.1;    % Sampling time
K = [1,1];  % Gain matrix

A = [0 1; -k/m -c/m];
B = [1 1/m]';
C = [1 0; 0 1]; % full state feedback
D = [0];

Q = eye(2,2);
R = eye(1,1);
Q(1,1) = Q(1,1)*20;
nn = [1; 0];

y_m = 0;
x_k = [1 ; 0];

C1 = [1 0; 0 1; 0 0];
D1 = [0 ; 0; 1];

T = 0.1; % sample period
to = 0;  % initial time
tf = 10; % final time
h = .01;

[PHI, GAMMA] = c2d(A,B,h); % discretize A and B
k = [1 1];
u_k = y_m - (k *x_k);

%%
hold_count = 0;
countmax = T/h;
frames = (tf-to)/h;
t = to;
y = (C1*x_k + D1*u_k)';

% Display the natural frequencies, damping ratios, time constants, and poles of sys
damp(A-B*k)

% closed loop continuous eigen values non SDR
[u,v] = eig(A-B*k);

for i=1:frames
    if hold_count == countmax
        u_k = y_m - (k*x_k);
        hold_count = 0;
    end

    x_k = PHI * x_k + GAMMA * u_k;
    y = [y; (C1*x_k + D1*u_k)'];
    t = [t; i*h];
    hold_count = hold_count + 1;

end

%%
figure(1)
subplot(2,1,1)
plot(t,y(:,1))
ylabel('x')
title('Mass Spring Damper Controller k=[1 1] T=.1')
grid on;

subplot(2,1,2)
plot(t,y(:,3))
xlabel('time')
ylabel('control')

%% initialize vars w/ to
T = .1;
hold_count = 0;
countmax = T/h;
frames = (tf-to)/h;
y_m = 0;
x_k = [1 ; 0];
k = [1 1];

[k,Qd,Rd,Nd,s,e] = lqrdjv(A,B,Q,R,nn,T);
u_k = y_m - (k *x_k);
t = to;
y = (C1*x_k + D1*u_k)';
% modal characteristics of closed loop SDR
damp(A-B*k)

% eig values SDR closed loop
[u2,v2] = eig(A-B*k);

% list needed variables for hw

disp('Q is:');
disp(Q);

disp('R is:');
disp(R);

disp('Qd is:');
disp(Qd);

disp('Rd is:');
disp(Rd);

disp('M is:');
disp(nn);

disp('k is:');
disp(k);


%%
for i=1:frames
    if hold_count == countmax
        u_k = y_m - (k*x_k);
        hold_count = 0;
    end


    x_k = PHI * x_k + GAMMA * u_k;
    y = [y; (C1*x_k + D1*u_k)'];
    t = [t; i*h];
    hold_count = hold_count + 1;

end

figure(2)
subplot(2,1,1)
plot(t,y(:,1))
ylabel('x')
title('Mass Spring Damper SDR Controller T=.1')
grid on;

subplot(2,1,2)
plot(t,y(:,3))
xlabel('time')
ylabel('control')

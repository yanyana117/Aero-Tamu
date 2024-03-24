%{
- Xingran Huang
- Aero 625 Problem 6
%}

clc
clear all
format short 


%% Initial Condition:
M1 = 0.18;
U1 = 205;   % feet/sec
H1 = 100;   % feet
alhpha_1 = 18.8;    % degree
q_bar = 50; % psf
% x_cg = 0.35*c_bar;

A = [ -0.132    0.324     -0.94     0.149   0      0       0;...
      -10.614   -1.179     1.0023   0       0      0       0; ...
       0.997    -0.00182  -0.259    0       0      0       0;...
       0         1         0.34     0       0      0       0;...
       0         0         1.0561   0       0      0       0;...
       0         0         0        0       0    -10.0000  0;...
       0         0         0        0       0      0      -10.0000]; 

B = [0.0069 0.0189; -5.935  1.203; -0.122   -0.614; 0    0; 0    0; 10 0; 0 10];
% D_matrix = [0];

% %% (6a)  
% disp(['(6a) I would like to choose " F-16A Fightcing Falcon, lat/d' ...
% ' model" as my model.' newline '     Due to our project need 2 control on it.'...
% ' And this model comes from "aircraft linear models - F22 - Ver1.ppt" page 8 email attached by Dr.Valasek.'])
% %% (6b)
% [V,D] = eig(A);  % D = eigenvalues,labda ; V = eigenvectors.
% [wn,zeta,p] = damp(A); % wn: the natural frequencies, zeta: damping ratio and p: poles of sys
% 
% disp('(6b) Eigenvalues (lambda_D) is: ')
% disp(diag(D))
% disp('(6b) Eigenvectors (V) is: ')
% disp(V)
% damp(A)
% 
% 
% %% (6c) Slection of sample period Ts
% % Take the highest frequency component from eigenvalues and equal to omega and set as equal to 2pi/T.
% wn_max = max(wn);    %
% Ts = pi/wn_max; % Sample time
% ws = 2*pi/Ts; % omega sample frequency
% 
% [phi,gamma] = c2d(A,B,Ts); % Phi is discrete A; 
% 
% fprintf('\r\n(6c) Slection of sample period Ts is %d. \r\n',Ts)
% 
% 
% %% (6d) Check system controllability
% C_new = ctrb(phi,gamma);
% unco = length(phi) - rank(C_new); 
% 
% if unco == 0 
%     fprintf('(6d) The system is controllable. Since the rank of the controllability matrix Co is equal to the number of states in the state-space model. \n')
% else 
%     fprintf('(6d) The system is uncontrollable. Rank of controllability matrix does not equal to the number of states in the state-space model. \n ')
% end
% 
% 
% %% (6d) Check system observability
% % You evaluate reachability, observability, etc using the DISCRETE versions of A, B.  Do not use the original continuous A, B. 
% % C = eye(5);
% 
% C = eye(size(phi)); % C should be identity matrix due to MIMO; C = [1 0 0 0 0] use for single input/output only.
% Ob = obsv(phi,C);  % Ob: observability matrix
% unobsv = length(phi) - rank(Ob);
% 
% if unobsv == 0 
%     fprintf('     The system is observable. Since the rank of the observability matrix ob is equal to the number of states in the state space model. \r\n' )
% else 
%     fprintf('      The system is unobservable. Rank of observability matrix does not equal to the number of states in the state-space model. \r\n ')
% end
% 
% fprintf('Update: Due on Oct.7th 5pm \r')
% 
% 
% %% Check system reachable (if A matrix is non singular)
% %  A square matrix is nonsingular iff its determinant is nonzero.
% [numRows,numCols] = size(phi);
% if det(phi) ~= 0 && numRows == numCols
%     fprintf('     The system is reachable. Since A matrix of the discrete system is nonsingular. Determinant = %d \n',det(A) )
% else
%     fprintf('      The system is unreachable. Since A matrix of the discrete system is singular. Determinant = %d \n ',det(A))
% end
% 
% 
% %% Check system constructable  
% % If O matrix is the full rank, system is constructable.
% O = obsv(phi,C);        % O: observability matrix 
% if rank(O) == size(O,2)
%     fprintf('     The system is constructable. Since observability matrix is full rank. \r\n') 
% else
%     fprintf('      The system is unconstructable. Since observability matrix is not full rank. \r\n') 
% end    
% 


% %% 
% C_identity = eye(size(phi)); 
% O_check = [C_identity; C_identity*phi; C_identity*phi^2;C_identity*phi^3 ];
% check_matrix = isequal(O,O_check)   
% % Returns logical 1 (true); returns logical 0 (false).
% 
% 
% sys = ss(A,B,C,D_matrix);
% sysd = c2d(sys,Ts,'foh')
% step(sys,'-',sysd,'--')

%% 
Q = eye(7,7); % Q matrix for adjusting weight on states
R = eye(2,2); % R matrix for adjusting weight on gains

nn = eye(7,2);
C1 = eye(11,7); % output coefficient matrix
C(8,6) = -10; 
C(9,7) = -10;
C = C1;

D1 = [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 10 0; 0 10; 1 0; 0 1]; % feed thru matrix
D = D1;

%%
y_m = 0;
x_k = [0 ; 0 ; 0 ; 10 ; 0; 0; 0];
x_k = x_k*(pi/180);
T = 0.05;     % sample period
to = 0;     % initial time
tf = 10;    % final time
h = .01;
[PHI, GAMMA] = c2d(A,B,h); % discretize A and B

[k,Qd,Rd,Nd,s,e] = lqrdjv(A,B,Q,R,nn,T);
u_k = y_m - (k *x_k);


hold_count = 0;
countmax = T/h;
frames = (tf-to)/h;
t = to;
y = (C1*x_k + D1*u_k)';

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

y = y*(180/pi);


%%
figure(1)
sgtitle('F16 Lateral Response to 10deg roll offset')

subplot(4,2,1)
plot(t,y(:,1))
title('Side Slip Angle Dynamics')
ylabel('Side Slip Angle ß')
grid on;

subplot(4,2,2)
plot(t,y(:,2))
title('Roll Rate Dynamics')
xlabel('time')
ylabel('Roll Rate p')
grid on;

subplot(4,2,3)
plot(t,y(:,3))
title('Yaw Rate Dynamics')
xlabel('time')
ylabel('Yaw Rate r')
grid on;

subplot(4,2,4)
plot(t,y(:,4))
title('Roll Angle Dynamics')
xlabel('time')
ylabel('Roll Angle ϕ')
grid on;

subplot(4,2,5)
plot(t,y(:,5))
title('Yaw Angle Dynamics')
xlabel('time')
ylabel('Yaw Angle Ψ')
grid on;

subplot(4,2,6);
plot(t,y(:,6))
title('Actuator Control')
xlabel('time')
ylabel('Actuator Angle [deg]')
grid on;
hold on;
plot(t,y(:,7))
xlabel('time')
legend("Aileron Control" , "Rudder Control")

subplot(4,2,7);
plot(t,y(:,8))
title('Actuator Control Rate')
xlabel('time')
ylabel('Actuator Rate [deg/s]')
grid on;
hold on;
plot(t,y(:,9))
xlabel('time')
legend("Aileron Control Rate" , "Rudder Control Rate")

subplot(4,2,8);
plot(t,y(:,10))
title('Actuator Commanded Value')
xlabel('time')
ylabel('Actuator Command [deg]')
grid on;
hold on;
plot(t,y(:,11))
xlabel('time')
legend("Commanded Aileron Angle" , "Commanded Rudder Angle")













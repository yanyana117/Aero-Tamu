%{
- Xingran Huang
- Aero 625 project

- Introduction From hw 3:
    (a) The Problem Statement; Include how you selected the model system.
        - Choose F-16A Fightcing Falcon, lat/d model.

    (b) Open-loop characteristics:
        eigenvalues 
        eigenvectors
        frequencies, damping rations and time constants as applicable
    (c) Selection of sample period T
    (d) Verification of controllability and observability

- Code: 
    [V,D] = eig(A): 
    Returns the eigenvectors and eigenvalues of A as symbolic matrices V and D. 
    The columns of V present eigenvectors of A. 
    The main diagonal of D present eigenvalues of A.
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

A = [ -0.132    0.324   -0.94  0.149    0;...
      -10.614 -1.179   1.0023   0     0; ...
       0.997       -0.00182  -0.259     0    0;...
       0         1       0.34       0         0;...
       0    0   1.0561  0   0]; 

B = [0.0069 0.0189; -5.935  1.203; -0.122   -0.614; 0    0; 0    0];

D_matrix = [0];

%% (6a)  
disp(['(6a) I would like to choose " F-16A Fightcing Falcon, lat/d' ...
' model" as my model.' newline '     Due to our project need 2 control on it.'...
' And this model comes from "aircraft linear models - F22 - Ver1.ppt" page 8 email attached by Dr.Valasek.'])
%% (6b)
[V,D] = eig(A);  % D = eigenvalues,labda ; V = eigenvectors.
[wn,zeta,p] = damp(A); % wn: the natural frequencies, zeta: damping ratio and p: poles of sys

disp('(6b) Eigenvalues (lambda_D) is: ')
disp(diag(D))
disp('(6b) Eigenvectors (V) is: ')
disp(V)

damp(A)


%% (6c) Slection of sample period Ts
% Take the highest frequency component from eigenvalues and equal to omega and set as equal to 2pi/T.
wn_max = max(wn);    %
Ts = pi/wn_max; % Sample time
ws = 2*pi/Ts; % omega sample frequency

[phi,gamma] = c2d(A,B,Ts); % Phi is discrete A; 

fprintf('\r\n(6c) Slection of sample period Ts is %d. \r\n',Ts)


%% (6d) Check system controllability
C_new = ctrb(phi,gamma);
unco = length(phi) - rank(C_new); 

if unco == 0 
    fprintf('(6d) The system is controllable. Since the rank of the controllability matrix Co is equal to the number of states in the state-space model. \n')
else 
    fprintf('(6d) The system is uncontrollable. Rank of controllability matrix does not equal to the number of states in the state-space model. \n ')
end


%% (6d) Check system observability
% You evaluate reachability, observability, etc using the DISCRETE versions of A, B.  Do not use the original continuous A, B. 
% C = eye(5);

C = eye(size(phi)); % C should be identity matrix due to MIMO; C = [1 0 0 0 0] use for single input/output only.
Ob = obsv(phi,C);  % Ob: observability matrix
unobsv = length(phi) - rank(Ob);

if unobsv == 0 
    fprintf('     The system is observable. Since the rank of the observability matrix ob is equal to the number of states in the state space model. \r\n' )
else 
    fprintf('      The system is unobservable. Rank of observability matrix does not equal to the number of states in the state-space model. \r\n ')
end

fprintf('Update: Due on Oct.7th 5pm \r')


%% Check system reachable (if A matrix is non singular)
%  A square matrix is nonsingular iff its determinant is nonzero.
[numRows,numCols] = size(phi);
if det(phi) ~= 0 && numRows == numCols
    fprintf('     The system is reachable. Since A matrix of the discrete system is nonsingular. Determinant = %d \n',det(A) )
else
    fprintf('      The system is unreachable. Since A matrix of the discrete system is singular. Determinant = %d \n ',det(A))
end


%% Check system constructable  
% If O matrix is the full rank, system is constructable.
O = obsv(phi,C);        % O: observability matrix 
if rank(O) == size(O,2)
    fprintf('     The system is constructable. Since observability matrix is full rank. \r\n') 
else
    fprintf('      The system is unconstructable. Since observability matrix is not full rank. \r\n') 
end    



%% 
C_identity = eye(size(phi)); 
O_check = [C_identity; C_identity*phi; C_identity*phi^2;C_identity*phi^3 ];
check_matrix = isequal(O,O_check)   % Returns logical 1 (true); returns logical 0 (false).


sys = ss(A,B,C,D_matrix);
sysd = c2d(sys,Ts,'foh')
% step(sys,'-',sysd,'--')




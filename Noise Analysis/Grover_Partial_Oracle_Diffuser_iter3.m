%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code - Noise Analysis

% Title: Performance of the Partitioned Hardware Efficient Quantum Search 
%        Algorithm with the Diagonalizable/Damping Noise Effect 

% Authors: A. Ahmadkhaniha1, Y. Mafi1, P. Kazemikhah, H. Aghababa1,
%         M.R. Kolahdouz1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Partial algorithm with Noisy diffuser and oracle %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clear
close all

ket0 = [1;0];
ket1 = [0;1];

I = eye(2);
H = 1/sqrt(2).*[1 1;1 -1];
X = [0 1;1 0];
Z = [1 0;0 -1];
CNOT = [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
CCCNOT=eye(16);
CCCNOT(15,15)=0;CCCNOT(16,15)=1;
CCCNOT(16,16)=0;CCCNOT(15,16)=1;

v = 13;
Oracle = eye(16);
Oracle(v,v) = -1;

Diffuser_Local_up = kron(kron(H,I),kron(I,I));
Diffuser_Local_up = kron(kron(X,Z),kron(I,I))*Diffuser_Local_up;
Diffuser_Local_up = kron(CNOT,kron(I,I))*Diffuser_Local_up;
Diffuser_Local_up = kron(kron(X,Z),kron(I,I))*Diffuser_Local_up;
Diffuser_Local_up = kron(kron(H,I),kron(I,I))*Diffuser_Local_up;

Diffuser_Local_down = kron(kron(I,I),kron(H,I));
Diffuser_Local_down = kron(kron(I,I),kron(X,Z))*Diffuser_Local_down;
Diffuser_Local_down = kron(kron(I,I),CNOT)*Diffuser_Local_down;
Diffuser_Local_down = kron(kron(I,I),kron(X,Z))*Diffuser_Local_down;
Diffuser_Local_down = kron(kron(I,I),kron(H,I))*Diffuser_Local_down;

Diffuser = kron(kron(H,H),kron(H,I));
Diffuser = kron(kron(X,X),kron(X,Z))*Diffuser;
Diffuser = CCCNOT*Diffuser;
Diffuser = kron(kron(X,X),kron(X,Z))*Diffuser;
Diffuser = kron(kron(H,H),kron(H,I))*Diffuser;

Psi = kron(kron(ket0,ket0),kron(ket0,ket0));    % Initial state
Psi = kron(kron(H,H),kron(H,H))*Psi;            % Hadamard operation

% Iteration #1
Psi = Oracle*Psi;                               % Oracle operation
Psi = Diffuser_Local_down*Psi;                  % Diffuser operation

% Iteration #2
Psi = Oracle*Psi;                               % Oracle operation
Psi = Diffuser*Psi;                             % Diffuser operation

% Iteration #3
Psi = Oracle*Psi;                               % Oracle operation
Psi = Diffuser_Local_up*Psi;                    % Diffuser operation


Psi_final = Psi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Noise Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms P_AD P_PD

E_AD_0 = [1 0;0 sqrt(1-P_AD)];
E_AD_1 = [0 sqrt(P_AD);0 0];

E_PD_0 = sqrt(1-P_PD).*I;
E_PD_1 = sqrt(P_PD).*[1 0;0 0];
E_PD_2 = sqrt(P_PD).*[0 0;0 1];

% Amplitude-Damping noise effect
E0 = kron(kron(E_AD_0,E_AD_0),kron(E_AD_0,E_AD_0));
E1 = kron(kron(E_AD_1,E_AD_1),kron(E_AD_1,E_AD_1));

% Phase-Damping noise effect
% E0 = kron(kron(E_PD_0,E_PD_0),kron(E_PD_0,E_PD_0));
% E1 = kron(kron(E_PD_1,E_PD_1),kron(E_PD_1,E_PD_1));
% E2 = kron(kron(E_PD_2,E_PD_2),kron(E_PD_2,E_PD_2));


Psi = kron(kron(ket0,ket0),kron(ket0,ket0));    % Initial state
rho = Psi*conj(Psi');

U1 = kron(kron(H,H),kron(H,H));                 % Hadamard operation
rho = U1*(rho*conj(U1'));

% Iteration #1
U2 = Oracle;                                    % Oracle operation
rho = U2*(rho*conj(U2'));
rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% AD Nosie operation
% rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1')+E2*rho*conj(E2'));% PD Nosie operation
U3 = Diffuser_Local_down;                         % Diffuser operation
rho = U3*(rho*conj(U3'));
rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% AD Nosie operation
% rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1')+E2*rho*conj(E2'));% PD Nosie operation

% Iteration #2
rho = U2*(rho*conj(U2'));           % Oracle operation
rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% AD Nosie operation
% rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1')+E2*rho*conj(E2'));% PD Nosie operation
U4 = Diffuser;                         % Diffuser operation
rho = U4*(rho*conj(U4'));
rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% AD Nosie operation
% rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1')+E2*rho*conj(E2'));% PD Nosie operation

% Iteration #2
rho = U2*(rho*conj(U2'));           % Oracle operation
rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% AD Nosie operation
% rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1')+E2*rho*conj(E2'));% PD Nosie operation
U4 = Diffuser_Local_down;                         % Diffuser operation
rho = U4*(rho*conj(U4'));
rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% AD Nosie operation
% rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1')+E2*rho*conj(E2'));% PD Nosie operation

% Fidelity calculation
Fidelity = simplify(conj(Psi_final')*rho*Psi_final);



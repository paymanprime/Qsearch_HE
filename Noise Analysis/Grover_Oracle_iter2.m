%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code - Noise Analysis

% Title: Performance of the Partitioned Hardware Efficient Quantum Search 
%        Algorithm with the Diagonalizable/Damping Noise Effect 

% Authors: A. Ahmadkhaniha1, Y. Mafi1, P. Kazemikhah, H. Aghababa1,
%         M.R. Kolahdouz1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Grover algorithm with Noisy oracle %%%%%%%%%%%%%%%%%%%%%%
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
CCCNOT=eye(16);
CCCNOT(15,15)=0;CCCNOT(16,15)=1;
CCCNOT(16,16)=0;CCCNOT(15,16)=1;

v = 13;  % The desired state
Oracle = eye(16);
Oracle(v,v) = -1;

Diffuser = kron(kron(H,H),kron(H,I));
Diffuser = kron(kron(X,X),kron(X,Z))*Diffuser;
Diffuser = CCCNOT*Diffuser;
Diffuser = kron(kron(X,X),kron(X,Z))*Diffuser;
Diffuser = kron(kron(H,H),kron(H,I))*Diffuser;

Psi = kron(kron(ket0,ket0),kron(ket0,ket0));    % Initial state
Psi = kron(kron(H,H),kron(H,H))*Psi;            % Hadamard operation
figure
subplot(1,2,1)
bar(Psi)
subplot(1,2,2)
bar(Psi.^2)
ylim([0 1])

Psi = Oracle*Psi;                               % Oracle operation
% figure
% subplot(1,2,1)
% bar(Psi)
% subplot(1,2,2)
% bar(Psi.^2)
% ylim([0 1])

Psi = Diffuser*Psi;                             % Diffuser operation
figure
subplot(1,2,1)
bar(Psi)
subplot(1,2,2)
bar(Psi.^2)
ylim([0 1])
Iteration #2
Psi = Oracle*Psi;                               % Oracle operation
Psi = Diffuser*Psi;                             % Diffuser operation
% figure
% bar(Psi.^2)
% ylim([0 1])
Psi_final = Psi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Noise Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms P_AD P_PD

% Kraus operators of amplitude-damping noise
E_AD_0 = [1 0;0 sqrt(1-P_AD)];
E_AD_1 = [0 sqrt(P_AD);0 0];

% Kraus operators of phase-damping noise
E_PD_0 = sqrt(1-P_PD).*I;
E_PD_1 = sqrt(P_PD).*[1 0;0 0];
E_PD_2 = sqrt(P_PD).*[0 0;0 1];

Psi = kron(kron(ket0,ket0),kron(ket0,ket0));    % Initial state
rho = Psi*conj(Psi');

U1 = kron(kron(H,H),kron(H,H));                 % Hadamard operation
rho = U1*(rho*conj(U1'));

U2 = Oracle;                                    % Oracle operation
rho = U2*(rho*conj(U2'));

% Amplitude-Damping noise effect
E0 = kron(kron(E_AD_0,E_AD_0),kron(E_AD_0,E_AD_0));
E1 = kron(kron(E_AD_1,E_AD_1),kron(E_AD_1,E_AD_1));

% Phase-Damping noise effect
% E0 = kron(kron(E_PD_0,E_PD_0),kron(E_PD_0,E_PD_0));
% E1 = kron(kron(E_PD_1,E_PD_1),kron(E_PD_1,E_PD_1));
% E2 = kron(kron(E_PD_2,E_PD_2),kron(E_PD_2,E_PD_2));

rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% Nosie operation

U3 = Diffuser;                                  % Diffuser operation
rho = U3*(rho*conj(U3'));

rho = U2*(rho*conj(U2'));                       % Oracle operation
rho = simplify(E0*rho*conj(E0')+E1*rho*conj(E1'));% Nosie operation
rho = U3*(rho*conj(U3'));                       % Diffuser operation

% Fidelity calculation
Fidelity = simplify(conj(Psi_final')*rho*Psi_final);




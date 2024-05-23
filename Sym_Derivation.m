clc, clear all, close all

syms q1 q2 q3 
syms dq1 dq2 dq3 
syms Ixx Iyy Izz Ixy Ixz Iyz

Rzyaw = [cos(q1) sin(q1) 0;-sin(q1) cos(q1) 0;0 0 1];
Rypitch = [cos(q2) 0 -sin(q2);0 1 0;sin(q2) 0 cos(q2)];
Rxroll = [1 0 0;0 cos(q3) sin(q3);0 -sin(q3) cos(q3)];

R = Rxroll * Rypitch *Rzyaw;
ws = R*[0;0;dq1] + Rxroll * Rypitch *[0;dq2;0] + Rxroll * [dq3;0;0];
I = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz];
Hs = I*ws;
Ts = 1/2 *ws.'*Hs;
Vs = 0;

% Ts = simplify(Ts)

% normHs = simplify(Hs(1)^2+Hs(2)^2+Hs(3)^2)

W = [0 sin(q3)/cos(q2) cos(q3)/cos(q2);0 cos(q3) -sin(q3);1 sin(q3)*tan(q2) cos(q3)*tan(q2)];

invW = simplify(inv(W))
% Lagrangian
L = Ts - Vs;

q = [q1; q2; q3];
dq = [dq1; dq2; dq3];

dL_dq = jacobian(L,q);
dL_ddq = jacobian(L,dq);

% Mass Matrix
M = jacobian (dL_ddq , dq);
B = jacobian (dL_ddq , q)*dq - dL_dq.';
Q = [0;0;0]; %No Active Force

M = simplify(M)
B = simplify(B)

% simplify(W.'*M*W)

B1 = simplify(cross(ws,Hs))

Hs_Global = simplify(R.'*Hs)
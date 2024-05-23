function dz = Sat_Dyn_MLag_Quat(t,z)
dz = z;
% z(1:4) = z(1:4)/norm(z(1:4));

global Ixx Iyy Izz Ixy Ixz Iyz

u1 = z(5); u2 = z(6); u3 = z(7);

M1 = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz];

W = [0 -u1 -u2 -u3;u1 0 u3 -u2;u2 -u3 0 u1;u3 u2 -u1 0];

B1 = [0 -u3 u2;u3 0 -u1;-u2 u1 0]*M1*z(5:7);

dz(1:4) = 1/2*W*z(1:4);
dz(5:7) = - M1\B1;
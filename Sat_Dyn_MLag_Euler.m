function dz = Sat_Dyn_MLag_Euler(t,z)
dz = z;

global Ixx Iyy Izz Ixy Ixz Iyz

q1 = z(1); q2 = z(2); q3 = z(3);
u1 = z(4); u2 = z(5); u3 = z(6);

M1 = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz];

W = [0 sin(q3)/cos(q2) cos(q3)/cos(q2);0 cos(q3) -sin(q3);1 sin(q3)*tan(q2) cos(q3)*tan(q2)];

B1 = [0 -u3 u2;u3 0 -u1;-u2 u1 0]*M1*z(4:6);

dz(1:3) = W*z(4:6);
dz(4:6) = - M1\B1;
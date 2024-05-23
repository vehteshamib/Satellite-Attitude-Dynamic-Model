function dz = Sat_Dyn_Lag_Euler(t,z)
dz = z;

global Ixx Iyy Izz Ixy Ixz Iyz

q1 = z(1); q2 = z(2); q3 = z(3);
dq1 = z(4); dq2 = z(5); dq3 = z(6);

M = [ Izz*cos(q2)^2*cos(q3)^2 - 2*Iyz*cos(q2)^2*cos(q3)*sin(q3) + Iyy*cos(q2)^2*sin(q3)^2 + 2*Ixz*cos(q2)*cos(q3)*sin(q2) + 2*Ixy*cos(q2)*sin(q2)*sin(q3) + Ixx*sin(q2)^2, Iyz*cos(q2) + Ixy*cos(q3)*sin(q2) - Ixz*sin(q2)*sin(q3) - 2*Iyz*cos(q2)*cos(q3)^2 + Iyy*cos(q2)*cos(q3)*sin(q3) - Izz*cos(q2)*cos(q3)*sin(q3), - Ixx*sin(q2) - Ixz*cos(q2)*cos(q3) - Ixy*cos(q2)*sin(q3)
    Iyz*cos(q2) + Ixy*cos(q3)*sin(q2) - Ixz*sin(q2)*sin(q3) - 2*Iyz*cos(q2)*cos(q3)^2 + Iyy*cos(q2)*cos(q3)*sin(q3) - Izz*cos(q2)*cos(q3)*sin(q3),                                                                                         Iyy*cos(q3)^2 + 2*Iyz*cos(q3)*sin(q3) + Izz*sin(q3)^2,                                 Ixz*sin(q3) - Ixy*cos(q3)
    - Ixx*sin(q2) - Ixz*cos(q2)*cos(q3) - Ixy*cos(q2)*sin(q3),                                                                                                                     Ixz*sin(q3) - Ixy*cos(q3),                                                       Ixx];

B = [Ixy*dq2^2*cos(q2)*cos(q3) - Iyz*dq2^2*sin(q2) - Ixy*dq3^2*cos(q2)*cos(q3) - Ixz*dq2^2*cos(q2)*sin(q3) + Ixz*dq3^2*cos(q2)*sin(q3) - Ixx*dq2*dq3*cos(q2) - 2*Ixz*dq1*dq2*cos(q3) - Iyy*dq2*dq3*cos(q2) + Izz*dq2*dq3*cos(q2) - 2*Ixy*dq1*dq2*sin(q3) + 2*Iyz*dq2^2*cos(q3)^2*sin(q2) + 2*Iyz*dq1*dq3*cos(q2)^2 - Iyy*dq2^2*cos(q3)*sin(q2)*sin(q3) + Izz*dq2^2*cos(q3)*sin(q2)*sin(q3) + 2*Ixx*dq1*dq2*cos(q2)*sin(q2) - 2*Iyy*dq1*dq2*cos(q2)*sin(q2) + 4*Ixz*dq1*dq2*cos(q2)^2*cos(q3) + 2*Iyy*dq2*dq3*cos(q2)*cos(q3)^2 - 2*Izz*dq2*dq3*cos(q2)*cos(q3)^2 + 4*Ixy*dq1*dq2*cos(q2)^2*sin(q3) - 4*Iyz*dq1*dq3*cos(q2)^2*cos(q3)^2 - 2*Ixz*dq1*dq3*cos(q2)*sin(q2)*sin(q3) + 2*Iyy*dq1*dq2*cos(q2)*cos(q3)^2*sin(q2) + 2*Iyy*dq1*dq3*cos(q2)^2*cos(q3)*sin(q3) - 2*Izz*dq1*dq2*cos(q2)*cos(q3)^2*sin(q2) - 2*Izz*dq1*dq3*cos(q2)^2*cos(q3)*sin(q3) + 2*Ixy*dq1*dq3*cos(q2)*cos(q3)*sin(q2) + 4*Iyz*dq2*dq3*cos(q2)*cos(q3)*sin(q3) + 4*Iyz*dq1*dq2*cos(q2)*cos(q3)*sin(q2)*sin(q3)
    Ixz*dq1^2*cos(q3) + Ixz*dq3^2*cos(q3) + Ixy*dq1^2*sin(q3) + Ixy*dq3^2*sin(q3) - 2*Iyz*dq2*dq3 - Ixx*dq1^2*cos(q2)*sin(q2) + Iyy*dq1^2*cos(q2)*sin(q2) + Ixx*dq1*dq3*cos(q2) - Iyy*dq1*dq3*cos(q2) + Izz*dq1*dq3*cos(q2) - 2*Ixz*dq1^2*cos(q2)^2*cos(q3) - 2*Ixy*dq1^2*cos(q2)^2*sin(q3) + 4*Iyz*dq2*dq3*cos(q3)^2 - 2*Ixz*dq1*dq3*cos(q3)*sin(q2) - 2*Iyy*dq2*dq3*cos(q3)*sin(q3) + 2*Izz*dq2*dq3*cos(q3)*sin(q3) - 2*Ixy*dq1*dq3*sin(q2)*sin(q3) - Iyy*dq1^2*cos(q2)*cos(q3)^2*sin(q2) + Izz*dq1^2*cos(q2)*cos(q3)^2*sin(q2) + 2*Iyy*dq1*dq3*cos(q2)*cos(q3)^2 - 2*Izz*dq1*dq3*cos(q2)*cos(q3)^2 - 2*Iyz*dq1^2*cos(q2)*cos(q3)*sin(q2)*sin(q3) + 4*Iyz*dq1*dq3*cos(q2)*cos(q3)*sin(q3)
    Iyz*dq2^2 - Iyz*dq1^2*cos(q2)^2 - 2*Iyz*dq2^2*cos(q3)^2 + Iyy*dq2^2*cos(q3)*sin(q3) - Izz*dq2^2*cos(q3)*sin(q3) - Ixx*dq1*dq2*cos(q2) + Iyy*dq1*dq2*cos(q2) - Izz*dq1*dq2*cos(q2) + 2*Iyz*dq1^2*cos(q2)^2*cos(q3)^2 + Ixz*dq1^2*cos(q2)*sin(q2)*sin(q3) + 2*Ixz*dq1*dq2*cos(q3)*sin(q2) + 2*Ixy*dq1*dq2*sin(q2)*sin(q3) - Iyy*dq1^2*cos(q2)^2*cos(q3)*sin(q3) + Izz*dq1^2*cos(q2)^2*cos(q3)*sin(q3) - 2*Iyy*dq1*dq2*cos(q2)*cos(q3)^2 + 2*Izz*dq1*dq2*cos(q2)*cos(q3)^2 - Ixy*dq1^2*cos(q2)*cos(q3)*sin(q2) - 4*Iyz*dq1*dq2*cos(q2)*cos(q3)*sin(q3)];

dz(1:3) = z(4:6);
dz(4:6) = - M\B;
clc, clear all, close all

global Ixx Iyy Izz Ixy Ixz Iyz

Ixx = 1000; Iyy = 800; Izz = 700; Ixy = 5; Ixz = 8; Iyz = 3;
I = [Ixx -Ixy -Ixz;-Ixy Iyy -Iyz;-Ixz -Iyz Izz];

z0 = [0;0;0;0.3;-0.5;0.4];
timespan = [0:0.01:10];
options=odeset('maxstep',10^-2);
% options=[];

t0 = clock;
[t,z] = ode45(@Sat_Dyn_Lag_Euler, timespan ,z0,options);
CPU_TIME_Lag = etime(clock,t0)

q1 = z(:,1); q2 = z(:,2); q3 = z(:,3);
dq1 = z(:,4); dq2 = z(:,5); dq3 = z(:,6);

Len = length(t);

Ts_Lag = t;
Hsx_Lag = t;
Hsy_Lag = t;
Hsz_Lag = t;

for i = 1:Len
    ws = [dq3(i) - dq1(i)*sin(q2(i))
        dq2(i)*cos(q3(i)) + dq1(i)*cos(q2(i))*sin(q3(i))
        dq1(i)*cos(q2(i))*cos(q3(i)) - dq2(i)*sin(q3(i))];
    
    R = angle2dcm(q1(i),q2(i),q3(i));    
    
    HsG = R'*I*ws;
    
    Ts_Lag(i) = 1/2*ws'*I*ws;
    Hsx_Lag(i) = HsG(1);
    Hsy_Lag(i) = HsG(2);
    Hsz_Lag(i) = HsG(3);    
end

% ------------------ Modified Lagrange ----------------

z0M = z0;
invW0 = [        -sin(z0(2)),        0, 1
    cos(z0(2))*sin(z0(3)),  cos(z0(3)), 0
    cos(z0(2))*cos(z0(3)), -sin(z0(3)), 0];
z0M(4:6) = invW0*z0(4:6);

t0 = clock;
[t,z] = ode45(@Sat_Dyn_MLag_Euler, timespan,z0M,options);
CPU_TIME_MLag = etime(clock,t0)

q1M = z(:,1); q2M = z(:,2); q3M = z(:,3);
u1 = z(:,4); u2 = z(:,5); u3 = z(:,6);

Ts_MLag = t;
Hsx_MLag = t;
Hsy_MLag = t;
Hsz_MLag = t;

for i = 1:Len
    ws = [u1(i);u2(i);u3(i)];
    
    R = angle2dcm(q1(i),q2(i),q3(i));      
    
    HsG = R'*I*ws;    
    
    Ts_MLag(i) = 1/2*ws'*I*ws;
    Hsx_MLag(i) = HsG(1);
    Hsy_MLag(i) = HsG(2);
    Hsz_MLag(i) = HsG(3); 
end

% --------------- Modified Lagrange + Quaternions ----------------
z0MQ = zeros(7,1);
z0MQ(1:4) = angle2quat(z0(1),z0(2),z0(3),'ZYX');
z0MQ(5:7) = z0M(4:6);

t0 = clock;
[t,z] = ode45(@Sat_Dyn_MLag_Quat, timespan,z0MQ,options);
CPU_TIME_MLag_Q = etime(clock,t0)

e0 = z(:,1); e1 = z(:,2); e2 = z(:,3); e3 = z(:,4);
u1 = z(:,5); u2 = z(:,6); u3 = z(:,7);

q1MQ = t;  q2MQ = t;  q3MQ = t;

for i = 1:Len
    [q1MQ(i),q2MQ(i),q3MQ(i)] = quat2angle([e0(i),e1(i),e2(i),e3(i)]);
end

Ts_MLagQ = t;
Hsx_MLagQ = t;
Hsy_MLagQ = t;
Hsz_MLagQ = t;
for i = 1:Len
    ws = [u1(i);u2(i);u3(i)];
    R = quat2dcm([e0(i),e1(i),e2(i),e3(i)]);
    
    HsG = R'*I*ws;
    
    Ts_MLagQ(i) = 1/2*ws'*I*ws;
    Hsx_MLagQ(i) = HsG(1);
    Hsy_MLagQ(i) = HsG(2);
    Hsz_MLagQ(i) = HsG(3); 
end

% -------------- Double Checking CPU TIME for Lagrange Method ------------

t0 = clock;
[t,z] = ode45(@Sat_Dyn_Lag_Euler, timespan,z0,options);
CPU_TIME_Lag = etime(clock,t0)


% ------------------------------- Plot Section ----------------------------
figure
hold on
plot(t,(Ts_Lag-Ts_Lag(1))/Ts_Lag(1)*100,'r-','linewidth',4)
plot(t,(Ts_MLag-Ts_MLag(1))/Ts_MLag(1)*100,'b-','linewidth',4)
plot(t,(Ts_MLagQ-Ts_MLagQ(1))/Ts_MLagQ(1)*100,'g--','linewidth',4)

legend('Lagrange','Modified Lagrange','Modified Lagrange + Quaternions')
set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',25,'fontweight','bold');
ylabel('Kinetic Energy Change(%)','fontsize',25,'fontweight','bold');


figure
hold on
plot(t,(Hsx_Lag-Hsx_Lag(1))/Hsx_Lag(1)*100,'r-','linewidth',4)
plot(t,(Hsx_MLag-Hsx_MLag(1))/Hsx_MLag(1)*100,'b-','linewidth',4)
plot(t,(Hsx_MLagQ-Hsx_MLagQ(1))/Hsx_MLagQ(1)*100,'g--','linewidth',4)

legend('Lagrange','Modified Lagrange','Modified Lagrange + Quaternions')

set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',15,'fontweight','bold');
ylabel('Change of Angular Momentum in x Direction(%)','fontsize',20,'fontweight','bold');

figure
hold on
plot(t,(Hsy_Lag-Hsy_Lag(1))/Hsy_Lag(1)*100,'r-','linewidth',4)
plot(t,(Hsy_MLag-Hsy_MLag(1))/Hsy_MLag(1)*100,'b-','linewidth',4)
plot(t,(Hsy_MLagQ-Hsy_MLagQ(1))/Hsy_MLagQ(1)*100,'g--','linewidth',4)

legend('Lagrange','Modified Lagrange','Modified Lagrange + Quaternions')

set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',15,'fontweight','bold');
ylabel('Change of Angular Momentum in y Direction(%)','fontsize',20,'fontweight','bold');

figure
hold on
plot(t,(Hsz_Lag-Hsz_Lag(1))/Hsz_Lag(1)*100,'r-','linewidth',4)
plot(t,(Hsz_MLag-Hsz_MLag(1))/Hsz_MLag(1)*100,'b-','linewidth',4)
plot(t,(Hsz_MLagQ-Hsz_MLagQ(1))/Hsz_MLagQ(1)*100,'g--','linewidth',4)

legend('Lagrange','Modified Lagrange','Modified Lagrange + Quaternions')

set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',15,'fontweight','bold');
ylabel('Change of Angular Momentum in z Direction(%)','fontsize',20,'fontweight','bold');


figure
hold on
plot(t,q1,'m-','linewidth',4)
plot(t,q2,'c-','linewidth',4)
plot(t,q3,'y-','linewidth',4)

plot(t,q1M,'r--','linewidth',4)
plot(t,q2M,'b--','linewidth',4)
plot(t,q3M,'k--','linewidth',4)

plot(t(1:10:end),q1MQ(1:10:end),'rs','linewidth',1)
plot(t(1:10:end),q2MQ(1:10:end),'bs','linewidth',1)
plot(t(1:10:end),q3MQ(1:10:end),'ks','linewidth',1)

legend('\psi','\theta','\phi')
set(gca,'fontsize',18,'fontweight','bold');
xlabel('Time (s)','fontsize',20,'fontweight','bold');
ylabel('Satellite Euler Angles (rad)','fontsize',20,'fontweight','bold');

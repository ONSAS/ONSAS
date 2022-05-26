function qelem = WOMV1(vprvect1, vprvect2, Udotdottp1k1, Udotdottp1k2,tl1, tl2, D, tnp1, dt)
% Computes the value of q for the element
% vpr1, vpr2: relative velocities at nodes 1 and 2 (vpr = Ucos(theta0))
% Explain Udotdott1, Udotdott2,tl1, tl2, D, CL0, tn, dt
global qvect; % to be replaced by a xls file
n = uint16(tnp1/dt - 1 % Current time
if n == 0 % First call
    %create xlsfile or txt file
    q01 = 2; q02 = 2; % Should be random number of magnitude O(0.001)
    %qnp1 = [q01;0;q02; 0] ; % qvect = [q0n1 q1n1 ...;dq0n1 dq1n1 ...;q0n2 q1n2 ...;dq0n2 dq1n2 ...;] 
    qvect = [q01;0;q02;0] ; % qvect = [q0n1 q1n1 ...;dq0n1 dq1n1 ...;q0n2 q1n2 ...;dq0n2 dq1n2 ...;] 
    qnp11 = q01;
    qnp12 = q02;
    qelem = (q01+q02)/2;
end
    vpr1 = norm(vprvect1);
    vpr2 = norm(vprvect2);
    %project Udotdottp1k1 and Udotdottp1k2 on tl1 and tl2 to have the
    %transverse acceleration
    ddY1 = dot(Udotdottp1k1, tl1); % Ydotdott of node 1 at time n
    ddY2 = dot(Udotdottp1k2, tl2); % Ydotdott of node 2 at time n
    % Calls computeq at both nodes
    [qnp11 dqnp11] = computeq(0, ddY1, D, tnp1, dt, vpr1);
    [qnp12 dqnp12] = computeq(1, ddY2, D, tnp1, dt, vpr2);
    % Updating qvect
    qvect(1:4, n+2) = [qnp11 dqnp11 qnp12 dqnp12];
    % Computing averaged qelem
    qelem = (qnp11+qnp12)/2

%function computing q at node i
function [qnp1 dqnp1] = computeq(i, ddYi, D, tnp1, dt, vpri)% at node i= 1,2
    %time increments
    N= 20; % Number of steps
    h= dt/N; % Time step
    t = tnp1 - dt:h:tnp1; % Interval on which ode45 solves the VdP equation
    % Initial parameters
    qvect
    qn = qvect(i+1,end);
    dqn = qvect(i+2,end);
    % WOM constants from Facchinetti et al
    A = 12; epsilon = 0.3; St = 0.16;
    % VdP oscillator constants
    cq = epsilon*(2*pi*St*vpri/D);
    kq = (2*pi*St*vpri/D)^2;
    [t,qode] = ode45(@(t, q) funcvanderpol(t,q, cq, kq, ddYi, D, A),t,[qn,dqn]);
    qnp1 = qode(end,1); 
    dqnp1 = qode(end,2);
end

function dq=funcvanderpol(t,q, c, k, ddY, D, A)
    qdot = q(2);
    qdotdot = c*(1-q(1)^2)*q(2)-k*q(1)+A*ddY/D; 
    dq = [qdot; qdotdot];        
end
end


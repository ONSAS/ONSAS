% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

function qelem = WOMV4(vprvect1, vprvect2, Udotdottp1k1, Udotdottp1k2,tl1, tl2, D, tnp1, dt,Kelem)
% Computes the value of q for the element Kelem with random initial
% conditions for q
% vpr1, vpr2: relative velocities at nodes 1 and 2 (vpr = Ucos(theta0))
% D: diameter
% tnp1: tn+dt
% Kelem: element id
global qvect; 
n = uint16(tnp1/dt - 1); % Current time
K = Kelem-1;
if n == 0 % First call
    qelem = (2*rand-1)*0.001;
    % Initialize qvect 
    qvect(1+K*2:2+K*2,1) = [qelem ; 0] ;  % qvect = [q0elem q1elem ...; dq0elem dq1elem ...;
else
    vpr1 = norm(vprvect1);
    vpr2 = norm(vprvect2);
    vprelem = (vpr1 + vpr2)/2;
    %project Udotdottp1k1 and Udotdottp1k2 on tl1 and tl2 to have the
    %transverse acceleration
    ddY1 = dot(Udotdottp1k1, tl1); % Ydotdott of node 1 at time n
    ddY2 = dot(Udotdottp1k2, tl2); % Ydotdott of node 2 at time n
    ddYelem = (ddY1 + ddY2)/2;
    % Calls computeq for the element
    qn = qvect(1+K*2,n);
    dqn = qvect(2+K*2,n);
    %[qnp1elem dqnp1elem ] = computeq(ddY2, D, tnp1, dt, vprelem, qn, dqn);
    [qnp1elem dqnp1elem ] = computeq(ddYelem, D, tnp1, dt, vprelem, qn, dqn);
    qelem = qnp1elem;
    % Updating qvect
    qvect(1+K*2:2+K*2, n+1) = [qnp1elem dqnp1elem];
end
%function computing q at node i
function [qnp1 dqnp1] = computeq(ddYelem, D, tnp1, dt, vprelem, qn, dqn)% at node i= 1,2
    %time increments
    N= 2; % Number of steps
    h= dt/N; % Time step
    t = tnp1 - dt:h:tnp1; % Interval on which ode45 solves the VdP equation
    % WOM constants from Facchinetti et al
    A = 12; epsilon = 0.3; St = 0.2;
    % VdP oscillator constants
    omegaf = 2*pi*St*vprelem/D; % Shedding pulsation
    cq = epsilon*omegaf;
    kq = omegaf^2;
    [t,qode] = ode45(@(t, q) funcvanderpol(t,q, cq, kq, ddYelem, D, A),t,[qn,dqn]);
    qnp1 = qode(end,1); 
    dqnp1 = qode(end,2);
end

function dq=funcvanderpol(t,qn, cq, kq, ddYnp1, D, A)
    qdot = qn(2);
    qdotdot = cq*(1-qn(1)^2)*qn(2)-kq*qn(1)+A*ddYnp1/D; 
    dq = [qdot; qdotdot];        
end
end
% Copyright 2023, ONSAS Authors (see documentation)
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
%
function pelem = WOM_IL(Uprvect1, Uprvect2, Udotdottp1k1, Udotdottp1k2, D, tnp1, dt,Kelem)
% Computes the value of q for the element Kelem with random initial
% conditions for q
% vpr1, vpr2: relative velocities at nodes 1 and 2 (vpr = Ucos(theta0))
% D: diameter
% tnp1: tn+dt
% Kelem: element id
% Cdmean = 1.2; Cdoscillating = 0.2; Ax = 12, epsilonx = 0.3
global pvect; 
n = uint16(tnp1/dt - 1); % Current time
K = Kelem-1;
if n == 0 % First call
    pelem = pvect(1+K*2,1) ;  % pvect = [p0elem p1elem ...; dp0elem dp1elem ...;
else
    Upr1 = norm(Uprvect1);
    Upr2 = norm(Uprvect2);
    Uelem = (Upr1 + Upr2)/2;
    %project Udotdottp1k1 and Udotdottp1k2 on tl1 and tl2 to have the
    tdrag = Uprvect1/Upr1; % constant along X
    %transverse acceleration
    ddX1 = dot(Udotdottp1k1, tdrag); % Xdotdott of node 1 at time n
    ddX2 = dot(Udotdottp1k2, tdrag); % Xdotdott of node 2 at time n
    ddXelem = (ddX1 + ddX2)/2;
    % Calls computep for the element
    pn = pvect(1+K*2,n);
    dpn = pvect(2+K*2,n);
    [pnp1elem dpnp1elem ] = computep(ddXelem, D, tnp1, dt, Uelem, pn, dpn);
    pelem = pnp1elem;
    % Updating pvect
    pvect(1+K*2:2+K*2, n+1) = [pnp1elem dpnp1elem];
end
%function computing p at node i
function [pnp1 dpnp1] = computep(ddXelem, D, tnp1, dt, U, pn, dpn)% at node i= 1,2
    %time increments
    N= 2; % Number of steps
    h= dt/N; % Time step
    t = tnp1 - dt:h:tnp1; % Interval on which ode45 solves the VdP epuation
    % WOM constants from Facchinetti et al
    %Ax = 12; epsilonx = 0.15; St = 0.2;
    % WOM constants from Wang et al 2018
    Ax = 96; epsilonx = 0.02; 
    St = 0.2;
    %St = 0.17;
    % VdP oscillator constants
    omegaf = 2*pi*St*U/D; % Shedding pulsation 
    cp = 2*epsilonx*omegaf;
    kp = 4*omegaf^2;
    [t,pode] = ode45(@(t, p) funcvanderpol(t,p, cp, kp, ddXelem, D, Ax),t,[pn,dpn]);
    pnp1 = pode(end,1); 
    dpnp1 = pode(end,2);
end

function dp=funcvanderpol(t,pn, cp, kp, ddXnp1, D, A)
    pdot = pn(2);
    pdotdot = cp*(1-pn(1)^2)*pn(2)-kp*pn(1)+A*ddXnp1/D; 
    dp = [pdot; pdotdot];
end
end

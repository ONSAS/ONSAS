% Copyright 2025, ONSAS Authors (see documentation)
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
% md This functions computes manually the aerodynamic loads submitted to the hole beam
function f = myLinearAeroLA(t, UsCell)
  numElements = 10;
  % geometric parameers
  l = 20;
  d = .5;
  % fluid parameters
  rhoA = 1.225;
  % inititialize external force
  f = zeros((numElements + 1) * 6, 1);
  % read the wnd velocity
  windVel = feval('windVelLA', 0, t);
  % the angle of incidence is
  betaRel = acos(dot([0 1 0], [0 0 1]));
  % the drag, lift and moment coefficients are:
  c_d = feval('dragCoefFunctionLA', betaRel);
  c_l = feval('liftCoefFunctionLA', betaRel);
  c_m = feval('momentCoefFunctionLA', -betaRel);
  % md Then the dynamic pressures $q_0$ defined above are expressed such that:
  q = 1 / 2 * rhoA * (windVel(3)^2 + windVel(2)^2);
  % md next the loads per unit of length are
  qz = q * c_d * d;
  qy = q * c_l * d;
  qm = q * c_m * d;

  % add forces into f vector
  for elem = 1:numElements
    dofsElem = (elem - 1) * 6 + 1:(elem - 1) * 6 + 12;
    lelem = l / numElements;
    % torsional laods
    Tx = qm * lelem / 2;
    dofsMx = dofsElem(2:6:end);
    % y Loads
    Fy = -qy * lelem / 2;
    dofsFY = dofsElem(3:6:end);
    Mz = qy * lelem^2 / 12;
    dofsMz = dofsElem(6:6:end);
    % z Loads
    Fz = qz * lelem / 2;
    dofsFz = dofsElem(5:6:end);
    My = qz * lelem^2 / 12;
    dofsMy = dofsElem(4:6:end);

    % add into f vectos
    f(dofsMx) = f(dofsMx) + Tx;
    f(dofsMy) = f(dofsMy) + [-My; My];
    f(dofsMz) = f(dofsMz) + [-Mz; Mz];
    f(dofsFY) = f(dofsFY) + Fy;
    f(dofsFz) = f(dofsFz) + Fz;
  end

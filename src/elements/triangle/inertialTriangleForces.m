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
% md Function for computation of nodal forces and tangent stiffness matrix for 2D 3 nodes triangle element with elastic-linear and linear isotropic hardening with von Mises flow rule models. The dofs are assumed to be on x-y.
% md
function [Fmase, Mmase ] = inertialTriangleForces(elemCoords, thickness, density, dotdotdispsElem)

  x = elemCoords(1:3:end)';
  y = elemCoords(2:3:end)';

  A = 0.5 * det([ones(1, 3); x'; y']); % element area

  assert(A >= 0, 'Element with negative area, check connectivity.');

  % md### compute inertial loads and mass matrix
  Mmase = sparse(9,9);
  Mmase([1 4 7],[1 4 7]) = density * A / 3.0 * speye(3, 3);
  Fmase = Mmase * dotdotdispsElem;

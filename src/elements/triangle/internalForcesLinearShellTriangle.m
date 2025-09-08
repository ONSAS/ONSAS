% Copyright 2024, ONSAS Authors (see documentation)
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
% Implementation of a triangular finite element with 6 dfos (3 translations and 3 rotations) per node for the analysis of linear elastic isotropic shells with constant thickness.
% The element is formed by the superposition of a plate element (DKT) and a plane stress element (CST) with
% with addition to artificial drilling (rotation about the axis normal to the element plane) stiffness.
%
function [fs, ks, fintLocCoord, Kefora, Kb] = internalForcesLinearShellTriangle(elemCoords, elemDisps, modelName, modelParams, thickness)

  % material and geometric parameters
  young_modulus = modelParams(1);
  poisson_ratio = modelParams(2);
  h = thickness;

  elemCoords = elemCoords';

  r1g = elemCoords(1:3);
  r2g = elemCoords(4:6);
  r3g = elemCoords(7:9);

  % elemDisps
  Ug = switchToTypeIndexing(elemDisps);

  [T, x02, x03, y03] = edgeLocalAxisShellTriangle(r1g, r2g, r3g);
  Te = blkdiag(T, T, T, T, T, T);
  Ul = Te * Ug;
  %
  global flag_OPT
  %
  [Kl, fintLocCoord] = localShellTriangle(x02, x03, y03, young_modulus, poisson_ratio, h, Ul, flag_OPT);

  Kg = Te' * Kl * Te;
  fg = Kg * Ug;

  ks = {switchToNodalIndexing(Kg)};
  fs = {switchToNodalIndexing(fg)};

end

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
function f = myLoadSpringMass(t, UsCell)

  % Force data
  k        = 39.47; % spring constant
  m        = 1; % mass of the system
  p0       = 40; % amplitude of applied load

  % md The free vibration motion parameters are:
  omegaN       = sqrt(k / m);
  % md The frequency of the sinusoidal external force is set:
  omegaBar     = 4 * omegaN;

  f = zeros(12, 1);
  f(7) = p0 * sin(omegaBar * t);

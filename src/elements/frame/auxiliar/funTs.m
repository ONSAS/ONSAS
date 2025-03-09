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
function [Dg] = funTs(t)

  nt = norm(t);
  I = eye(3, 3);

  if nt == 0
    Dg = I;
  else
    a = 2 * (sin(nt / 2) / nt)^2;
    b = (1 - sin(nt) / nt) / nt^2;
    M = skew(t);
    Dg = I + a * M + b * M * M;
  end

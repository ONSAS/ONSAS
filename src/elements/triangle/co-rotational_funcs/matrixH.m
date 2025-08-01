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

function [H] = matrixH(t)

  spin_t = skew(t);
  I = eye(3);
  v = norm(t);

  eta = (2 * sin(v) - v * (1 + cos(v))) / (2 * v^2 * sin(v));
  H = I - 0.5 * spin_t + eta * spin_t * spin_t;

end

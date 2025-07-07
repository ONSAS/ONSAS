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
function [Ts] = Ts(t)

  psi = norm(t);
  I = eye(3, 3);

  if psi == 0
    Ts = I;
  else
    u = t / psi;
    a = sin(psi) / psi;
    Ts = a * I + (1-a) * u * u' + 1/2 * ( sin(psi/2) / (psi/2) )^2 * skew(t);
  end

end

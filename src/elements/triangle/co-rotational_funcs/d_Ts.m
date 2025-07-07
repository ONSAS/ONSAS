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
function [d_Ts] = d_Ts(t, v)

  psi = norm(t);
  I = eye(3, 3);

  if psi == 0
    d_Ts = 0.5 * skew(v);
  else
    u = t / psi;
    a = sin(psi) / psi;
    b = (sin(psi/2) / (psi/2));
    sum1 = -(a - b^2) * cross(u,v) * u';
    sum2 = 0.5 * b^2 * skew(v);
    sum3 = (cos(psi)-a) * 1 / psi * (v * u'- (u' * v) * u * u');
    sum4 = (1-a) * 1 / psi * (u * v' -2 * (u' * v) * u * u' + (u' * v) * I); 
    d_Ts = sum1 + sum2 + sum3 + sum4;
  end

end

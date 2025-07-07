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

function [Kki] = matrixKki(q, m) % ok
  % Eq. (41) of 10.1016/j.cma.2006.10.006

  q0s = 1.0 - q(1)^2 - q(2)^2 - q(3)^2;
  if q0s < 0
    q0s;
    error('q0s <0 - ojo');
  end
  q0 = sqrt(q0s);
  A = q(1) * m(1) + q(2) * m(2) + q(3) * m(3);

  H = zeros(3, 3);
  H(1, 1) = (q0s + q(1)^2) * A;
  H(2, 2) = (q0s + q(2)^2) * A;
  H(3, 3) = (q0s + q(3)^2) * A;
  H(1, 2) = q0s * (q(1) * m(2) - q(2) * m(1) - q0 * m(3)) + q(1) * q(2) * A;
  H(1, 3) = q0s * (q(1) * m(3) - q(3) * m(1) + q0 * m(2)) + q(1) * q(3) * A;
  H(2, 1) = q0s * (q(2) * m(1) - q(1) * m(2) + q0 * m(3)) + q(2) * q(1) * A;
  H(2, 3) = q0s * (q(2) * m(3) - q(3) * m(2) - q0 * m(1)) + q(2) * q(3) * A;
  H(3, 1) = q0s * (q(3) * m(1) - q(1) * m(3) - q0 * m(2)) + q(3) * q(1) * A;
  H(3, 2) = q0s * (q(3) * m(2) - q(2) * m(3) + q0 * m(1)) + q(3) * q(2) * A;

  Kki = (2 / q0^3) * H;

end

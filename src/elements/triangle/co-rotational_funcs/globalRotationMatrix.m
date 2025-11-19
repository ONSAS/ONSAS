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

function [R] = globalRotationMatrix(q) % ok
  % Eq. (35) of 10.1016/j.cma.2006.10.006
  q1 = q(1);
  q2 = q(2);
  q3 = q(3);
  q0 = sqrt(1.0 - q1^2 - q2^2 - q3^2);

  R  = 2 * [[(q0^2 + q1^2 - 0.5),  (q1 * q2 - q0 * q3),    (q1 * q3 + q0 * q2)]
            [(q1 * q2 + q0 * q3),  (q0^2 + q2^2 - 0.5),    (q2 * q3 - q0 * q1)]
            [(q1 * q3 - q0 * q2),  (q2 * q3 + q0 * q1),    (q0^2 + q3^2  - 0.5)]];
end

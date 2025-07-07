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

function [v] = rotationVector(R, flag_first_mod) % ok
  % Eq. (13) of 10.1016/j.cma.2006.10.006

  if flag_first_mod == 0
    v = zeros(3, 1);
    v(1) = .5 * (R(3, 2) - R(2, 3));
    v(2) = .5 * (R(1, 3) - R(3, 1));
    v(3) = .5 * (R(2, 1) - R(1, 2));
    if norm(v) ~= 0
      % asin(norm(v)) / norm(v)
      v = asin(norm(v)) / norm(v) * v;
    end
  else
    v = logar(R);
  end
end

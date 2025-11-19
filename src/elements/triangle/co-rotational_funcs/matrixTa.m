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

function [Ta] = matrixTa(R)
  % Eq. (15) of 10.1016/j.cma.2006.10.006
  Ta = 0.5 * [[R(2, 2) + R(3, 3),  -R(1, 2),            -R(1, 3)]
              [-R(2, 1),            R(1, 1) + R(3, 3),  -R(2, 3)]
              [-R(3, 1),           -R(3, 2),             R(1, 1) + R(2, 2)]
             ];
end

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

% This function applies a permutation to switch a vector/matrix of displacements/rotations from type to nodal
function O = switchToNodalIndexing(B)
  O = zeros(size(B));

  assert(mod(size(B, 1), 6) == 0);
  permutIndxs = [];
  for i = 1:(size(B, 1) / 6)
    permutIndxs = [permutIndxs ([1:2:5] + 6 * (i - 1)) ([2:2:6] + 6 * (i - 1))];
  end

  if size(B, 2) > 1
    O(permutIndxs, permutIndxs) = B;
  else
    O(permutIndxs) = B;
  end
end

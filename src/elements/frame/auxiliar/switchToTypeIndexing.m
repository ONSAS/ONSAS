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

% This function applies a permutation to switch a vector/matrix of displacements/rotations from nodal to type
function B = switchToTypeIndexing(O)

  assert( mod(size(O,1),6)==0)
  permutIndxs = [ ];
  for i=1:(size(O,1)/6)
    permutIndxs = [ permutIndxs ([1:2:5]+6*(i-1)) ([2:2:6]+6*(i-1)) ] ;
  end

  if size(O,2)>1
	B = O(permutIndxs, permutIndxs)
  else
	B = O(permutIndxs);
  end
end
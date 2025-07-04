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

function matFint = getInternalForces(forcesStruct, elems, fieldNames)

  nElems  = length(elems);
  nFields = length(fieldNames);

  matFint = zeros(nElems, nFields);

  for i = 1:nElems
    for j = 1:nFields
      matFint(i, j) = getfield(forcesStruct(elems(i)), fieldNames{j});
    end
  end

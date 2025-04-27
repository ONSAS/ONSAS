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
% This function provides information regarding the degrees of freedom
% of each element.
% Inputs:
%    - the elemType string
% Output:
%    - the number of nodes of the element
%    - the vector of entries of the nodal degrees of freedom
%      Assuming that the dofs per node are sorted, for node 1, as:
%      [ u_x^1 \theta_x^1 u_y^1 \theta_y^1 u_z^1 \theta_z^1 ]
function [numNodes, nodalDofsEntries] = elementTypeDofs(elemType)

  if strcmp(elemType, 'node')
    numNodes = 1;
    nodalDofsEntries = (1:2:6)';

  elseif strcmp(elemType, 'truss') || strcmp(elemType, 'edge')
    numNodes = 2;
    nodalDofsEntries = (1:2:6)';

  elseif strcmp(elemType, 'frame')
    numNodes = 2;
    nodalDofsEntries = (1:6)';

  elseif strcmp(elemType, 'tetrahedron')
    numNodes = 4;
    nodalDofsEntries = (1:2:6)';

  elseif strcmp(elemType, 'triangle')
    numNodes = 3;
    nodalDofsEntries = [1 3 5]'; % only x-y displacements

  elseif strcmp(elemType, 'triangle-plate')
    numNodes = 3;
    nodalDofsEntries = [5 2 4]'; % assumed plate surface in x-y

  elseif strcmp(elemType, 'triangle-shell')
    numNodes = 3;
    nodalDofsEntries = [1 2 3 4 5 6]'; % assumed plate surface in x-y

  end

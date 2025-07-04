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
function [vtkNodes, vtkConec, vtkNodalDisps, vtkInternalForces] ...
   = shellVtkData(Nodes, Conec, elemCrossSecParams, U, internalForces)

  vtkNodes        = [];
  vtkConec        = [];
  vtkNodalDisps   = [];

  vtkInternalForcesNames = fieldnames(internalForces);

  nelem            = size(Conec, 1);

  % thickness
  tz = elemCrossSecParams{2};

  Nodes = Nodes + reshape(U(1:2:end), 3, size(Nodes, 1))';

  indMx = 0;
  indMy = 0;
  indMxy = 0;
  vtkInternalForces = cell(length(vtkInternalForcesNames), 1);
  for i = 1:length(vtkInternalForcesNames)
    if strcmp(vtkInternalForcesNames{i}, 'Mx')
      indMx  = i;
    elseif strcmp(vtkInternalForcesNames{i}, 'My')
      indMy  = i;
    elseif strcmp(vtkInternalForcesNames{i}, 'Mxy')
      indMxy = i;
    else
      vtkInternalForces{i} = zeros(nelem, 1);
    end
  end

  for i = 1:nelem

    % nodes and degrees of freedom of current element
    nodesElem  = Conec(i, 1:3);
    dofsElem   = nodes2dofs(nodesElem, 6);

    crossVector = cross( ...
                        Nodes(nodesElem(2), :) - Nodes(nodesElem(1), :), ...
                        Nodes(nodesElem(3), :) - Nodes(nodesElem(1), :) ...
                       );
    normalVector = crossVector / norm(crossVector);

    coordsThickness = tz * .5 * [ones(3, 1); -ones(3, 1)] * normalVector;

    Nodesvtk = [Nodes(nodesElem, :);   ...
                Nodes(nodesElem, :)] + ...
               coordsThickness;

    % column vector with displacements of the dofs of the current element
    dispsElem  = U(dofsElem);

    aux  = reshape(dispsElem(1:2:end)', [3, 3])';
    rots = reshape(dispsElem(2:2:end)', [3, 3])';

    Dispsvtk = ...
     [(expon(rots(1, :)) * coordsThickness(1, :)')'; ...
      (expon(rots(2, :)) * coordsThickness(2, :)')'; ...
      (expon(rots(3, :)) * coordsThickness(3, :)')'; ...
      (expon(rots(1, :)) * coordsThickness(4, :)')'; ...
      (expon(rots(2, :)) * coordsThickness(5, :)')'; ...
      (expon(rots(3, :)) * coordsThickness(6, :)')'] - ...
     coordsThickness + ...
     [aux; aux];

    Conecvtk = [13  (nodes2dofs(i, 6) - 1)'];

    vtkNodes             = [vtkNodes;     Nodesvtk];
    vtkConec             = [vtkConec;     Conecvtk];
    vtkNodalDisps        = [vtkNodalDisps; Dispsvtk];

    if indMx > 0
      vtkInternalForces{indMx}  = [vtkInternalForces{indMx};  internalForces(i).Mx * ones(size(Conecvtk, 1), 1)];
    end
    if indMy > 0
      vtkInternalForces{indMy}  = [vtkInternalForces{indMy};  internalForces(i).My * ones(size(Conecvtk, 1), 1)];
    end
    if indMxy > 0
      vtkInternalForces{indMxy} = [vtkInternalForces{indMxy}; internalForces(i).Mxy * ones(size(Conecvtk, 1), 1)];
    end

  end % for elements

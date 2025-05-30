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
  = trussVtkData(Nodes, Conec, elemCrossSecParams, U, internalForces)

  vtkNodes        = [];
  vtkConec        = [];

  vtkNodalDisps   = [];
  vtkInternalForcesNames   = fieldnames(internalForces);

  nelem = size(Conec, 1);

  indNx = 0;
  vtkInternalForces = cell(length(vtkInternalForcesNames), 1);

  for i = 1:length(vtkInternalForcesNames)
    if strcmp(vtkInternalForcesNames{i}, 'Nx')
      indNx = i;
    else
      vtkInternalForces{i} = zeros(nelem, 1);
    end
  end

  counterNodes     = 0;

  for i = 1:nelem
    nodesElem    = Conec(i, 1:2);
    dofsElem     = nodes2dofs(nodesElem, 6);

    coordsElemNodes = reshape(Nodes(nodesElem(:), :)', 6, 1);

    % length of current element
    elemLength = norm(coordsElemNodes(4:6) - coordsElemNodes(1:3));

    % q section tengo
    [iniNodes, midNodes, endNodes, sectPar] = crossSectionVtkSolidConnec(elemCrossSecParams);

    dispsElem          = U(dofsElem);
    thetaLocIniSubElem = zeros(3, 1);
    thetaLocEndSubElem = zeros(3, 1);

    [R0, Rr, locDisp] = elementBeamRotData(coordsElemNodes, dispsElem);

    dispLocIniSubElem  = zeros(3, 1);
    dispLocEndSubElem  = [locDisp(1); 0; 0];

    coordLocSubElem = [0; elemLength];

    [Nodesvtk, Conecvtk, Dispsvtk] = vtkBeam2SolidConverter(coordsElemNodes, ...
                                                            dispsElem, coordLocSubElem, dispLocIniSubElem, dispLocEndSubElem, thetaLocIniSubElem, thetaLocEndSubElem, sectPar, Rr, R0);

    Conecvtk(:, 2:end) = Conecvtk(:, 2:end) + counterNodes;
    vtkNodes             = [vtkNodes;     Nodesvtk];
    vtkConec             = [vtkConec;     Conecvtk];
    vtkNodalDisps        = [vtkNodalDisps; Dispsvtk];
    if indNx > 0
      vtkInternalForces{indNx} = [vtkInternalForces{indNx}; internalForces(i).Nx * ones(size(Conecvtk, 1), 1)];
    end
    counterNodes = counterNodes + size(Nodesvtk, 1);
  end % for plot element subdivision

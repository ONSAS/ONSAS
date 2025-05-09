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
   = frameVtkData(Nodes, Conec, elemCrossSecParams, U, internalForces)

  vtkNodes        = [];
  vtkConec        = [];

  vtkNodalDisps   = [];
  vtkInternalForcesNames = fieldnames(internalForces);

  nelem = size(Conec, 1);

  indNx = 0;
  indMy = 0;
  indMz = 0;
  vtkInternalForces = cell(length(vtkInternalForcesNames), 1);
  for i = 1:length(vtkInternalForcesNames)
    if strcmp(vtkInternalForcesNames{i}, 'Nx')
      indNx = i;
    elseif strcmp(vtkInternalForcesNames{i}, 'My')
      indMy = i;
    elseif strcmp(vtkInternalForcesNames{i}, 'Mz')
      indMz = i;
    else
      vtkInternalForces{i} = zeros(nelem, 1);
    end
  end

  nPlotSubElements = 10; % number of plot subsegments
  counterNodes     = 0;

  for i = 1:nelem

    % nodes and degrees of freedom of current element
    nodesElem  = Conec(i, 1:2);
    dofsElem   = nodes2dofs(nodesElem, 6);

    % column vector with the coordinates of the nodes of the element
    coordsElemNodes = reshape(Nodes(nodesElem(:), :)', 6, 1);

    % computes connectivity of vtk element associated with frame cross-section
    [iniNodes, midNodes, endNodes, sectPar] = crossSectionVtkSolidConnec(elemCrossSecParams);

    % length of current element
    elemLength = norm(coordsElemNodes(4:6) - coordsElemNodes(1:3));

    % column vector with displacements of the dofs of the current element
    dispsElem  = U(dofsElem);

    % column vector with discretization in subelements, using local coordinates (r1,r2,r3)
    xsloc  = linspace(0, elemLength, nPlotSubElements + 1)';

    % interpolation functions evaluation matrix with columns
    %    [  N1(lin1)                          N2(lin2)            N3(v1) N4(theta1) N5(v2) N6(theta2) ]
    interFuncLinear = [(elemLength - xsloc) / elemLength  xsloc / elemLength];
    interFuncCubic  = bendingInterFuns(xsloc, elemLength, 0);
    interFuncQuad   = bendingInterFuns(xsloc, elemLength, 1);

    [R0, Rr, locDisp] = elementBeamRotData(coordsElemNodes, dispsElem);

    ul              = locDisp(1);
    thetaLocIniElem = locDisp(2:4);
    thetaLocEndElem = locDisp(5:7);

    % interpolation of displacements
    valsLocDispXSubElements = interFuncLinear * [0;  ul];
    valsLocDispYSubElements = interFuncCubic  * [0;  thetaLocIniElem(3); 0;  thetaLocEndElem(3)];
    valsLocDispZSubElements = interFuncCubic   * [0; -thetaLocIniElem(2); 0; -thetaLocEndElem(2)];

    % interpolation of angles
    valsLocThetaXSubElements = interFuncLinear * [thetaLocIniElem(1);    thetaLocEndElem(1)];
    valsLocThetaYSubElements = interFuncQuad   * [0; thetaLocIniElem(2); 0; thetaLocEndElem(2)];
    valsLocThetaZSubElements = interFuncQuad   * [0; thetaLocIniElem(3); 0; thetaLocEndElem(3)];

    for j = 1:nPlotSubElements

      dispLocIniSubElem = [valsLocDispXSubElements(j); ...
                           valsLocDispYSubElements(j); ...
                           valsLocDispZSubElements(j)];

      dispLocEndSubElem = [valsLocDispXSubElements(j + 1); ...
                           valsLocDispYSubElements(j + 1); ...
                           valsLocDispZSubElements(j + 1)];

      % 2-compute local rotations for the initial and final sections
      thetaLocIniSubElem = [valsLocThetaXSubElements(j); ...
                            valsLocThetaYSubElements(j); ...
                            valsLocThetaZSubElements(j)];
      %
      thetaLocEndSubElem = [valsLocThetaXSubElements(j + 1); ...
                            valsLocThetaYSubElements(j + 1); ...
                            valsLocThetaZSubElements(j + 1)];

      coordLocSubElem = [xsloc(j); xsloc(j + 1)];

      [Nodesvtk, Conecvtk, Dispsvtk] = vtkBeam2SolidConverter(coordsElemNodes, ...
                                                              dispsElem, coordLocSubElem, dispLocIniSubElem, dispLocEndSubElem, thetaLocIniSubElem, thetaLocEndSubElem, sectPar, Rr, R0);

      Conecvtk(:, 2:end) = Conecvtk(:, 2:end) + counterNodes;
      vtkNodes             = [vtkNodes;     Nodesvtk];
      vtkConec             = [vtkConec;     Conecvtk];
      vtkNodalDisps        = [vtkNodalDisps; Dispsvtk];

      if indNx > 0
        vtkInternalForces{indNx} = [vtkInternalForces{indNx}; internalForces(i).Nx * ones(size(Conecvtk, 1), 1)];
      end
      if indMy > 0
        vtkInternalForces{indMy} = [vtkInternalForces{indMy}; internalForces(i).My * ones(size(Conecvtk, 1), 1)];
      end
      if indMz > 0
        vtkInternalForces{indMz} = [vtkInternalForces{indMz}; internalForces(i).Mz * ones(size(Conecvtk, 1), 1)];
      end

      counterNodes = counterNodes + (size(Conecvtk, 2) - 1);

    end % for plot points

  end % for elements

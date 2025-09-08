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
% md function that constructs the assembled Fext vector for one given BC

function fext = elem2NodalLoads (Conec, indBC, elements, boundaryCond, Nodes)

  % declare output fext
  nnodes      = size(Nodes, 1);
  fext        = zeros(6 * nnodes, 1);

  % get element indexes with current BC
  elemsWithBC = find(Conec(:, 3) == indBC);

  % extract BC load base vals and coord system
  loadCoordSys = boundaryCond.loadsCoordSys;
  loadvals     = boundaryCond.loadsBaseVals;

  loadedNodes = [];

  % loop in elements with current BC
  for elemInd = 1:length(elemsWithBC)

    elem      = elemsWithBC(elemInd);
    elemInd   = Conec(elem, 2);
    elemType  = elements(elemInd).elemType;

    % md nodal loads
    if strcmp(elemType, 'node') % node

      if strcmp(loadCoordSys, 'global')
        nodes     = Conec(elem, 3 + 1);
      else
        error(' only global flag in load by now.');
      end
      elemNodeLoadsMatrix = loadvals;

      % md truss
    elseif strcmp(elemType, 'truss')   %
      error(' not yet.');

      % md frame
    elseif strcmp(elemType, 'frame')   %
      error(' not yet.');

      % md edge
    elseif strcmp(elemType, 'edge')   %
      nodes          = Conec(elem, 3 + (1:2));

      % vector from node 1 to node 2
      directionVector = Nodes(nodes(2), :) - Nodes(nodes(1), :);
      % length of the edge
      lengthElem = norm(directionVector);
      % thickness of the edge element
      thickness  = elements(elemInd).elemCrossSecParams;

      % check oriented vector is defined into x-y plane
      % assert(abs(directionVector(3)) < eps * 10, 'edge must be defined into x-y plane');
      % if abs(Nodes(nodes(1), 3)) > eps * 10 || abs(Nodes(nodes(2), 3)) > eps * 10 
        % error('edge must be defined into x-y plane');
      % end
      
      % % factor of the load for each node
      factor = lengthElem * thickness * 0.5;

      if strcmp(loadCoordSys, 'global')
        qx = loadvals(1) ;
        qy = loadvals(2) ;
        xaxis = [1,0,0];
        
        normalVector = cross(xaxis,directionVector);
        nvector = normalVector / norm(normalVector);
        if norm(normalVector) == 0 
          angle = 0;
        else
          angle = asin(normalVector/( norm(xaxis)*norm(directionVector)*nvector )); % cross(a,b) = |a|*|b|*sin(theta)*n
        end
        q_perp = qx*sin(angle) - qy*cos(angle);
                
        % Aedge = lengthElem * thickness; % edge area
        Mz_q = q_perp*thickness*lengthElem^2/12; % 
        % Mz_q = qx*thickness*lengthElem^2/8; % 
        % Mz_q=0;
        Fx = loadvals(1) * factor;
        Mx = loadvals(2) * factor;
        Fy = loadvals(3) * factor;
        My = loadvals(4) * factor;
        Fz = loadvals(5) * factor;
        Mz = loadvals(6) * factor;
        
        elemNodeLoadsMatrix = [Fx Mx Fy My Fz Mz+Mz_q; Fx Mx Fy My Fz Mz-Mz_q] ;
        % elemNodeLoadsMatrix
      elseif strcmp(loadCoordSys, 'local')
        % consider a 90 degrees rotation of the oriented vector of the line element
        % tangent unitary vector
        tangUniVec = directionVector / lengthElem;
        % normal unitary vector
        normalUniVec = cross([0 0 1], tangUniVec);
        % tension vector
        assert(norm(loadvals([2 4 5 6])) < eps);
        tensionVec = (loadvals(1) * tangUniVec  + loadvals(3) * normalUniVec) * factor;
        % nodal forces in global coordinates
        Fx = tensionVec(1);
        Fy = tensionVec(2);
        Fz = 0;
        Mx = 0;
        My = 0;
        Mz = 0;
        elemNodeLoadsMatrix = ones(length(nodes), 1) * [Fx Mx Fy My Fz Mz];
      
      end % if global/local system

      assert(size(elemNodeLoadsMatrix, 2) == 6, 'error, maybe missing thickness');

      % md triangle tension
    elseif strcmp(elemType, 'triangle') || strcmp(elemType, 'triangle-plate') || strcmp(elemType, 'triangle-shell')   %

      nodes = Conec(elem, 3 + (1:3));

      % area of the element with direction
      crossVector = cross( ...
                          Nodes(nodes(2), :) - Nodes(nodes(1), :), ...
                          Nodes(nodes(3), :) - Nodes(nodes(1), :) ...
                         );

      areaElem = norm(crossVector) / 2;

      if strcmp(loadCoordSys, 'global')

        Fx = loadvals(1) * areaElem / 3;
        Fy = loadvals(3) * areaElem / 3;
        Fz = loadvals(5) * areaElem / 3;

        assert(sum(loadvals([2 4 6]) == 0) == 3, ...
               'error only pressure loads, not moments, create an issue!');

      elseif strcmp(loadCoordSys, 'local') % local coordinates load

        normalVector = crossVector / norm(crossVector);

        dofsaux = nodes2dofs(nodes, 6);
        dofs    = dofsaux(1:2:length(dofsaux));

        %  x, y and z components of the normal tension
        Fx = normalVector(1) * loadvals(5) / 3  * areaElem;
        Fy = normalVector(2) * loadvals(5) / 3  * areaElem;
        Fz = normalVector(3) * loadvals(5) / 3  * areaElem;

        assert(sum(loadvals([1 2 3 4 6]) == 0) == 5, ...
               'error only normal pressure loads in local coords, create an issue!');

      else
        loadCoordSys;
        error(' loadsCoordSys field must be local or global.');
      end % if global/local system

      elemNodeLoadsMatrix = ones(length(nodes), 1) * [Fx 0 Fy 0 Fz 0];

    end % if elemTypes
    % mdadd loads to matrix of loaded nodes
    loadedNodes = [loadedNodes; ...
                   nodes'  elemNodeLoadsMatrix];
  end % for elements

  % md convert to assembled fext vector
  if exist('loadedNodes') ~= 0
    for i = 1:size(loadedNodes, 1)
      aux = nodes2dofs (loadedNodes(i, 1), 6);
      fext(aux) = fext(aux) + loadedNodes(i, 2:7)';
    end
  end


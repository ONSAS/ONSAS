% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro
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


function [ U, Udot, Udotdot ] = initialCondsProcessing( mesh, initialConds )

% Build nodes MEBI matrix by selecting the first Conec index null 
conecCellMat = myCell2Mat( mesh.conecCell )                ;
conecNodes = conecCellMat(find( conecCellMat(:,1)==0 ),:)  ;

% Build nodes MEBI matrix
conecElems = conecCellMat(find( conecCellMat(:,1) ~=0 ),:) ;

% Compute number of nodes
nNodes = size( mesh.nodesCoords, 1 ) ;

% Create velocity and displacements vectors
U       = zeros( 6*nNodes,   1 ) ;
Udot    = zeros( 6*nNodes,   1 ) ;
Udotdot = zeros( 6*nNodes,   1 ) ;

% Nodes initial conditions

% Process nodal displacement initial condition (IC)
if isfield( initialConds, 'nonHomogeneousInitialCondU0' )

  % Compute different initial conditions loaded 
  initialCondsTypes  = unique( conecNodes( :, 4) ) ;

  %md loop over the types of initial conditions added in the mesh
  for indIC = 1:length( initialCondsTypes )

    % number of current IC processed
    ICnum = initialCondsTypes( indIC ) ;

    % locate nodes with that index
    indexNodesIC = find( conecNodes( :, 4 ) == indIC ) ;

    % load initial condition value and dof
    ICindex = initialConds.nonHomogeneousInitialCondU0( indIC, : ) ;

    % compute dofs index to add the IC into the displacement vector
    dofsNodesIC = 6 * (indexNodesIC - 1) + ICindex(1) ;

    % add the initial condition displacements vector
    U(dofsNodesIC, 1 ) = ICindex(2)
  end
end
% Process nodal velocity initial condition 
if isfield( initialConds, 'nonHomogeneousInitialCondUdot0' )

  % Compute different initial conditions loaded 
  initialCondsTypes  = unique( conecNodes( :, 4) ) ;

  %md loop over the types of initial conditions added in the mesh
  for indIC = 1:length( initialCondsTypes )

    % number of current IC processed
    ICnum = initialCondsTypes( indIC ) ;

    % locate nodes with that index
    indexNodesIC = find( conecNodes( :, 4 ) == indIC ) ;

    % load initial condition value and dof
    ICindex = initialConds.nonHomogeneousInitialCondUdot0( indIC, : ) ;

    % compute dofs index to add the IC into the displacement vector
    dofsNodesIC = 6 * (indexNodesIC - 1) + ICindex(1) ;

    % add the initial condition displacements vector
    Udot(dofsNodesIC, 1 ) = ICindex(2)
  end
end
% Process nodal aceleration initial condition 
if isfield( initialConds, 'nonHomogeneousInitialCondUdotdot0' )

  % Compute different initial conditions loaded 
  initialCondsTypes  = unique( conecNodes( :, 4) ) ;

  %md loop over the types of initial conditions added in the mesh
  for indIC = 1:length( initialCondsTypes )

    % number of current IC processed
    ICnum = initialCondsTypes( indIC ) ;

    % locate nodes with that index
    indexNodesIC = find( conecNodes( :, 4 ) == indIC ) ;

    % load initial condition value and dof
    ICindex = initialConds.nonHomogeneousInitialCondUdotdot0( indIC, : ) ;

    % compute dofs index to add the IC into the displacement vector
    dofsNodesIC = 6 * (indexNodesIC - 1) + ICindex(1) ;

    % add the initial condition displacements vector
    Udotdot(dofsNodesIC, 1 ) = ICindex(2)
  end
end

% Add elements initial condition here
% ---------------------------------------------------

% computation of initial acceleration for some cases
% ---------------------------------------------------

% --- initial tangent matrices ---
%~ [ mats ] = assembler ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut, dynamicAnalysisBoolean, 2, Udotdott, massMatType ) ;

%~ systemDeltauMatrix = mats{1}

%~ stop
%~ if dynamicAnalysisBoolean == 1,

  %~ massMat    = mats{2} ;

  %~ % --- computation of initial Udotdott for truss elements only!!!
  %~ Fext = computeFext( constantFext, variableFext, loadFactors(1), userLoadsFilename ) ;

  %~ Udotdott (neumdofs) = massMat( neumdofs, neumdofs ) \ ( Fext(neumdofs) -Fintt( neumdofs ) ) ;

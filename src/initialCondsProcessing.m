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


% This script declares several matrices and vectors required for the analysis. In this script, the value of important magnitudes, such as internal and external forces, displacements, and velocities are computed for step/time 0.


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
  % delete null initial conditions
  initialCondsTypes(initialCondsTypes==0) = [] ;
  %md loop over the types of initial conditions added in the mesh
  for indIC = 1:length( initialCondsTypes ) ;

    % number of current IC processed
    ICnum = initialCondsTypes( indIC ) ;

    % locate nodes with that index of initial condition into conecNodes matrix
    indexNodesIC = find( conecNodes( :, 4 ) == ICnum ) ;

    % get the nodes that has that IC
    nodesIC = conecNodes(indexNodesIC,5) ;
    
    % add the imposedDisps for all the nodes with the same IC
    % this could be vectorized! using the Dofs relative to the dof of the node dofNodes = nodes2dofs(nodesIC, 6) and the dofNodes(1:dofsIC:end) and so on
    for nodeIC = nodesIC'
      % compute dofs of nodes with that IC 
      dofsNodesIC = nodes2dofs(nodeIC, 6) ;

      % load initial condition values
      dofsICindex = initialConds(indIC).nonHomogeneousInitialCondU0.Dofs ; 
      valsICindex = initialConds(indIC).nonHomogeneousInitialCondU0.Vals ; 

      % add the initial condition displacements vector
      U( dofsNodesIC(dofsICindex), 1 ) = valsICindex ;
    end
  end
end
% Process nodal velocity initial condition 
if isfield( initialConds, 'nonHomogeneousInitialCondUdot0' )

  % Compute different initial conditions loaded 
  initialCondsTypes  = unique( conecNodes( :, 4) ) ; 
  % delete null initial conditions
  initialCondsTypes(initialCondsTypes==0) = [] ;
  %md loop over the types of initial conditions added in the mesh
  for indIC = 1:length( initialCondsTypes ) ;

    % number of current IC processed
    ICnum = initialCondsTypes( indIC ) ;

    % locate nodes with that index of initial condition into conecNodes matrix
    indexNodesIC = find( conecNodes( :, 4 ) == ICnum ) ;

    % get the nodes that has that IC
    nodesIC = conecNodes(indexNodesIC,5) ;
    
    % add the imposedDisps for all the nodes with the same IC
    % this could be vectorized! using the Dofs relative to the dof of the node dofNodes = nodes2dofs(nodesIC, 6) and the dofNodes(1:dofsIC:end) and so on
    for nodeIC = nodesIC'
      % compute dofs of nodes with that IC 
      dofsNodesIC = nodes2dofs(nodeIC, 6) ;

      % load initial condition values
      dofsICindex = initialConds(indIC).nonHomogeneousInitialCondUdot0.Dofs ; 
      valsICindex = initialConds(indIC).nonHomogeneousInitialCondUdot0.Vals ; 

      % add the initial condition displacements vector
      Udot( dofsNodesIC(dofsICindex), 1 ) = valsICindex ;
    end
  end
end

% Add initial acceleration for some numerical method ?

% Add elements initial condition here
% ---------------------------------------------------
% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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


function [ U, Udot, Udotdot ] = initialCondsProcessing( mesh, initialConds, elements )

  % Extract mesh parameters
  Conec  = myCell2Mat( mesh.conecCell ) ;
  Nodes  = mesh.nodesCoords ;
  nNodes = size( Nodes,1) ;

  % Create kinematic vectors
  U       = zeros( 6*nNodes, 1 ) ;
  Udot    = zeros( 6*nNodes, 1 ) ;
  Udotdot = zeros( 6*nNodes, 1 ) ;

  % displacements initial conditions
  if isfield( initialConds, 'nonHomogeneousUVals' )
    
    % Compute different initial conditions loaded 
    initialCondsTypes  = unique( Conec( :, 4) ) ;
 
    % delete null initial conditions
    if initialCondsTypes(1) == 0
      initialCondsTypes(1) = [] ;
    end
    
    %md loop over the types of initial conditions added in the mesh
    for indIC = 1:length( initialCondsTypes ) ;

      % number of current IC processed
      indIC = initialCondsTypes( indIC ) ;
     
      %md find the elements with the current initial condition
      elemsWithIC = find( Conec(:,4) == indIC ) ;
     
      %md values and imposed dofs of current IC
      impoUDofs = initialConds(indIC).nonHomogeneousUDofs ;
      impoUVals = initialConds(indIC).nonHomogeneousUVals ;
     
      % compute the imposed dofs and vals for the elements with that IC 
      [ nonHomUDiriVals, diriDofs, nonHomoDiriDofs ] = elem2NodalDisps ( Conec, indIC, elemsWithIC, elements, impoUDofs, impoUVals, Nodes ) ;
      
      % add the initial condition velocity vector
      U( diriDofs, 1 ) = U( diriDofs, 1 ) + nonHomUDiriVals ;
    
    end %for types of IC
  
  end %if disp IC

  % Process velocity initial condition 
  % displacements initial conditions
  if isfield( initialConds, 'nonHomogeneousUdotVals' )
    
    % Compute different initial conditions loaded 
    initialCondsTypes  = unique( Conec( :, 4) ) ;
 
    % delete null initial conditions
    if initialCondsTypes(1) == 0
      initialCondsTypes(1) = [] ;
    end
    
    %md loop over the types of initial conditions added in the mesh
    for indIC = 1:length( initialCondsTypes ) ;

      % number of current IC processed
      indIC = initialCondsTypes( indIC ) ;
     
      %md find the elements with the current initial condition
      elemsWithIC = find( Conec(:,4) == indIC ) ;
     
      %md values and imposed dofs of current IC
      impoUDofs = initialConds(indIC).nonHomogeneousUdotDofs ;
      impoUVals = initialConds(indIC).nonHomogeneousUdotVals ;
     
      % compute the imposed dofs and vals for the elements with that IC 
      [ nonHomUDiriVals, ~, nonHomoDiriDofs ] = elem2NodalDisps ( Conec, indIC, elemsWithIC, elements, impoUDofs, impoUVals, Nodes ) ; 
      
      % add the initial condition velocity vector
      Udot( nonHomoDiriDofs, 1 ) = Udot( nonHomoDiriDofs, 1 ) + nonHomUDiriVals ;
    
    end %for types of IC
  
  end %if disp IC

end %end function

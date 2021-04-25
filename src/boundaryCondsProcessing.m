% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
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

function [ Conec, Nodes, factorLoadsFextCell, loadFactorsFuncCell, diriDofs, neumDofs, KS, userLoadsFilename ] = boundaryCondsProcessing ( mesh, ...
                        materials, ...     % M
                        elements, ...      % E
                        boundaryConds, ... % B
                        initialConds ) 		 % I

Conec = myCell2Mat( mesh.conecCell ) ;
Nodes = mesh.nodesCoords ;

nnodes = size( Nodes,1);

% compute num elements and BC we have 
elementTypes   = unique( Conec( :, 2) ) ;
boundaryTypes  = unique( Conec( :, 3) ) ;

if boundaryTypes(1) == 0, % removes elements without BCs
  boundaryTypes(1)=[];
end

if elementTypes(1)  == 0, % checks if all elements have a type
  error('all elements must be defined');
end

factorLoadsFextCell = {};
loadFactorsFuncCell = {};
diriDofs            = [];

for indBC = 1:length( boundaryTypes )
  
  % loads verification
  % ------------------
  if ~isempty( boundaryConds.loadsCoordSys{ boundaryTypes(indBC) } ) % if load applied
    
    factorLoadsFextCell{ boundaryTypes(indBC) }  = elem2NodalLoads ( Conec, boundaryTypes(indBC), elements, boundaryConds, Nodes ) ;
    
    loadFactorsFuncCell{ boundaryTypes(indBC) }  = boundaryConds.loadTimeFact{ boundaryTypes(indBC) } ;
  end % if load
  
  % displacement verification
  % -------------------------
  if ~isempty( boundaryConds.impoDispDofs{ boundaryTypes(indBC) } ),
    [ nonHomDiriVals, bcDiriDofs, nonHomDiriDofs ]  = elem2NodalDisps ( Conec, boundaryTypes(indBC), elements, boundaryConds, Nodes ) ; 
    
  end % if: disp dofs
  
  diriDofs = [ diriDofs; bcDiriDofs ] ;
  
end % for: elements with boundary condition assigned

diriDofs = unique( diriDofs) ;


% construction of neumandofs
% --------------------------

neumDofs = zeros( 6*nnodes, 1 ) ; % maximum possible vector
elemsToRemove = [] ;

% loop for construction of vector of dofs
for elemNum = 1:length( elementTypes )
  elementsNums = find( Conec( :, 2 ) == elemNum ) ;
  
  elemType = elements.elemType{elemNum} ;
   
  if strcmp( elemType, 'node') 
    elemsToRemove = [ elemsToRemove ; elementsNums ] ;
  
  elseif length( elementsNums ) > 0
    
    [numNodes, dofsStep] = elementTypeInfo ( elemType ) ;
    nodes    = Conec( elementsNums, (4+1):(4+numNodes) ) ; 
    dofs     = nodes2dofs( nodes, 6)'       ;
    dofs     = dofs(1:dofsStep:end)         ;
  
    neumDofs ( dofs ) = dofs ;
  end
end

if length( diriDofs)>0
  neumDofs( diriDofs ) = 0 ;
end
neumDofs = unique( neumDofs ) ;
if neumDofs(1) == 0,
  neumDofs(1)=[] ;
end

Conec( elemsToRemove, :  ) = [] ;

if isfield( boundaryConds,'userLoadsFilename' )
  userLoadsFilename = boundaryConds.userLoadsFilename ;
else
  userLoadsFilename =[];
end
% ----------------------------------------------------------------------


% ----------------------
KS        = sparse( 6*nnodes, 6*nnodes );  

%~ for i=1:size(nodalSprings,1)
  %~ aux = nodes2dofs ( nodalSprings (i,1) , 6 ) ;
  %~ for k=1:6
    %~ %
    %~ if nodalSprings(i,k+1) == inf,
      %~ fixeddofs = [ fixeddofs; aux(k) ] ;
    %~ elseif nodalSprings(i,k+1) > 0,
      %~ KS( aux(k), aux(k) ) = KS( aux(k), aux(k) ) + nodalSprings(i,k+1) ;
    %~ end
  %~ end
%~ end




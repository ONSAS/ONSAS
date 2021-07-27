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

% Script for verification of the input variables definition. Default
% values ar assigned.

function  [ materialsParams, elementsParams, loadsParams, crossSecsParams, springsParams, loadFactorsFunc, Conec ] = structsToCellsConv( materials, elements, boundaryConds, initialConds, Conec )

% materials
nMats = length( materials.hyperElasModel )
materialsParams = {} ;
for i=1:nMats
  if strcmp( materials.hyperElasModel, 'SVK' )
    materialsParams{i,1} = [ 2 materials.hyperElasParams{i} ]
  end
end

% elements and cross
nElem = length( elements.elemType )
elementsParams = {} ;
crossSecsParams = {};
vecOldCross = [];
ncros = 0;
for i=1:nElem
  if strcmp( elements.elemType{i}, 'node' )
    elementsParams{i,1} = 1 ;
    vecOldCross(i)= 0 ;
  elseif strcmp( elements.elemType{i}, 'truss' )
    elementsParams{i,1} = 2 ;
    ncros = ncros+1 ;
    crossSecsParams{ncros,1} = elements.elemTypeGeometry{i} ;
    vecOldCross(i)= ncros;
  end
end

elementsParams


% loads and springs
nBConds =  length( boundaryConds.loadCoordSys )
loadsParams = {} ;
springsParams = {};
nloads = 0;
vecOldLoads = [];
vecOldSprings = [];
nsprings = 0 ;
for i=1:nBConds
  auxLoad = zeros(1, 8 ) ;
  auxSpri = zeros(1, 6 ) ;
  if strcmp( boundaryConds.loadCoordSys{i},'global' )
    auxLoad(1) = 1 ;
    nloads = nloads +1;
  elseif strcmp( boundaryConds.loadCoordSys{i},'local' )
    auxLoad(1) = 0 ;
    nloads = nloads+1;
  elseif ~isempty(boundaryConds.loadCoordSys{i})
    error('load coord case error')
  end

  if is_function_handle( boundaryConds.loadTimeFact{i} )
    loadFactorsFunc = boundaryConds.loadTimeFact{i}
    auxLoad(2) = 1 ;
    aux = boundaryConds.loadBaseVals{i} ;
    auxLoad( [3:8] ) = aux ;
    loadsParams{ nloads,1} = auxLoad ;
    vecOldLoads(i)= nloads ;
  end

  if ~isempty( boundaryConds.impoDispDofs{i} ) ;
    auxSpri(boundaryConds.impoDispDofs{i} ) = boundaryConds.impoDispVals{i} ;
    nsprings = nsprings+1;
    springsParams{nsprings,1} = auxSpri ;
    vecOldSprings(i)= nsprings ;
  end
end
loadsParams


for j=1:size(Conec,1)
j
oldConec = Conec{j,1};
  newconec(1) = oldConec(1);

  newconec(2) = oldConec(2);
  newconec(4) = oldConec(2);

  if oldConec(3)>0
    newconec(3) = vecOldLoads(   oldConec(3) ) ;
    newconec(5) = vecOldSprings( oldConec(3) ) ;
  end
  newconec(6:(6+length(oldConec)-5)) = oldConec(5:end)
end

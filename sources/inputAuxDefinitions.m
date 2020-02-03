% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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


tic
% ---------------------------------------------------

ndofpnode = 6; 

strucsize = strucSize(Nodes) ;

%~ absoluteimperfection = imperfactor * strucsize ;

if nonLinearAnalysisBoolean == 0 && dynamicAnalysisBoolean == 0
  controlDofInfo = [ ] ;
else 
  % control dof info
  controlDof       = nodes2dofs( controlDofInfo(1), 6)(controlDofInfo(2)) ;
  controlDofFactor = controlDofInfo(3) ; 
end


%~ releasesDofs = [];
%~ for i=1:nelems
  %~ if Conec(i,7)==1
    %~ releasesDofs = [ releasesDofs; nodes2dofs( Conec(i,1:2),ndofpnode)(2:2:end) ];
    %~ Releases = [ Releases; i ones(1,4) ] ;
  %~ end
%~ end

%~ releasesDofs = unique( releasesDofs);

coordsElemsMat = zeros(nelems,1) ;

for i=1:nelems
  % obtains nodes and dofs of element
  nodeselem = Conec(i, find(Conec(i,1:4)>0) )' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
  for j=1:length(nodeselem)
    coordsElemsMat( i, (j-1)*6+[1:2:5] ) = Nodes( nodeselem(j),:) ;
  end
end

%~ coordsElemsMat
%~ stop


hyperElasParamsMat = [] ;
for i=1:size(hyperElasParams,1)
  hyperElasParamsMat (i,1:length(hyperElasParams{i})) = hyperElasParams{i} ;
end



% ---------------- load vectors assembly -----------------------
variableFext = zeros( ndofpnode*nnodes , 1 );
constantFext = zeros( ndofpnode*nnodes , 1 );

if exist( 'nodalVariableLoads' ) ~= 0
  for i=1:size(nodalVariableLoads,1)
    aux = nodes2dofs ( nodalVariableLoads(i,1), ndofpnode ) ;
    variableFext( aux ) = variableFext( aux ) + nodalVariableLoads(i,2:7)' ;
  end
end

if exist( 'nodalConstantLoads' ) ~= 0
  for i=1:size(nodalConstantLoads,1)
    aux = nodes2dofs ( nodalConstantLoads(i,1), ndofpnode ) ;
    constantFext( aux ) = constantFext( aux ) + nodalConstantLoads(i,2:7)' ;
  end
end


% ------------------------------------------------------------
if exist( 'nodalConstantLoads' ) ~= 0 || exist( 'nodalVariableLoads' ) ~= 0
  [maxNorm2F, visualloadfactor] = visualLoadFac(strucsize, variableFext, constantFext, nnodes) ;
else
  visualloadfactor = 1 ;
end

% -------------------- indexes computation --------------------
indexesElems = zeros(nelems,1) ;
ntruss = 0 ;
nbeam = 0 ;
ntet = 0 ;
nplate = 0 ;
trussElem = [] ;
beamElem = [] ;
beamNodes = [] ;
tetElem = [] ;
plateElem = [] ;

for i = 1:nelems
  if Conec(i,7) == 1
    nodeselem = Conec(i,1:2) ;
    ntruss = ntruss + 1 ;
    indexesElems(i) = nbeam+1 ;
    trussElem = [ trussElem ; i ] ;
  elseif Conec(i,7) == 2
    nodeselem = Conec(i,1:2) ;
    nbeam = nbeam + 1 ;
    indexesElems(i) = nbeam ;
    beamElem = [ beamElem ; i ] ;
    beamNodes = [ beamNodes ; nodeselem' ] ;
  elseif Conec(i,7) == 3
    ntet = ntet + 1 ;
    indexesElems(i) = ntet ;
    tetElem = [ tetElem ; i ] ;
  elseif Conec(i,7) == 4
    nplate = nplate + 1 ;
    indexesElems(i) = nplate ;
    plateElem  = [ plateElem ; i ] ;
  end
end

% Output parameters

cellStress = [] ;
matNts = [] ;
matUts = [] ;
% ------------------------------------------------------------
tInputAuxDefs = toc;

%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.


%This script declares several matrices and vectors required for the analysis. In this script, the value of important magnitudes, such as internal forces, are computed for step/time 0. TEST
% TEST

% ----------- fixeddofs and spring matrix computation ---------
fixeddofs = [] ;
KS      = sparse(ndofpnode*nnodes,ndofpnode*nnodes);  

for i=1:size(nodalSprings,1)
  aux = nodes2dofs ( nodalSprings (i,1) , ndofpnode ) ;
  for k=1:ndofpnode
    %
    if nodalSprings(i,k+1) == inf
      fixeddofs = [ fixeddofs; aux(k) ] ;
    else
      KS(aux(k), aux(k) ) = nodalSprings(i,k+1) ;
    end
  end
end

diridofs = fixeddofs ;
diridofs = [ diridofs ; releasesDofs] ;
neumdofs = (1:(ndofpnode*nnodes))';
neumdofs(diridofs) = [];
% -------------------------------------------------------------



loadFactors     = 0 ;
itersPerTime    = 0 ;
itersPerTimeVec = 0 ;
controlDisps    = 0 ;

timesVec = [ 0] ;

factorescriticos = [] ;

% create velocity and displacements vectors
Ut      = zeros( ndofpnode*nnodes,   1 ) ;  
Udotdott= zeros( ndofpnode*nnodes,   1 ) ;

if exist( 'nonHomogeneousInitialCondU0') ~=0
  for i=1:size(nonHomogeneousInitialCondU0,1)
    dofs= nodes2dofs(nonHomogeneousInitialCondU0(i,1),ndofpnode);
    Ut( dofs (nonHomogeneousInitialCondU0(i,2)))=nonHomogeneousInitialCondU0(i,3);
    
  end
end


if exist( 'nonHomogeneousInitialCondUdot0') ~=0 
  if dynamicAnalysisBoolean == 1
    Udott   = zeros( ndofpnode*nnodes,   1 ) ;  
    for i=1:size(nonHomogeneousInitialCondUdot0,1)
       dofsdot= nodes2dofs(nonHomogeneousInitialCondUdot0(i,1),ndofpnode);
       Udott( dofsdot (nonHomogeneousInitialCondUdot0(i,2)))=nonHomogeneousInitialCondUdot0(i,3);
    end
  else
    error('Velocity initial conditions set for static analysis.');
  end
end


Utm1    = Ut ;  

Utp1    = zeros( ndofpnode*nnodes,   1 ) ;  

FintGt  = zeros( ndofpnode*nnodes,   1 ) ;  

matUts = Ut ;

matNts = zeros(nelems,1) ;

dispsElemsMat = zeros(nelems,2*ndofpnode) ;
for i=1:nelems
  % obtains nodes and dofs of element
  nodeselem = Conec(i,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , ndofpnode ) ;
  dispsElemsMat( i, : ) = Ut(dofselem)' ;
end

Stresst   = zeros(nelems,1) ;
Strainst  = zeros(nelems,1) ;
dsigdepst = zeros(nelems,1) ;

stopTimeIncrBoolean = 0 ;

currLoadFactor  = 0 ;
currTime        = 0 ;
timeIndex       = 1 ;
convDeltau      = zeros(nnodes*ndofpnode,1) ;

stopCritPar = 0;

loadFactors( timeIndex,1) = currLoadFactor ;
controlDisps(timeIndex,1) = Ut(controlDof)*controlDofFactor ;

FintGt = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut ,1 ) ;

factor_crit = 0;
nKeigpos=0;
nKeigneg=0;

if dynamicAnalysisBoolean == 0
  nextLoadFactor  = currLoadFactor + targetLoadFactr / nLoadSteps ;

else 
  deltaT         = numericalMethodParams(2)        ;
  nextLoadFactor = loadFactorsFunc(currTime+deltaT);

end

% stores model data structures
modelCompress

printsOutputScreen

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


% This script declares several matrices and vectors required for the analysis. In this script, the value of important magnitudes, such as internal and external forces, displacements, and velocities are computed for step/time 0.


function [ modelCurrState, modelProperties, BCsCurrState, auxIO, controlDisps, loadFactors, stopTimeIncrBoolean, dispsElemsMat ] ...
  = initialDefinitions( ...
  Conec, nnodes, nodalSprings, ndofpnode, nonHomogeneousInitialCondU0 ...
  , nonHomogeneousInitialCondUdot0, dynamicAnalysisBoolean, controlDofsAndFactors ...
  , secGeomProps, coordsElemsMat, hyperElasParamsMat, numericalMethodParams ...
  , loadFactorsFunc, booleanConsistentMassMat, nodalDamping, booleanScreenOutput ...
  , constantFext, variableFext, userLoadsFilename, stabilityAnalysisBoolean ...
  , problemName, outputDir, nLoadSteps ...
  )

nelems = size(Conec,1) ;

% ----------- fixeddofs and spring matrix computation ---------
[ neumdofs, diridofs, KS] = computeBCDofs(nnodes, Conec, nelems, nodalSprings ) ;
% -------------------------------------------------------------

itersPerTime    = 0 ;
itersPerTimeVec = 0 ;

timesVec = [ 0 ] ;

factorescriticos = [] ;

% create velocity and displacements vectors
Ut      = zeros( ndofpnode*nnodes,   1 ) ;  
Udott   = zeros( ndofpnode*nnodes,   1 ) ;  
Udotdott= zeros( ndofpnode*nnodes,   1 ) ;

if length( nonHomogeneousInitialCondU0 ) > 0
  for i=1:size(nonHomogeneousInitialCondU0,1) % loop over rows of matrix
    dofs= nodes2dofs(nonHomogeneousInitialCondU0(i,1), ndofpnode ) ;
    Ut( dofs (nonHomogeneousInitialCondU0(i,2)))=nonHomogeneousInitialCondU0(i,3);
  end 
end % if nonHomIniCond


if length( nonHomogeneousInitialCondUdot0 ) > 0 
  if dynamicAnalysisBoolean == 1
    for i=1:size(nonHomogeneousInitialCondUdot0,1)
      dofs = nodes2dofs(nonHomogeneousInitialCondUdot0(i,1), ndofpnode ) ;
      Udott( dofs(nonHomogeneousInitialCondUdot0(i,2)))=nonHomogeneousInitialCondUdot0(i,3);
    end
  else
    error('Velocity initial conditions set for static analysis.');
  end
end


% computation of initial acceleration for some cases
% --------------------------------------------------- 


Utp1   = zeros( ndofpnode*nnodes,   1 ) ;  

Fintt  = zeros( ndofpnode*nnodes,   1 ) ;  

matUts = Ut ;


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

stopCritPar = 0 ;



% --- load factors and control displacements ---
loadFactors     = 0 ;
loadFactors( timeIndex,1) = currLoadFactor ;

controlDisps    = 0 ;
controlDisps(timeIndex, :) = Ut( controlDofsAndFactors(:,1) ) ...
                              .* controlDofsAndFactors(:,2) ;
% ----------------------------------------------

if dynamicAnalysisBoolean == 0,
  nextLoadFactor  = currLoadFactor + numericalMethodParams(5) / nLoadSteps ;

else 
  deltaT         = numericalMethodParams(2)         ;
  nextLoadFactor = loadFactorsFunc(currTime+deltaT) ;
end


% --- initial force vectors ---
if length(nodalDamping)>0
  dampingMat = speye( nnodes*6, nnodes*6 ) * nodalDamping   ;
else
  dampingMat = [] ;
end

[ fs, Strainst, Stresst ] = assembler ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut , dynamicAnalysisBoolean , 1, Udotdott, booleanConsistentMassMat ) ;

Fintt = fs{1} ;
Fmast = fs{2} ;

  systemDeltauMatrix          = computeMatrix( Conec, secGeomProps, coordsElemsMat, ...
    hyperElasParamsMat, KS, Ut, neumdofs, numericalMethodParams, [], ...
    dampingMat, booleanConsistentMassMat, Udotdott );

% ----------------------------



% --- initial tangent matrices ---
%~ [ mats ] = assembler ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut, dynamicAnalysisBoolean, 2, Udotdott, booleanConsistentMassMat ) ;


%~ systemDeltauMatrix = mats{1} 

%~ stop
%~ if dynamicAnalysisBoolean == 1,

  %~ massMat    = mats{2} ;

  %~ % --- computation of initial Udotdott for truss elements only!!!
  %~ Fext = computeFext( constantFext, variableFext, loadFactors(1), userLoadsFilename ) ;

  %~ Udotdott (neumdofs) = massMat( neumdofs, neumdofs ) \ ( Fext(neumdofs) -Fintt( neumdofs ) ) ;  

%~ else
  %~ dampingMat = [] ;
  %~ massMat    = [] ;
%~ end

factorCrit = 0 ;
nKeigpos   = 0 ;
nKeigneg   = 0 ;
% ----------------------------


%~ indselems12 = find( ( Conec(:,7) == 1) | ( Conec(:,7) == 2) ) ;
Areas = secGeomProps(Conec(:,6),1) ;
currentNormalForces = Stresst(:,1) .* Areas ;

matNts = currentNormalForces ;

% stores model data structures
modelCompress

% --- prints headers and time0 values ---
printSolverOutput( outputDir, problemName, timeIndex, 0 ) ;

fprintf( '|-------------------------------------------------|\n' ) ;
fprintf( '| TimeSteps progress: 1|                   |%4i  |\n                        ', nLoadSteps)

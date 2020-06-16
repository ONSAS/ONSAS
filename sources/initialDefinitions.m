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
  Conec, nNodes, nodalSprings, nonHomogeneousInitialCondU0 ...
  , nonHomogeneousInitialCondUdot0, controlDofsAndFactors ...
  , crossSecsParams, coordsElemsMat, materialsParamsMat, numericalMethodParams ...
  , loadFactorsFunc, booleanConsistentMassMat, nodalDamping, booleanScreenOutput ...
  , constantFext, variableFext, userLoadsFilename, stabilityAnalysisBoolean ...
  , problemName, outputDir ...
  )

nElems    = size(Conec,1) ;

% ----------- fixeddofs and spring matrix computation ---------
[ neumdofs, diridofs, KS] = computeBCDofs( nNodes, Conec, nElems, nodalSprings ) ;
% -------------------------------------------------------------

% create velocity and displacements vectors
Ut       = zeros( 6*nNodes,   1 ) ;  
Udott    = zeros( 6*nNodes,   1 ) ;  
Udotdott = zeros( 6*nNodes,   1 ) ;

if length( nonHomogeneousInitialCondU0 ) > 0
  for i=1:size(nonHomogeneousInitialCondU0,1) % loop over rows of matrix
    dofs= nodes2dofs(nonHomogeneousInitialCondU0(i,1), 6 ) ;
    Ut( dofs ( nonHomogeneousInitialCondU0(i,2))) = ...
      nonHomogeneousInitialCondU0(i,3);
  end 
end % if nonHomIniCond

dispsElemsMat = zeros(nElems,2*6) ;
for i=1:nElems
  % obtains nodes and dofs of element
  nodeselem = Conec(i,1:2)' ;
  dofselem  = nodes2dofs( nodeselem , 6 ) ;
  dispsElemsMat( i, : ) = Ut(dofselem)' ;
end


if length( nonHomogeneousInitialCondUdot0 ) > 0
  if numericalMethodParams(1) >= 3
    for i=1:size(nonHomogeneousInitialCondUdot0, 1)
      dofs = nodes2dofs( nonHomogeneousInitialCondUdot0(i, 1), 6 ) ;
      Udott( dofs( nonHomogeneousInitialCondUdot0(i,2))) = ...
        nonHomogeneousInitialCondUdot0(i,3);
    end
  else
    error( ' velocity initial conditions set for a static analysis method');  
  end
end


% computation of initial acceleration for some cases
% --------------------------------------------------- 

stopTimeIncrBoolean = 0 ;

currTime        = 0 ;
timeIndex       = 1 ;
convDeltau      = zeros(nNodes*6,1) ;

timeStepIters    = 0 ;
timeStepStopCrit = 0 ;


% --- load factors and control displacements ---
currLoadFactor           = loadFactorsFunc( currTime ) ;
loadFactors              = currLoadFactor ; % initialize

controlDisps    = 0 ;
controlDisps(timeIndex, :) = Ut( controlDofsAndFactors(:,1) ) ...
                              .* controlDofsAndFactors(:,2) ;
% ----------------------------------------------

[ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
  = extractMethodParams( numericalMethodParams ) ;
           
nextLoadFactor = loadFactorsFunc ( currTime + deltaT ) ;

% --- initial force vectors ---
dampingMat = speye( nNodes*6, nNodes*6 ) * nodalDamping   ;

[ fs, Strainst, Stresst ] = assembler ( ...
  Conec, crossSecsParams, coordsElemsMat, materialsParamsMat, KS, Ut, 1, Udotdott, booleanConsistentMassMat ) ;

stop
Fintt = fs{1} ;
Fmast = fs{2} ;

systemDeltauMatrix     = computeMatrix( ...
  Conec, crossSecsParams, coordsElemsMat, materialsParamsMat, KS, Ut, ...
  neumdofs, numericalMethodParams, dampingMat, booleanConsistentMassMat, Udotdott );

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

[ nKeigpos, nKeigneg ] = stabilityAnalysis ( [], systemDeltauMatrix, currLoadFactor, nextLoadFactor )
% ----------------------------


%~ indselems12 = find( ( Conec(:,7) == 1) | ( Conec(:,7) == 2) ) ;
Areas = crossSecsParams (Conec(:,6),1) ;
currentNormalForces = Stresst(:,1) .* Areas ;

matNts = currentNormalForces ;

% stores model data structures
modelCompress

% --- prints headers and time0 values ---
printSolverOutput( outputDir, problemName, timeIndex, 0 ) ;

fprintf( '|-------------------------------------------------|\n' ) ;
fprintf( '| TimeSteps progress: 1|                   |%4i  |\n                        ', nLoadSteps)

printSolverOutput( ...
  outputDir, problemName, timeIndex, [ 2 currLoadFactor 0 0 nKeigpos nKeigneg ] ) ;


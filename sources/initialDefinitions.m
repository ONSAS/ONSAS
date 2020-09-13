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


% This script declares several matrices and vectors required for the analysis. In this script, the value of important magnitudes, such as internal and external forces, displacements, and velocities are computed for step/time 0.


function [ modelCurrSol, modelProperties, BCsData ] ...
  = initialDefinitions( ...
  Conec, nNodes, nodalSprings, nonHomogeneousInitialCondU0 ...
  , nonHomogeneousInitialCondUdot0 ...
  , crossSecsParamsMat, coordsElemsMat, materialsParamsMat, numericalMethodParams ...
  , loadFactorsFunc, nodalDispDamping, booleanScreenOutput ...
  , constantFext, variableFext, userLoadsFilename, stabilityAnalysisBoolean ...
  , problemName, outputDir, finalTime, elementsParamsMat ...
  )

nElems    = size(Conec, 1 ) ;

% ----------- fixeddofs and spring matrix computation ---------
[ neumdofs, diridofs, KS] = computeBCDofs( nNodes, Conec, nElems, nodalSprings, elementsParamsMat ) ;
% -------------------------------------------------------------

% create velocity and displacements vectors
U       = zeros( 6*nNodes,   1 ) ;  
Udot    = zeros( 6*nNodes,   1 ) ;  
Udotdot = zeros( 6*nNodes,   1 ) ;

% adds non homogeneous initial conditions 
if length( nonHomogeneousInitialCondU0 ) > 0
  for i = 1 : size( nonHomogeneousInitialCondU0, 1 ) % loop over rows of matrix
    dofs = nodes2dofs(nonHomogeneousInitialCondU0(i, 1 ), 6 ) ;
    U( dofs ( nonHomogeneousInitialCondU0 (i, 2 ) ) ) = ...
      nonHomogeneousInitialCondU0 ( i, 3 ) ;
  end 
end % if nonHomIniCond

if length( nonHomogeneousInitialCondUdot0 ) > 0
  if numericalMethodParams(1) >= 3
    for i=1:size(nonHomogeneousInitialCondUdot0, 1)
      dofs = nodes2dofs( nonHomogeneousInitialCondUdot0(i, 1), 6 ) ;
      Udot( dofs( nonHomogeneousInitialCondUdot0(i, 2 ))) = ...
        nonHomogeneousInitialCondUdot0(i, 3 );
    end
  else
    warning(' velocity initial conditions set for a static analysis method' ) ;  
  end
end


%~ dispsElemsMat = zeros( nElems, 4*6) ;
%~ for i=1:nElems
  %~ % obtains nodes and dofs of element
  %~ nodeselem = Conec(i,1:2)' ;
  %~ dofselem  = nodes2dofs( nodeselem , 6 ) ;
  %~ dispsElemsMat( i, : ) = U(dofselem)' ;
%~ end

% computation of initial acceleration for some cases
% --------------------------------------------------- 

currTime        = 0 ;
timeIndex       = 1 ;
convDeltau      = zeros( nNodes*6, 1 ) ;

timeStepIters    = 0 ;
timeStepStopCrit = 0 ;


% --- load factors and control displacements ---
currLoadFactor   = loadFactorsFunc( currTime ) ;
loadFactors      = currLoadFactor ; % initialize

% ----------------------------------------------

[ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
  = extractMethodParams( numericalMethodParams ) ;
           
nextLoadFactor = loadFactorsFunc ( currTime + deltaT ) ;

%~ fs = assembler ( ...
  %~ Conec, crossSecsParams, coordsElemsMat, materialsParamsMat, KS, U, 1, Udot, ...
  %~ Udotdot, nodalDispDamping, solutionMethod, booleanConsistentMassMat, booleanCSTangs ) ;

%~ Fint = fs{1} ;  Fvis =  fs{2};  Fmas = fs{3} ;

systemDeltauMatrix     = computeMatrix( ...
  Conec, crossSecsParamsMat, coordsElemsMat, materialsParamsMat, KS, U, ...
  neumdofs, numericalMethodParams, nodalDispDamping, Udot, Udotdot, ...
  elementsParamsMat ) ;

Stress = assembler ( ...
  Conec, crossSecsParamsMat, coordsElemsMat, materialsParamsMat, KS, U, 3, Udot, ...
  Udotdot, nodalDispDamping, solutionMethod, elementsParamsMat ) ;
% ----------------------------

factorCrit = 0 ; nKeigneg = 0 ; nKeigpos = 0 ;
%%%%%%%%%%%%%%%%%
%~ [ nKeigpos, nKeigneg ] = stabilityAnalysis ( ...
  %~ [], systemDeltauMatrix, currLoadFactor, nextLoadFactor ) ;
%%%%%%%%%%%%%%%%%
% ----------------------------


% stores model data structures
modelCompress

% --- prints headers and time0 values ---
printSolverOutput( outputDir, problemName, timeIndex, 0 ) ;

fprintf( '|-------------------------------------------------|\n' ) ;
fprintf( '| TimeSteps progress: 1|                   |%4i  |\n                        ', nLoadSteps)

printSolverOutput( ...
  outputDir, problemName, timeIndex, [ 2 currLoadFactor 0 0 nKeigpos nKeigneg ] ) ;


% --- initial tangent matrices ---
%~ [ mats ] = assembler ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Ut, dynamicAnalysisBoolean, 2, Udotdott, booleanConsistentMassMat ) ;

%~ systemDeltauMatrix = mats{1} 

%~ stop
%~ if dynamicAnalysisBoolean == 1,

  %~ massMat    = mats{2} ;

  %~ % --- computation of initial Udotdott for truss elements only!!!
  %~ Fext = computeFext( constantFext, variableFext, loadFactors(1), userLoadsFilename ) ;

  %~ Udotdott (neumdofs) = massMat( neumdofs, neumdofs ) \ ( Fext(neumdofs) -Fintt( neumdofs ) ) ;  

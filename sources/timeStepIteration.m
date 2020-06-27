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

%% This functions performs the iteration for the computations of the state in
% the next time step using the numerical method and parameters provided by the
% user.

function  [ modelCurrSol, BCsData ] = timeStepIteration( modelCurrSol, BCsData, modelProperties ) ;

startPreviousVels = 0 ; % 0 recommended

% ----   extracts variables  ----
modelExtract

% -----   pre-iteration definitions     ----------
nElems     = size( Conec, 1 ) ; ndofpnode = 6;

[ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
      = extractMethodParams( numericalMethodParams ) ;
% --------------------------------------------------------------------


% --------------------------------------------------------------------
% ----       iteration in displacements or load-displacements     ----
% --------------------------------------------------------------------

KTtred = systemDeltauMatrix ;

% assign time t
Ut = U ; Udott = Udot ; Udotdott = Udotdot ;

%~ Fintt = Fint ; Fmast = Fmas ; Fvist = Fvis ;

% --- start iteration with previous displacements ---
Utp1k       = Ut       ;   % initial guess

if startPreviousVels == 1
  Udottp1k    = Udott    ;
  Udotdottp1k = Udotdott ;
else
  [ Udottp1k, Udotdottp1k ] = updateTime( ...
    Ut, Udott, Udotdott, Utp1k, numericalMethodParams, currTime ) ;
end

if solutionMethod == 2
  nextLoadFactor = currLoadFactor ; % initial guess for next load factor
end

% --- compute RHS for initial guess ---
[ systemDeltauRHS, FextG ]  = computeRHS( ...
  Conec, crossSecsParams, coordsElemsMat, ...
  materialsParamsMat, KS, constantFext, variableFext, userLoadsFilename, ...
  currLoadFactor, nextLoadFactor, numericalMethodParams, neumdofs, nodalDispDamping, ...
  booleanConsistentMassMat, booleanCSTangs, ...
  Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k ) ;
% ---------------------------------------------------

booleanConverged = 0                              ;
dispIters        = 0                              ;
currDeltau       = zeros( length( neumdofs ), 1 ) ;

while  booleanConverged == 0
  dispIters = dispIters + 1 ;

  % --- solve system ---
  [ deltaured, nextLoadFactor ] = computeDeltaU ( systemDeltauMatrix, systemDeltauRHS, dispIters, convDeltau(neumdofs), numericalMethodParams, nextLoadFactor , currDeltau ) ;
  % ---------------------------------------------------

  % --- updates: model variables and computes internal forces ---
  [Utp1k, currDeltau] = updateUiter(Utp1k, deltaured, neumdofs, solutionMethod, currDeltau ) ;

  % --- update next time magnitudes ---
  [ Udottp1k, Udotdottp1k, nextTime ] = updateTime( ...
    Ut, Udott, Udotdott, Utp1k, numericalMethodParams, currTime ) ;
  % ---------------------------------------------------
  
  % --- system matrix ---
  systemDeltauMatrix          = computeMatrix( Conec, crossSecsParams, coordsElemsMat, ...
    materialsParamsMat, KS, Utp1k, neumdofs, numericalMethodParams, ...
    nodalDispDamping, booleanConsistentMassMat, Udott, Udotdott, booleanCSTangs ) ;
  % ---------------------------------------------------

  % --- new rhs ---
  [ systemDeltauRHS, FextG ]  = computeRHS( Conec, crossSecsParams, coordsElemsMat, ...
    materialsParamsMat, KS, constantFext, variableFext, ...
    userLoadsFilename, currLoadFactor, nextLoadFactor, numericalMethodParams, ...
    neumdofs, nodalDispDamping, booleanConsistentMassMat, booleanCSTangs, ...
    Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k ) ;
  % ---------------------------------------------------

  % --- check convergence ---
  [booleanConverged, stopCritPar, deltaErrLoad ] = convergenceTest( numericalMethodParams, [], FextG(neumdofs), deltaured, Utp1k(neumdofs), dispIters, [], systemDeltauRHS ) ;
  % ---------------------------------------------------

  % --- prints iteration info in file ---
  printSolverOutput( outputDir, problemName, timeIndex, [ 1 dispIters deltaErrLoad norm(deltaured) ] ) ;

end % iteration while
% --------------------------------------------------------------------

Utp1       = Utp1k ;
Udottp1    = Udottp1k ;
Udotdottp1 = Udotdottp1k ;

% computes KTred at converged Uk
KTtp1red = systemDeltauMatrix ;

% --------------------------------------------------------------------

Stresstp1 = assembler ( Conec, crossSecsParams, coordsElemsMat, materialsParamsMat, KS, Utp1, 3, Udottp1, Udotdottp1, nodalDispDamping, solutionMethod, booleanConsistentMassMat, booleanCSTangs ) ;

if stabilityAnalysisBoolean == 1
  [ nKeigpos, nKeigneg, factorCrit ] = stabilityAnalysis ( KTtred, KTtp1red, currLoadFactor, nextLoadFactor ) ;
else
  [ nKeigpos, nKeigneg ] = stabilityAnalysis ( KTtred, KTtp1red, currLoadFactor, nextLoadFactor ) ;
  factorCrit = 0;
end

% prints iteration info in file
printSolverOutput( ...
  outputDir, problemName, timeIndex+1, [ 2 nextLoadFactor dispIters stopCritPar nKeigpos nKeigneg ] ) ;


% --- stores next step values ---
U          = Utp1 ;
Udot       = Udottp1  ;
Udotdot    = Udotdottp1 ;
convDeltau = Utp1 - Ut ;
%
Stress     = Stresstp1 ;

timeIndex  = timeIndex + 1 ;

currTime   = nextTime ;
if solutionMethod == 2
  currTime = nextLoadFactor ;
end

timeStepStopCrit = stopCritPar ;
timeStepIters = dispIters ;

modelCompress
% -------------------------------------


% ==============================================================================
%
% ==============================================================================
function [ Udottp1, Udotdottp1, nextTime ] = updateTime(Ut,Udott,Udotdott, Uk, numericalMethodParams, currTime )

  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
      = extractMethodParams( numericalMethodParams ) ;

  nextTime   = currTime + deltaT                          ;

if solutionMethod == 3 || solutionMethod == 4

  if solutionMethod == 4
      deltaNW = (1-2*alphaHHT)/2 ;
      AlphaNW = (1-alphaHHT^2)/4 ;
  end
  
  [a0NM, a1NM, a2NM, a3NM, a4NM, a5NM, a6NM, a7NM ] = coefsNM( AlphaNW, deltaNW, deltaT ) ;
  
  Udotdottp1 = a0NM*(Uk-Ut) - a2NM*Udott - a3NM*Udotdott;
  Udottp1    = Udott + a6NM*Udotdott + a7NM*Udotdottp1    ;
  
else
  Udotdottp1 = [] ;
  Udottp1    = [] ;
end


% ==============================================================================
%
% ==============================================================================
function [Uk, currDeltau] = updateUiter(Uk, deltaured, neumdofs, solutionMethod, currDeltau ) 

  Uk ( neumdofs ) = Uk(neumdofs ) + deltaured ;

  currDeltau      = currDeltau    + deltaured ;

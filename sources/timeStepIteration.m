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

function  [ modelCurrState, BCsCurrState, auxIO ]  = timeStepIteration( modelCurrState, BCsCurrState, modelProperties, auxIO ) ;

% ----   extracts variables  ----
modelExtract

materialsParamsMat
stop

% -----   pre-iteration definitions     ----------
nelems     = size(Conec,1) ; ndofpnode = 6;

[ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
      = extractMethodParams( numericalMethodParams ) ;

% --------------------------------------------------------------------
  

% --------------------------------------------------------------------
% ----       iteration in displacements or load-displacements     ----
% --------------------------------------------------------------------

KTt = systemDeltauMatrix ;

% --- start iteration with previous displacements ---

Utp1k       = Ut     ;   % initial guess
Finttp1k    = Fintt  ;
Fmastp1k    = Fmast  ;
Udottp1k    = Udotdott ;
Udotdottp1k = Udotdott ;

if solutionMethod == 2
  nextLoadFactor = currLoadFactor ; % initial guess
end

% --- compute RHS for initial guess ---
[ systemDeltauRHS, FextG ]  = computeRHS( Conec, crossSecsParams, coordsElemsMat, materialsParams, KS, Utp1k, constantFext, variableFext, userLoadsFilename, currLoadFactor, nextLoadFactor, numericalMethodParams, neumdofs, Finttp1k, dampingMat, Ut, Udott, Udotdott, Fintt, Fmast ) ;
% ---------------------------------------------------


booleanConverged = 0 ;
dispIters        = 0 ;
currDeltau       = zeros( length(neumdofs), 1 ) ;

while  booleanConverged == 0
  dispIters = dispIters + 1 ;

%~ systemDeltauMatrix
%~ systemDeltauRHS
%~ st
  % --- solve system ---
  [ deltaured, nextLoadFactor ] = computeDeltaU ( systemDeltauMatrix, systemDeltauRHS, dispIters, convDeltau(neumdofs), numericalMethodParams, nextLoadFactor , currDeltau ) ;
  % ---------------------------------------------------

materialsParams
  % --- updates: model variables and computes internal forces ---
  [Utp1k, currDeltau] = updateUiter(Utp1k, deltaured, neumdofs, solutionMethod, currDeltau ) ;

  % --- update next time magnitudes ---
  [ Utp1k, Udottp1k, Udotdottp1k, nextTime ] = updateTime( ...
    Ut, Udott, Udotdott, Fintt, Utp1k, Finttp1k, numericalMethodParams, currTime ) ;

  Fs = assembler ( Conec, crossSecsParams, coordsElemsMat, materialsParams, KS, Utp1k, 1, Udotdottp1k, booleanConsistentMassMat ) ;
  % ---------------------------------------------------
  
  % --- system matrix ---
  systemDeltauMatrix          = computeMatrix( Conec, crossSecsParams, coordsElemsMat, ...
    materialsParams, KS, Utp1k, neumdofs, numericalMethodParams, ...
    dampingMat, booleanConsistentMassMat, Udotdott );
  % ---------------------------------------------------

  % --- new rhs ---
  [ systemDeltauRHS, FextG ]  = computeRHS( Conec, crossSecsParams, coordsElemsMat, ...
    materialsParams, KS, Utp1k, constantFext, variableFext, ...
    userLoadsFilename, currLoadFactor, nextLoadFactor, numericalMethodParams, ...
    neumdofs, Fs{1}, massMat, dampingMat, Ut, Udott, Udotdott, Fintt ) ;
  % ---------------------------------------------------

  % --- check convergence ---
  [booleanConverged, stopCritPar, deltaErrLoad ] = convergenceTest( numericalMethodParams, FintGk(neumdofs), FextG(neumdofs), deltaured, Uk(neumdofs), dispIter, [], systemDeltauRHS ) ;
  % ---------------------------------------------------

  % --- prints iteration info in file ---
  printSolverOutput( outputDir, problemName, timeIndex, [ 1 dispIter deltaErrLoad norm(deltaured) ] ) ;

end % iteration while
% --------------------------------------------------------------------
% --------------------------------------------------------------------

% computes KTred at converged Uk
[ KTtp1 ] = assembler( Conec, crossSecsParams, coordsElemsMat, materialsParams, KS, Uk, [], 2, Udotdott, booleanConsistentMassMat ) ;

factor_crit = 0;

[FintGk, Strainsk, Stressk ] = assembler ( Conec, crossSecsParams, coordsElemsMat, materialsParams, KS, Uk, [], 1, Udotdott, booleanConsistentMassMat ) ;

if stabilityAnalysisBoolean == 1
  [ factor_crit, nKeigpos, nKeigneg ] = stabilityAnalysis ( KTtm1( neumdofs, neumdofs ), KTt( neumdofs, neumdofs ), currLoadFactor, nextLoadFactor ) ;
else
  factor_crit = 0; nKeigpos=0; nKeigneg=0;
end


% prints iteration info in file
printSolverOutput( ...
  outputDir, problemName, timeIndex+1, [ 2 nextLoadFactor dispIter stopCritPar nKeigpos nKeigneg ] ) ;

% --- stores next step as Ut and Ft ---

% -------------------------------------
currTime  = nextTime ;
if solutionMethod == 2
  currTime = nextLoadFactor ;
end

Ut        = Utp1 ;
FintGt    = FintGtp1 ;
Udott     = Udottp1 ;
Udotdott  = Udotdottp1 ;
timeIndex = timeIndex + 1 ;


modelCompress




% ==============================================================================
%
% ==============================================================================
function [ Utp1, Udottp1, Udotdottp1, nextTime ] = updateTime(Ut,Udott,Udotdott, FintGt, Uk, FintGk, numericalMethodParams, currTime )

  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
      = extractMethodParams( numericalMethodParams ) ;

  Utp1       = Uk                                         ;
  nextTime   = currTime + deltaT                          ;

if solutionMethod == 3 || solutionMethod == 4

  if solutionMethod == 4
      deltaNW = (1-2*alphaHHT)/2 ;
      AlphaNW = (1-alphaHHT^2)/4 ;
  end
  
  [a0NM, a1NM, a2NM, a3NM, a4NM, a5NM, a6NM, a7NM ] = coefsNM( AlphaNW, deltaNW, deltaT ) ;
  
  Udotdottp1 = a0NM*(Utp1-Ut) - a2NM*Udott - a3NM*Udotdott;
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

  if solutionMethod == 2
    currDeltau      = currDeltau    + deltaured ;
  end

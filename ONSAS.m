% ==============================================================================
% --------     ONSAS: an Open Non-linear Structural Analysis System     --------
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
% ==============================================================================

Taux = cputime() ;

% ==============================================================================
% ----------------------------     Input       ---------------------------------
ONSASversion ='0.1.10';

if isunix, dirSep = '/'; else, dirSep = '\'; end
addpath( [ pwd dirSep 'sources'] ) ;
addpath( [ pwd dirSep 'user'   ] ) ;

% verifies the definition of input variables and sets default values
inputVarsVerification

% ==============================================================================

% ==============================================================================
% ----------------------------    Analysis     ---------------------------------


% Initial computations: sets initial state.
[ modelCurrSol, modelProperties, BCsData ] =  initialDefinitions( ...
  Conec, nNodes, nodalSprings, nonHomogeneousInitialCondU0 ...
  , nonHomogeneousInitialCondUdot0 ...
  , crossSecsParamsMat, coordsElemsMat, materialsParamsMat, numericalMethodParams ...
  , loadFactorsFunc, booleanConsistentMassMat, nodalDispDamping, booleanScreenOutput ...
  , constantFext, variableFext, userLoadsFilename, stabilityAnalysisBoolean ...
  , problemName, outputDir, booleanCSTangs, finalTime, elementsParamsMat  ) ;

% --- increment step analysis ---
stopTimeIncrBoolean = 0 ;
while ( stopTimeIncrBoolean == 0 )

  % -----   computes the model state at the next load/time step   -----
  [ modelNextSol, BCsData ] = timeStepIteration ...
  ( modelCurrSol, BCsData, modelProperties );

  % --- checks stopping criteria and stores model state
  storeAndPlotResults

  % --- update state data --- 
  modelCurrSol           = modelNextSol ;
  
  BCsData.currLoadFactor = BCsData.nextLoadFactor                   ;
  BCsData.nextLoadFactor = loadFactorsFunc( modelCurrSol.currTime + deltaT ) ;

  % ----   evals stop time incr crit    --------
  stopTimeIncrBoolean = modelCurrSol.currTime >= finalTime ;
    
end

fprintf( '\n| end time-step%4i - (max,avg) iters: (%3i,%5.2f) | \n ',...
  modelCurrSol.timeIndex, max( itersPerTimeVec ) , mean( itersPerTimeVec(2:end) ) );
 

% if analytical solution is provided, numerical results are validated. 
if analyticSolFlag > 0
  [numericalVals, analyticVals] = analyticSolVerif ...
    ( analytSol, analyticFunc, loadFactors, controlDisps, timesVec, ...
    analyticCheckTolerance, analyticSolFlag, problemName, printFlag, outputDir, plotParamsVector );
end
% ==============================================================================


% ==============================================================================
% ----------------------------     Output      ---------------------------------

% plots and/or visualization files are generated
%~ if plotParamsVector(1) > 0
  %~ outputPlots( matUs, coordsElemsMat, plotParamsVector, ...
    %~ Conec, Nodes, constantFext, variableFext, controlDisps, ...
    %~ deformedScaleFactor, printFlag, ...
    %~ outputDir, problemName, loadFactors, sectPar, ...
    %~ deltaT, cellStress, plotsViewAxis, booleanScreenOutput ) ;
%~ end

% report with results is generated
if reportBoolean
  outputReport
end
totalTime = cputime()-Taux ;
fprintf([ '|-------------------------------------------------|\n'])
fprintf(  '|  ONSAS finished in: %7.1e seconds /%5.2f mins |\n', totalTime, totalTime/60 )
fprintf([ '|=================================================|\n\n\n'])

% Copyright 2023, ONSAS Authors (see documentation)
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
% 
% This functions performs the iteration for the computations of the state in
% the next time step using the numerical method and parameters provided by the
% user.
function modelNextSol = timeStepIteration( modelCurrSol, modelProperties, BCsData ) ;

% assign current time (t) variables
% ---------------------------------
Ut         = modelCurrSol.U ; Udott = modelCurrSol.Udot ; Udotdott = modelCurrSol.Udotdot ;
convDeltau = modelCurrSol.convDeltau ;
currLoadFactorsVals = modelCurrSol.currLoadFactorsVals ;

% stability analysis
% -----------------------------------------------------------
stabilityAnalysisFlag = modelProperties.analysisSettings.stabilityAnalysisFlag ;
if ~(stabilityAnalysisFlag==0)
  error(' stability analysis pending: see issue https://github.com/ONSAS/ONSAS/issues/351');  
  % reduced tangent matrix of previous time for nonlinear buckling analysis
  KTtred = modelCurrSol.systemDeltauMatrix ;
end

% update time and set candidate displacements and derivatives
% -----------------------------------------------------------
if isempty( modelProperties.analysisSettings.Utp10 )
  Utp1k = Ut ;
else
  error('Add case for several times.')
end

[ Udottp1k, Udotdottp1k, nextTime ] = updateTime( ...
  Ut, Udott, Udotdott, Utp1k, modelProperties.analysisSettings, modelCurrSol.currTime ) ;

% compute RHS for initial guess Utp1 and in next time step
% --------------------------------------------------------
if strcmp( modelProperties.analysisSettings.methodName, 'arcLength') == 1
  nextLoadFactorsVals = currLoadFactorsVals ;
  args = argsAL(modelProperties.analysisSettings, length(convDeltau), BCsData.neumDofs, modelCurrSol.timeIndex) ;
else
  args = [] ;
  nextLoadFactorsVals = [] ;
end

% current system variables
% ----------------------
systemDeltauRHS    = modelCurrSol.systemDeltauRHS    ;
systemDeltauMatrix = modelCurrSol.systemDeltauMatrix ;
previousStateCell  = modelCurrSol.previousStateCell  ;

% --- assemble system of equations ---
[ systemDeltauMatrix, systemDeltauRHS, FextG, ~, nextLoadFactorsVals ] = system_assembler( modelProperties, BCsData, Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k, nextTime, nextLoadFactorsVals, previousStateCell ) ;

booleanConverged = false ;
dispIters        = 0     ;
currDeltau       = zeros( length( BCsData.neumDofs ), 1 ) ;

while  booleanConverged == 0

  %fprintf(' ============== new iteration ====================\n')
  dispIters = dispIters + 1 ;

  % solve system
  [ deltaured, nextLoadFactorsVals ] = computeDeltaU( systemDeltauMatrix, systemDeltauRHS, dispIters, convDeltau(BCsData.neumDofs), modelProperties.analysisSettings, nextLoadFactorsVals , currDeltau, modelCurrSol.timeIndex, BCsData.neumDofs, args ) ;

  % updates: model variables and computes internal forces ---
  [Utp1k, currDeltau] = updateUiter(Utp1k, deltaured, BCsData.neumDofs, currDeltau ) ;

  % --- update next time magnitudes ---
  [ Udottp1k, Udotdottp1k, nextTime ] = updateTime( ...
    Ut, Udott, Udotdott, Utp1k, modelProperties.analysisSettings, modelCurrSol.currTime ) ;

  % --- assemble system of equations ---
  [ systemDeltauMatrix, systemDeltauRHS, FextG, ~, nextLoadFactorsVals, fnorms ] = system_assembler( modelProperties, BCsData, Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k, nextTime, nextLoadFactorsVals, previousStateCell ) ;

  % --- check convergence ---
  [ booleanConverged, stopCritPar, deltaErrLoad, normFext ] = convergenceTest( modelProperties.analysisSettings, FextG(BCsData.neumDofs), deltaured, Utp1k(BCsData.neumDofs), dispIters, systemDeltauRHS(:,1) ) ;
  % ---------------------------------------------------

  %
% paso tiempo     iters     norm rhs  norm fext norm fint   vis = mas = fs aero   ther  
  fnorms = [ modelCurrSol.timeIndex; dispIters; deltaErrLoad; normFext; fnorms; modelCurrSol.currTime  ] ;

  

  % --- prints iteration info in file ---
  printSolverOutput( modelProperties.outputDir, modelProperties.problemName, [ 1 norm(nextLoadFactorsVals) dispIters deltaErrLoad norm(deltaured) ], fnorms ) ;

  # stop

end % iteration while
% --------------------------------------------------------------------

Utp1       = Utp1k ;
Udottp1    = Udottp1k ;
Udotdottp1 = Udotdottp1k ;

% computes KTred at converged Uk
KTtp1red = systemDeltauMatrix ;

% compute stress at converged state
[~, Stresstp1, ~, matFint, strain_vec, acum_plas_strain_vec ] = assembler ( modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData.KS, Utp1, Udottp1, Udotdottp1, modelProperties.analysisSettings, [ 0 1 0 1 ], modelProperties.nodalDispDamping, nextTime, previousStateCell ) ;

printSolverOutput( modelProperties.outputDir, modelProperties.problemName, [ 2 (modelCurrSol.timeIndex)+1 nextTime dispIters stopCritPar ] ,[]) ;

if stabilityAnalysisFlag == 2
  [ nKeigpos, nKeigneg, factorCrit ] = stabilityAnalysis ( KTtred, KTtp1red, currLoadFactor, nextLoadFactor ) ;
elseif stabilityAnalysisFlag == 1
  [ nKeigpos, nKeigneg ] = stabilityAnalysis ( KTtred, KTtp1red, currLoadFactor, nextLoadFactor ) ;
  factorCrit = 0;
else
  nKeigpos = 0;  nKeigneg = 0; factorCrit = 0 ;
end

% --- stores next step values ---
U          = Utp1 ;
Udot       = Udottp1  ;
Udotdot    = Udotdottp1 ;
convDeltau = Utp1 - Ut ;
%
Stress     = Stresstp1 ;

timeIndex  = modelCurrSol.timeIndex + 1 ;

currTime   = nextTime ;

timeStepStopCrit = stopCritPar ;
timeStepIters = dispIters ;

for i = 1:size(Stress,1)
	previousStateCell(i,1) = {Stress(i,:)} ;
end

previousStateCell(:,2) = strain_vec ;
previousStateCell(:,3) = acum_plas_strain_vec ;

modelNextSol = construct_modelSol( timeIndex, currTime, U , Udot, ...
                                   Udotdot, Stress, convDeltau, ...
                                   nextLoadFactorsVals, systemDeltauMatrix, ...
                                   systemDeltauRHS, timeStepStopCrit, timeStepIters, matFint, previousStateCell ) ;

% ==============================================================================
% ==============================================================================
function [ Udottp1, Udotdottp1, nextTime ] = updateTime(Ut, Udott, Udotdott, Uk, analysisSettings, currTime )

  nextTime   = currTime + analysisSettings.deltaT                     ;

  if strcmp( analysisSettings.methodName, 'newmark') || strcmp( analysisSettings.methodName, 'alphaHHT' )

    deltaT = analysisSettings.deltaT ;

    if strcmp( analysisSettings.methodName, 'alphaHHT' )
      alphaHHT = analysisSettings.alphaHHT ;
      deltaNM  = (1-2*alphaHHT)/2 ;
      alphaNM  = (1-alphaHHT)^2/4 ;
    else
      deltaNM  = analysisSettings.deltaNM ;
      alphaNM  = analysisSettings.alphaNM ;
    end

    Udotdottp1 = 1.0/( alphaNM * (deltaT)^2 ) * ( Uk - Ut ) - 1.0/( alphaNM * deltaT ) * Udott - ( 1.0/ ( alphaNM * 2 ) - 1 ) * Udotdott ;

    Udottp1    = Udott + ( ( 1 - deltaNM ) * Udotdott + deltaNM * Udotdottp1 ) * deltaT    ;

  else
    Udotdottp1 = Udotdott ;
    Udottp1    = Udott ;
  end

% ==============================================================================
% update Uiter
% ==============================================================================

function [Uk, currDeltau] = updateUiter(Uk, deltaured, neumdofs, currDeltau )
  Uk( neumdofs ) = Uk( neumdofs ) + deltaured ;
  currDeltau     = currDeltau     + deltaured ;

function vec = antiSkew( mat )
  vec = [ mat(3,2); mat(1,3); mat(2,1) ] ;

function args = argsAL(analysisSettings, len, neumDofs, timeIndex)
  arcLengthNorm = zeros( len ) ;
  arcLengthNorm(1:2:end) = 1 ;
  arcLengthNorm = arcLengthNorm(neumDofs) ;
  if length( analysisSettings.incremArcLen ) > 1
    incremArcLen = analysisSettings.incremArcLen(timeIndex) ;
  else	
    incremArcLen = analysisSettings.incremArcLen ;
  end
  args = {arcLengthNorm; incremArcLen} ;
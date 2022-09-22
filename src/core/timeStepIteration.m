% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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

function  modelNextSol = timeStepIteration( modelCurrSol, modelProperties, BCsData ) ;

% assign current time (t) variables
% ---------------------------------
Ut         = modelCurrSol.U ; Udott = modelCurrSol.Udot ; Udotdott = modelCurrSol.Udotdot ;
KTtred     = modelCurrSol.systemDeltauMatrix ;
convDeltau = modelCurrSol.convDeltau ;
currLoadFactorsVals = modelCurrSol.currLoadFactorsVals ;

% update time and set candidate displacements and derivatives
% -----------------------------------------------------------
if isempty( modelProperties.analysisSettings.Utp10 )
  Utp1k       = Ut       ;
else
  error('add case for several times')
end

[ Udottp1k, Udotdottp1k, nextTime ] = updateTime( ...
  Ut, Udott, Udotdott, Utp1k, modelProperties.analysisSettings, modelCurrSol.currTime ) ;

% current tangent matrix
% ----------------------
systemDeltauMatrix = KTtred ;

% compute RHS for initial guess Utp1 and in next time step
% --------------------------------------------------------
if strcmp( modelProperties.analysisSettings.methodName, 'arcLength')==1
  nextLoadFactorsVals = currLoadFactorsVals ;
else
  nextLoadFactorsVals = [] ;
end

previous_state_mat = modelCurrSol.previous_state_mat ;

systemDeltauRHS    = modelCurrSol.systemDeltauRHS    ;
systemDeltauMatrix = modelCurrSol.systemDeltauMatrix ;

% --- assemble system of equations ---
[ systemDeltauMatrix, systemDeltauRHS, FextG, ~, nextLoadFactorsVals ] = system_assembler( modelProperties, BCsData, Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k, nextTime, nextLoadFactorsVals, previous_state_mat ) ;

booleanConverged = false                              ;
dispIters        = 0                              ;
currDeltau       = zeros( length( BCsData.neumDofs ), 1 ) ;

%global timeInd
%	timeInd = modelCurrSol.timeIndex ;

while  booleanConverged == 0

  %fprintf(' ============== new iteration ====================\n')
  dispIters = dispIters + 1 ;

  % solve system
  [ deltaured, nextLoadFactorsVals ] = computeDeltaU( systemDeltauMatrix, systemDeltauRHS, dispIters, convDeltau, modelProperties.analysisSettings, nextLoadFactorsVals , currDeltau, modelCurrSol.timeIndex, BCsData.neumDofs ) ;

  % updates: model variables and computes internal forces ---
  [Utp1k, currDeltau] = updateUiter(Utp1k, deltaured, BCsData.neumDofs, currDeltau ) ;

  % --- update next time magnitudes ---
  [ Udottp1k, Udotdottp1k, nextTime ] = updateTime( ...
    Ut, Udott, Udotdott, Utp1k, modelProperties.analysisSettings, modelCurrSol.currTime ) ;

  % --- assemble system of equations ---
  [ systemDeltauMatrix, systemDeltauRHS, FextG, ~, nextLoadFactorsVals ] = system_assembler( modelProperties, BCsData, Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k, nextTime, nextLoadFactorsVals, previous_state_mat ) ;

  % --- check convergence ---
  [ booleanConverged, stopCritPar, deltaErrLoad ] = convergenceTest( modelProperties.analysisSettings, [], FextG(BCsData.neumDofs), deltaured, Utp1k(BCsData.neumDofs), dispIters, [], systemDeltauRHS ) ;
  % ---------------------------------------------------

  % --- prints iteration info in file ---
  printSolverOutput( modelProperties.outputDir, modelProperties.problemName, [ 1 norm(nextLoadFactorsVals) dispIters deltaErrLoad norm(deltaured) ] ) ;

end % iteration while
% --------------------------------------------------------------------

Utp1       = Utp1k ;
Udottp1    = Udottp1k ;
Udotdottp1 = Udotdottp1k ;

% computes KTred at converged Uk
KTtp1red = systemDeltauMatrix ;

% compute stress at converged state
[~, Stresstp1, ~, matFint, strain_vec, acum_plas_strain_vec ] = assembler ( modelProperties.Conec, modelProperties.elements, modelProperties.Nodes, modelProperties.materials, BCsData.KS, Utp1, Udottp1, Udotdottp1, modelProperties.analysisSettings, [ 0 1 0 1 ], modelProperties.nodalDispDamping, nextTime, previous_state_mat ) ;

printSolverOutput( modelProperties.outputDir, modelProperties.problemName, [ 2 (modelCurrSol.timeIndex)+1 nextTime dispIters stopCritPar ] ) ;


% --- (temporary) computation and storage of separated assembled matrices ---
%~ mats  = assembler(  Conec, crossSecsParamsMat, coordsElemsMat, materialsParamsMat, KS, Utp1,   2, Udott, Udotdott, nodalDispDamping, solutionMethod, elementsParamsMat ) ;
%~ ktout = mats{1};

%~ if isunix
  %~ save  'Ktp1.dat' ktout ;
  %~ status = system('tail -n +7 Ktp1.dat > aux.dat' );
  %~ status = system(['mv aux.dat Ktp1_' sprintf('%04i', timeIndex) '.dat'] ) ;
%~ end

%~ if solutionMethod > 2
  %~ dampingMat = mats{2} ;
  %~ massMat    = mats{3} ;

  %~ if isunix
    %~ save  'dampingMattp1.dat' dampingMat ;
    %~ status = system('tail -n +7 dampingMattp1.dat > aux.dat' );
    %~ status = system( ['mv aux.dat dampingMattp1_' sprintf('%04i', timeIndex) '.dat'] ) ;

    %~ save  'massMattp1.dat' massMat ;
    %~ status = system('tail -n +7 massMattp1.dat > aux.dat' );
    %~ status = system( [ 'mv aux.dat massMattp1_' sprintf('%04i', timeIndex) '.dat' ] ) ;
  %~ end

%~ end
% --------------------------------------------------------------------


% %%%%%%%%%%%%%%%%
%~ stabilityAnalysisFlag = stabilityAnalysisBoolean ;
stabilityAnalysisFlag = 0 ;
% %%%%%%%%%%%%%%%%

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

previous_state_mat = [ Stress(:,1) strain_vec acum_plas_strain_vec ] ;

modelNextSol = construct_modelSol( timeIndex, currTime, U , Udot, ...
                                   Udotdot, Stress, convDeltau, ...
                                   nextLoadFactorsVals, systemDeltauMatrix, ...
                                   systemDeltauRHS, timeStepStopCrit, timeStepIters, matFint, previous_state_mat ) ;



% ==============================================================================
%
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

%test
% ==============================================================================
%
% ==============================================================================
function [Uk, currDeltau] = updateUiter(Uk, deltaured, neumdofs, currDeltau )

  oddNeumDofsInds  = find( mod ( neumdofs , 2)==1 ) ;
  evenNeumDofsInds = find( mod ( neumdofs , 2)==0 ) ;

  Uk( neumdofs(oddNeumDofsInds ) ) = Uk( neumdofs(oddNeumDofsInds ) ) + deltaured(oddNeumDofsInds ) ;

  nNodes = length( Uk) / 6 ;

  deltauComplete = zeros( size( Uk)) ;
  deltauComplete( neumdofs ) = deltaured ;

  for i=1:nNodes
    nodeDofs = nodes2dofs( i , 6 ) ;
    nodeAngDofs = nodeDofs(2:2:6)  ;

    %~ updateA = antiSkew( logm( expm( skew( deltauComplete ( nodeAngDofs ) ) ) * ...
                                         %~ expm( skew( Uk             ( nodeAngDofs ) ) ) ...
                                       %~ ) ) ;
    updateB = deltauComplete ( nodeAngDofs ) + Uk             ( nodeAngDofs ) ;
    %~ updateC = logar( expon( deltauComplete ( nodeAngDofs ) ) * ...
                                %~ expon( Uk             ( nodeAngDofs ) ) ) ;

    Uk ( nodeAngDofs ) = updateB ;

  end

  currDeltau      = currDeltau    + deltaured ;

function vec = antiSkew( mat )
  vec = [ mat(3,2) mat(1,3) mat(2,1) ]' ;

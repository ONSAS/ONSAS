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

%% This functions performs the iteration for the computations of the state in
% the next time step using the numerical method and parameters provided by the
% user.

function  [ modelCurrSol, BCsData ] = timeStepIteration( modelCurrSol, BCsData, modelProperties, Utp10 ) ;

startPreviousVels = 0 ; % (0 or 1) - 0 strongly recommended

% ----   extracts variables  ----
modelExtract

if cppSolverBoolean
  cppInterface

else
  
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
  
  % --- start iteration with previous displacements ---
  if isempty( Utp10 )
    Utp1k       = Ut       ;   % initial guess
  else
    Utp1k       = Utp10    ;
  end
  
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
    Conec, crossSecsParamsMat, coordsElemsMat, ...
    materialsParamsMat, KS, constantFext, variableFext, userLoadsFilename, ...
    currLoadFactor, nextLoadFactor, numericalMethodParams, neumdofs, nodalDispDamping, ...
    Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k, elementsParamsMat ) ;
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
    systemDeltauMatrix          = computeMatrix( Conec, crossSecsParamsMat, coordsElemsMat, ...
      materialsParamsMat, KS, Utp1k, neumdofs, numericalMethodParams, ...
      nodalDispDamping, Udott, Udotdott, elementsParamsMat ) ;
    % ---------------------------------------------------
  
    % --- new rhs ---
    [ systemDeltauRHS, FextG ]  = computeRHS( Conec, crossSecsParamsMat, coordsElemsMat, ...
      materialsParamsMat, KS, constantFext, variableFext, ...
      userLoadsFilename, currLoadFactor, nextLoadFactor, numericalMethodParams, ...
      neumdofs, nodalDispDamping, ...
      Ut, Udott, Udotdott, Utp1k, Udottp1k, Udotdottp1k, elementsParamsMat ) ;
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
  Stresstp1 = assembler ( Conec, crossSecsParamsMat, coordsElemsMat, materialsParamsMat, KS, Utp1, 3, Udottp1, Udotdottp1, nodalDispDamping, solutionMethod, elementsParamsMat ) ;



  % --- (temporary) computation and storage of separated assembled matrices ---
  mats  = assembler(  Conec, crossSecsParamsMat, coordsElemsMat, materialsParamsMat, KS, Utp1,   2, Udott, Udotdott, nodalDispDamping, solutionMethod, elementsParamsMat ) ;
  ktout = mats{1}; 
  
  if isunix
    save  'Ktp1.dat' ktout ;
    status = system('tail -n +7 Ktp1.dat > aux.dat' );
    status = system(['mv aux.dat Ktp1_' sprintf('%04i', timeIndex) '.dat'] ) ; 
  end
  
  if solutionMethod > 2
    dampingMat = mats{2} ;
    massMat    = mats{3} ;

    if isunix
      save  'dampingMattp1.dat' dampingMat ;
      status = system('tail -n +7 dampingMattp1.dat > aux.dat' );
      status = system( ['mv aux.dat dampingMattp1_' sprintf('%04i', timeIndex) '.dat'] ) ; 
  
      save  'massMattp1.dat' massMat ;
      status = system('tail -n +7 massMattp1.dat > aux.dat' );
      status = system( [ 'mv aux.dat massMattp1_' sprintf('%04i', timeIndex) '.dat' ] ) ; 
    end

  end
  % --------------------------------------------------------------------  

  
  % %%%%%%%%%%%%%%%%
  stabilityAnalysisFlag = stabilityAnalysisBoolean ;
  % %%%%%%%%%%%%%%%%
  
  if stabilityAnalysisFlag == 2
    [ nKeigpos, nKeigneg, factorCrit ] = stabilityAnalysis ( KTtred, KTtp1red, currLoadFactor, nextLoadFactor ) ;
  elseif stabilityAnalysisFlag == 1
    [ nKeigpos, nKeigneg ] = stabilityAnalysis ( KTtred, KTtp1red, currLoadFactor, nextLoadFactor ) ;
    factorCrit = 0;
  else
    nKeigpos = 0;  nKeigneg = 0; factorCrit = 0 ;
  end
  
  % prints iteration info in file
  printSolverOutput( ...
    outputDir, problemName, timeIndex+1, [ 2 nextLoadFactor dispIters stopCritPar nKeigpos nKeigneg ] ) ;
  
end % if solver C++


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
function [ Udottp1, Udotdottp1, nextTime ] = updateTime(Ut, Udott, Udotdott, Uk, numericalMethodParams, currTime )

  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
  stopTolIts,     targetLoadFactr, nLoadSteps,    ...
  incremArcLen, deltaT, deltaNW, AlphaNW, alphaHHT, finalTime ] ...
      = extractMethodParams( numericalMethodParams ) ;

  nextTime   = currTime + deltaT                          ;

if solutionMethod == 3 || solutionMethod == 4

  if solutionMethod == 4
    deltaNW = (1-2*alphaHHT)/2 ;
    AlphaNW = (1-alphaHHT)^2/4 ;
  end
  
  Udotdottp1 = 1.0/( AlphaNW * (deltaT)^2 ) * ( Uk - Ut ) - 1.0/( AlphaNW * deltaT ) * Udott - ( 1.0/ ( AlphaNW * 2 ) - 1 ) * Udotdott ;

  Udottp1    = Udott + ( ( 1-deltaNW ) * Udotdott + deltaNW * Udotdottp1 ) * deltaT    ;
  
else
  Udotdottp1 = Udotdott ;
  Udottp1    = Udott ;
end


% ==============================================================================
%
% ==============================================================================
function [Uk, currDeltau] = updateUiter(Uk, deltaured, neumdofs, solutionMethod, currDeltau ) 

  oddNeumDofsInds  = find( mod ( neumdofs , 2)==1 ) ;
  evenNeumDofsInds = find( mod ( neumdofs , 2)==0 ) ;

  Uk ( neumdofs(oddNeumDofsInds ) ) = Uk( neumdofs(oddNeumDofsInds ) ) + deltaured (oddNeumDofsInds ) ;

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







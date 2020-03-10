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
% ------------------------------------------------------------------------------

% function for iterative resolution of nonlinear equations.

function ...
%  outputs ---
[ nextLoadFactor, dispIter, stopCritPar, factor_crit, nKeigpos, nKeigneg, Uk, FintGk, Stressk, Strainsk, systemDeltauMatrix ] ...
  = iterativeMethods( ...
% inputs ---
  % constant data
  Conec, secGeomProps, coordsElemsMat, neumdofs, nnodes, hyperElasParamsMat, ...
  numericalMethodParams, constantFext, variableFext, KS, userLoadsFilename , bendStiff, ...
  % model variable data
  Uk, Stressk, Strainsk, FintGk, currLoadFactor, nextLoadFactor, ...
  % specific iterative methods variables 
  convDeltau, stabilityAnalysisBoolean, booleanScreenOutput ) ;
% ------------------------------------------------------------------------------

  % --------------------------------------------------------------------
  % -----------      pre-iteration definitions     ---------------------
  nelems    = size(Conec,1) ; ndofpnode = 6;
  
  booleanConverged = 0 ;
  dispIter         = 0 ;

  % parameters for the Arc-Length iterations
  currDeltau      = zeros( length(neumdofs), 1 ) ;
  
  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
    stopTolIts,     targetLoadFactr, nLoadSteps,    ...
    incremArcLen, deltaT, deltaNW, AlphaNW, finalTime ] ...
        = extractMethodParams( numericalMethodParams ) ;

  % current stiffness matrix for buckling analysis
  if stabilityAnalysisBoolean == 1
    [~, KTtm1 ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 2 ) ;
  end
  % --------------------------------------------------------------------
  % --------------------------------------------------------------------
  

  % --------------------------------------------------------------------
  % --- iteration in displacements (NR) or load-displacements (NR-AL) --
  if  booleanScreenOutput
    fprintf(' iter  normResLoad\n----------------------\n' ) ;
  end
  while  booleanConverged == 0
    dispIter += 1 ;

  auxT = time();
    % system matrix
    systemDeltauMatrix          = computeMatrix( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, neumdofs, solutionMethod, bendStiff);
tiempoComputeMatrix = time() - auxT ;

auxT = time();    
    % system rhs
    [ systemDeltauRHS, FextG ]  = computeRHS( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, dispIter, constantFext, variableFext, userLoadsFilename, currLoadFactor, nextLoadFactor, solutionMethod, neumdofs, FintGk)  ;
tiempoComputeRHS = time() - auxT;

opa = time();
    % computes deltaU
    [deltaured, currLoadFactor] = computeDeltaU ( systemDeltauMatrix, systemDeltauRHS, dispIter, convDeltau(neumdofs), numericalMethodParams, currLoadFactor , currDeltau );
tiempoSystemSolve = time() - opa;
    
    % updates: model variables and computes internal forces
    Uk ( neumdofs ) = Uk(neumdofs ) + deltaured ;
    if solutionMethod == 2
      currDeltau      = currDeltau    + deltaured ;
    end
    [FintGk, ~ ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 1 ) ;

    % check convergence
    [booleanConverged, stopCritPar, deltaErrLoad ] = convergenceTest( numericalMethodParams, FintGk(neumdofs), FextG(neumdofs), deltaured, Uk(neumdofs), dispIter ) ;

    if  booleanScreenOutput
      fprintf(' %3i %12.3e \n' , dispIter, deltaErrLoad ) ;
    end
  end
  if  booleanScreenOutput
    fprintf('----  iters ended --------\n' ) ;
  end

  % --------------------------------------------------------------------
  % --------------------------------------------------------------------



  % computes KTred at converged Uk
  [~, KTt ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 2 ) ;

  if solutionMethod == 2;    
    nextLoadFactor = currLoadFactor ;
  end

  factor_crit = 0;

  [FintGk, ~, Strainsk, Stressk ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 1 ) ;

  if stabilityAnalysisBoolean == 1
    [ factor_crit, nKeigpos, nKeigneg ] = stabilityAnalysis ( KTtm1( neumdofs, neumdofs ), KTt( neumdofs, neumdofs ), currLoadFactor, nextLoadFactor ) ;
  else
    factor_crit = 0; nKeigpos=0; nKeigneg=0;
  end

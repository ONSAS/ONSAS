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


% function for iteration of Newton-Raphson or Newton-Raphson-Arc-Length.

function ...
%  outputs ---
[ nextLoadFactor, dispIter, stopCritPar, factor_crit, nKeigpos, nKeigneg, Uk, FintGk, Stressk, Strainsk, systemDeltauMatrix ] ...
  = analysisNRAndNRAL( ...
% inputs ---
  % constant data
  Conec, secGeomProps, coordsElemsMat, neumdofs, nnodes, hyperElasParamsMat, ...
  numericalMethodParams, constantFext, variableFext, KS, userLoadsFilename , bendStiff, ...
  % model variable data
  Uk, Stressk, Strainsk, FintGk, currLoadFactor, nextLoadFactor, ...
  % specific iterative methods variables 
  convDeltau, stabilityAnalysisBoolean ) ;
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
  [~, KTtm1 ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 2 ) ;
  % --------------------------------------------------------------------
  % --------------------------------------------------------------------
  

  % --------------------------------------------------------------------
  % --- iteration in displacements (NR) or load-displacements (NR-AL) --
  while  booleanConverged == 0
    dispIter += 1 ;

    % system matrix
    systemDeltauMatrix = computeMatrix( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, neumdofs, solutionMethod, bendStiff);

%~ dispIter
%~ full(systemDeltauMatrix)
    
    % system rhs
    [systemDeltauRHS, FextG]    = computeRHS( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, dispIter, constantFext, variableFext, userLoadsFilename, currLoadFactor, nextLoadFactor, solutionMethod, neumdofs, FintGk)  ;

%~ systemDeltauRHS

    % computes deltaU
    [deltaured, currLoadFactor] = computeDeltaU ( systemDeltauMatrix, systemDeltauRHS, dispIter, convDeltau(neumdofs), numericalMethodParams, currLoadFactor , currDeltau );
    
    % updates: model variables and computes internal forces
    Uk ( neumdofs ) = Uk(neumdofs ) + deltaured ;
    if solutionMethod == 2
      currDeltau      = currDeltau    + deltaured ;
    end
    [FintGk, ~ ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 1 ) ;

    % check convergence
    [booleanConverged,stopCritPar] = convergenceTest( numericalMethodParams, FintGk(neumdofs), FextG(neumdofs), deltaured, Uk(neumdofs), dispIter ) ;
 
  end

  % computes KTred at converged Uk
  [~, KTt ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 2 ) ;


  if solutionMethod == 2;    
    nextLoadFactor = currLoadFactor ;
  end


  factor_crit = 0;

  [FintGk, ~, Strainsk, Stressk ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 1 ) ;

  if stabilityAnalysisBoolean == 1
    [ factor_crit, nKeigpos, nKeigneg ] = stabilityAnalysis ( KTtm1( neumdofs, neumdofs ), KTt( neumdofs, neumdofs ), currLoadFactor, nextLoadFactor ) ;
  end
  % -----------------------------------





% ======================================================================
% ======================================================================





% ======================================================================
function systemDeltauMatrix = computeMatrix( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, neumdofs, solutionMethod , bendStiff)

  % computes static tangent matrix
  [~, KT ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, bendStiff, 2 ) ;

  % performs one iteration
  if solutionMethod == 1 || solutionMethod == 2
    systemDeltauMatrix = KT ( neumdofs, neumdofs ) ;
  end
    
    

% ======================================================================

function [systemDeltauRHS, FextG] = computeRHS( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, dispIter, constantFext, variableFext, userLoadsFilename, currLoadFactor, nextLoadFactor, solutionMethod, neumdofs, FintGk) 

  if strcmp( userLoadsFilename , '')
    FextUser = zeros(size(constantFext)) ;
  else
    if solutionMethod == 2
      error('user load not valid for this implementation of arc-length');
    end
    FextUser = feval( userLoadsFilename, nextLoadFactor)  ;
  end

  if solutionMethod == 1

    %~ if (dispIter==1),
      FextG  = variableFext * nextLoadFactor + constantFext  + FextUser ;
    %~ end

    Resred          = FintGk(neumdofs) - FextG(neumdofs) ;
    systemDeltauRHS = - ( Resred ) ;

  elseif solutionMethod == 2

    FextG  = variableFext * currLoadFactor + constantFext ;

    Resred = FintGk(neumdofs) - FextG(neumdofs)  ;

    % incremental displacement
    systemDeltauRHS = [ -Resred  variableFext(neumdofs) ] ;

  end
    


% ======================================================================

function [deltaured, currLoadFactor] = computeDeltaU ( systemDeltauMatrix, systemDeltauRHS, dispIter, redConvDeltau, numericalMethodParams, currLoadFactor, currDeltau  )

  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
    stopTolIts,     targetLoadFactr, nLoadSteps,    ...
    incremArcLen, deltaT, deltaNW, AlphaNW, finalTime ] ...
        = extractMethodParams( numericalMethodParams ) ;
  
  convDeltau = redConvDeltau ;  

  if solutionMethod == 1
    % incremental displacement
    deltaured = systemDeltauMatrix \ systemDeltauRHS ;
  
  elseif solutionMethod == 2
  
    aux = systemDeltauMatrix \ systemDeltauRHS ;
    
    deltauast = aux(:,1) ;  deltaubar = aux(:,2) ;
    
    if dispIter == 1
      if norm(convDeltau)==0
        deltalambda = targetLoadFactr / nLoadSteps ;
      else
        aux = sign( convDeltau' * deltaubar ) ;
        deltalambda =   incremArcLen * aux / ( sqrt( deltaubar' * deltaubar ) ) ;
      end
    else
      ca =    deltaubar' * deltaubar ;
      cb = 2*(currDeltau + deltauast)' * deltaubar ;
      cc = (currDeltau + deltauast)' * (currDeltau + deltauast) - incremArcLen^2 ; 
      disc = cb^2 - 4 * ca * cc ;
      if disc < 0
        disc, error( 'negative discriminant'); 
      end
      sols = -cb/(2*ca) + sqrt(disc) / (2*ca)*[-1 +1]' ;
      
      vals = [ ( currDeltau + deltauast + deltaubar * sols(1) )' * currDeltau;
               ( currDeltau + deltauast + deltaubar * sols(2) )' * currDeltau ] ;
     
      deltalambda = sols( find( vals == max(vals) ) ) ;
    end
    
    currLoadFactor = currLoadFactor + deltalambda(1) ;
    
    deltaured = deltauast + deltalambda(1) * deltaubar ;
  end



% ======================================================================

% --- nonlinear buckling analysis as in section 6.8.2 from Bathe, FEM Procedures 2nd edition. ---

function [ factor_crit, nKeigpos, nKeigneg] = stabilityAnalysis ( KTtm1red, KTtred, currLoadFactor, nextLoadFactor );  

  [a,b] = eig( KTtred ) ;
  Keigvals = diag(b) ; 
  nKeigpos = length( find(Keigvals >  0 ) ) ;
  nKeigneg = length( find(Keigvals <= 0 ) ) ;

  [vecgamma, gammas ] = eig( KTtred, KTtm1red ) ;
  
  gammas = diag( gammas);
 
  if length( find( gammas >  0 ) ) > 0,
  
    gamma_crit  = min ( gammas ( find( gammas >  0 ) ) ) ;
    if gamma_crit ~= 1 
      lambda_crit = 1 / ( 1 - gamma_crit )  ;               
      factor_crit = currLoadFactor + lambda_crit * (nextLoadFactor - currLoadFactor) ;
    else
      factor_crit = 0 ;
    end
  else
    factor_crit = 0;
  end




% ======================================================================
% ======================================================================
function [ booleanConverged, stopCritPar ] = convergenceTest( ...
  numericalMethodParams, redFint, redFext, redDeltaU, redUk, dispIter ) 

  [ solutionMethod, stopTolDeltau,   stopTolForces, ...
    stopTolIts,     targetLoadFactr, nLoadSteps,    ...
    incremArcLen, deltaT, deltaNW, AlphaNW, finalTime ] ...
        = extractMethodParams( numericalMethodParams ) ;

  normaUk       = norm( redUk )               ;
  normadeltau   = norm( redDeltaU         )   ;
  deltaErrLoad  = norm( redFint - redFext )   ;
  normFext      = norm( redFext )             ;
  
  logicDispStop = ( normadeltau  < ( normaUk  * stopTolDeltau ) )  ;
  logicForcStop = ( deltaErrLoad < ( normFext * stopTolForces ) )  ;
                
  if logicForcStop
    stopCritPar = 1 ;      booleanConverged = 1 ;

  elseif logicDispStop
    stopCritPar = 2 ;      booleanConverged = 1 ;

  elseif ( dispIter >= stopTolIts )
    warning('displacements iteration stopped by max iterations.');
    stopCritPar = 3 ;      booleanConverged = 1 ;
  else
    booleanConverged = 0;  stopCritPar = [];
  end
  


% ======================================================================
% ======================================================================
function [ solutionMethod, stopTolDeltau,   stopTolForces, ...
           stopTolIts,     targetLoadFactr, nLoadSteps,    ...
           incremArcLen, deltaT, deltaNW, AlphaNW, finalTime ] ...
           = extractMethodParams( numericalMethodParams ) 

  solutionMethod   = numericalMethodParams(1) ;
  
  if ( solutionMethod == 1) || ( solutionMethod == 2)

    % ----- resolution method params -----
    stopTolDeltau    = numericalMethodParams(2) ;
    stopTolForces    = numericalMethodParams(3) ;
    stopTolIts       = numericalMethodParams(4) ;
    targetLoadFactr  = numericalMethodParams(5) ;
    nLoadSteps       = numericalMethodParams(6) ;
  
    if solutionMethod ==2
      incremArcLen     = numericalMethodParams(7) ;
    else
      incremArcLen = [] ;
    end
    
    deltaT = []; finalTime = []; deltaNW = []; AlphaNW = [] ;
  
  else
    deltaT         = numericalMethodParams(2)        ;
    finalTime      = numericalMethodParams(3)        ;
    stopTolDeltau  = numericalMethodParams(4)        ;
    stopTolForces  = numericalMethodParams(5)        ;
    stopTolIts     = numericalMethodParams(6)        ;
    deltaNW        = numericalMethodParams(7)        ;
    AlphaNW        = numericalMethodParams(8)        ;

    targetLoadFactr = [] ; nLoadSteps = []; incremArcLen = [] ;    
    
  end
  
% ======================================================================

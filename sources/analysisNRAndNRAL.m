%function for iteration of Newton-Raphson or Newton-Raphson-Arc-Length.

%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

% ------------------------------------------------------------------------------
function ...
%  outputs ---
[ nextLoadFactor, dispIter, stopCritPar, factor_crit, nKeigpos, nKeigneg, Uk, FintGk, Stressk, Strainsk, Dsigdeps ] ...
  = analysisNRAndNRAL( ...
% inputs ---
  % constant data
  Conec, secGeomProps, coordsElemsMat, neumdofs, nnodes, hyperElasParamsMat, ...
  numericalMethodParams, constantFext, variableFext, KS, userLoadsFilename , ...
  % model variable data
  Uk, Stressk, Strainsk, dsigdepsk, FintGk, currLoadFactor, nextLoadFactor, ...
  % specific iterative methods variables 
  convDeltau ) ;
% ------------------------------------------------------------------------------

  % ----- resolution method params -----
  solutionMethod   = numericalMethodParams(1) ;
  stopTolDeltau    = numericalMethodParams(2) ;
  stopTolForces    = numericalMethodParams(3) ;
  stopTolIts       = numericalMethodParams(4) ;
  targetLoadFactr  = numericalMethodParams(5) ;
  nLoadSteps       = numericalMethodParams(6) ;

  if solutionMethod ==2
    incremArcLen     = numericalMethodParams(7) ;
  end
  % --------------------------------------

  ndofpnode = 6;           nelems   = size(Conec,1) ;
  
  iterDispConverged = 0 ;  dispIter = 0 ;

  % parameters for the Arc-Length iterations
  currDeltau      = zeros( length(neumdofs), 1 ) ;
  convDeltau      = convDeltau(neumdofs);

  [~, KTtm1 ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, 2 ) ;

  % ----------------------------------
  % iteration in displacements (NR) or load-displacements (NR-AL)
  while ( iterDispConverged == 0 )
    dispIter += 1 ;

    % computes tangent matrix
    [~, KT ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, 2 ) ;

    % performs one iteration
    if solutionMethod == 1

      if (dispIter==1),

        if strcmp( userLoadsFilename , '')
          FextUser = zeros(size(constantFext));
        else
          FextUser = feval( userLoadsFilename, nextLoadFactor)  ;
        end
        FextG  = variableFext * nextLoadFactor + constantFext  + FextUser ;
      end

      NRIter  % performs one newton-raphson iteration

    elseif solutionMethod == 2
      NRALIter     % performs one newton-raphson-arc-length iteration
    end
    % --------------------------

    % updates model variables and computes internal forces
    [FintGk, ~ ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk,1 ) ;
   
    % --- stopping criteria verification ---
    deltaErrLoad  = norm(FintGk(neumdofs) - FextG(neumdofs) )        ;
    normFext      = norm(FextG(neumdofs) )                           ;
    logicDispStop = ( normadeltau  < ( normaUk  * stopTolDeltau ) )  ;
    logicForcStop = ( deltaErrLoad < ( normFext * stopTolForces ) )  ;
                  
    if logicForcStop
      stopCritPar = 1 ;      iterDispConverged = 1 ;
  
    elseif logicDispStop
      stopCritPar = 2 ;      iterDispConverged = 1 ;
  
    elseif ( dispIter >= stopTolIts )
      warning('displacements iteration stopped by max iterations.');
      stopCritPar = 3 ;      iterDispConverged = 1 ;
    end
    % -------------------------
    
  end

  % computes KTred at converged Uk
  [~, KTt ] = assemblyFintVecTangMat( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk, 2 ) ;


  if solutionMethod == 2;    
    nextLoadFactor = currLoadFactor ;
  end


  factor_crit = 0;

  % -----------------------------------
  % buckling analysis

  [FintGk, ~, Strainsk, Stressk, Dsigdeps ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk,1 ) ;


  KTtred   = KTt  ( neumdofs, neumdofs );
  KTtm1red = KTtm1( neumdofs, neumdofs );

  [a,b] = eig( KTtred ) ;
  Keigvals = diag(b) ; 
  nKeigpos = length( find(Keigvals >  0 ) ) ;
  nKeigneg = length( find(Keigvals <= 0 ) ) ;

  [vecgamma, gammas ] = eig( KTtm1red, KTtred ) ;
  
  gammas = diag( gammas);
  
  if length( find( gammas >  0 ) ) > 0
    gamma_crit = min ( gammas ( find( gammas >  0 ) ) ) ;
    lambda_crit  = 1 / ( 1 - gamma_crit ) ;
    factor_crit = lambda_crit * currLoadFactor ;

  else
    factor_crit = 0;
  end
  % -----------------------------------


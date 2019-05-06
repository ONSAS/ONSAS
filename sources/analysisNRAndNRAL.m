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
    logicDispStop = ( normadeltau < ( normaUk * stopTolDeltau ) )    ;
    logicForcStop = (  deltaErrLoad < ( normFext * stopTolForces ) ) ;
                  
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

  if solutionMethod == 2;    
    nextLoadFactor = currLoadFactor ;
  end


  nKeigpos =  0 ;
  nKeigneg =  0 ;
    factor_crit = 0;

  % -----------------------------------
  % buckling analysis
  
  % computes KTred at converged Uk
  %~ [~, KT, KL0 ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk,2 ) ;

  [FintGk, ~, Strainsk, Stressk, Dsigdeps ] = assemblyFintVecTangMat ( Conec, secGeomProps, coordsElemsMat, hyperElasParamsMat, KS, Uk,1 ) ;

  %~ KTred  = KT  ( neumdofs, neumdofs );

  %~ [a,b] = eig( KTred ) ;
  %~ Keigvals = diag(b) ; 
  %~ nKeigpos = length( find(Keigvals >  0 ) );
  %~ nKeigneg = length( find(Keigvals <= 0 ) );

  %~ KL0red = KL0 ( neumdofs, neumdofs );
  %~ [a, lambtech ] = eig( KTred ,  KL0red ) ;    
  %~ lambtech = diag(lambtech) ;
  %~ %
  %~ if length( find( lambtech >  0 ) ) > 0
    %~ lambdatech_crit = min ( lambtech ( find( lambtech >  0 ) ) ) ;
    %~ lambda_crit  = 1 / ( 1 - lambdatech_crit ) ;
    %~ factor_crit = lambda_crit * currLoadFactor ;
  %~ else
    %~ factor_crit = 0;
  %~ end

  % linearized according to Bathe
  %~ if (loadIter == 1)
    %~ KG0redBathe = KGred / currLoadFactor ;
    %~ [a,b] = eig( KL0red , - KG0redBathe ) ;
    %~ factor_crit_lin_Bathe = min( diag(b) );
  %~ end
  % -----------------------------------


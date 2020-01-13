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



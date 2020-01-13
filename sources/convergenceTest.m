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
  

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

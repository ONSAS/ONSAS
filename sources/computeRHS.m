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
    

function [cl, cd, cm] = uBEMinterpAeroParams(aeroCoefs, elemIDsection, betaRel)

    aoast    = aeroCoefs{1} ;
    liftCoef = aeroCoefs{2} ;
    dragCoef = aeroCoefs{3} ;
    momCoef  = aeroCoefs{4} ;

    id = elemIDsection  ;

    [m, n] = size(aoast) ;

    % Compute Cl, Cd and Cm of elem section
    if n == 1
        cl    =  interp1( aoast(:,1), liftCoef(:,1),  rad2deg(betaRel), 'linear', 'extrap' );
        cd    =  interp1( aoast(:,1), dragCoef(:,1),  rad2deg(betaRel), 'linear', 'extrap' );
        cm    =  interp1( aoast(:,1),  momCoef(:,1),  rad2deg(betaRel), 'linear', 'extrap' );
    else
        cl    =  interp1( aoast(:,id), liftCoef(:,id),  rad2deg(betaRel), 'linear', 'extrap' );
        cd    =  interp1( aoast(:,id), dragCoef(:,id),  rad2deg(betaRel), 'linear', 'extrap' );
        cm    =  interp1( aoast(:,id),  momCoef(:,id),  rad2deg(betaRel), 'linear', 'extrap' );
    end
end
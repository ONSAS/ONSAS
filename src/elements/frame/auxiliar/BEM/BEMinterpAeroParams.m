function [cl, cd, cm] = BEMinterpAeroParams(aeroCoefs, elemIDsection, betaRel)

    aoast    = aeroCoefs(:,1) ;
    liftCoef = aeroCoefs(:,2) ;
    dragCoef = aeroCoefs(:,3) ;
    momCoef  = aeroCoefs(:,4) ;

    id = elemIDsection  ;

    % Compute Cl, Cd and Cm of elem section
    cl    =  interp1( aoast(:,id), liftCoef(:,id),  rad2deg(betaRel), 'linear', 'extrap' );
    cd    =  interp1( aoast(:,id), dragCoef(:,id),  rad2deg(betaRel), 'linear', 'extrap' );
    cm    =  interp1( aoast(:,id), momCoef(:,id) ,  rad2deg(betaRel), 'linear', 'extrap' );
end
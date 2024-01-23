function [inducedVeln, inducedIntVeln, inducedQSVeln] = uBEMinducedVelocity(lift, VpiRel, inducedVeln1, ...
                            inducedQSVeln1, inducedIntVeln1, radio, Rrot, Rhub, betaRel, rho, ...
                            deltaT, L2, Rroof, Rr, DWMbool) 

% Project induced wake componentes into deformed coordiantes
inducedVeln1    = Rroof'*Rr'*inducedVeln1 ;
inducedQSVeln1  = Rroof'*Rr'*inducedQSVeln1 ;
inducedIntVeln1 = Rroof'*Rr'*inducedIntVeln1 ;

V0 = norm( VpiRel ) ;

a  = dot( inducedVeln1', [0, -1, 0] )/V0;
if(a<0.333)
    fglau=1;
else
    fglau=0.25*(5-3*a);
end

% Apply prandtl tip and hub correction factor
F     = uBEMprandtlFactor(radio, Rrot, Rhub, rad2deg(betaRel));
aux1  = norm( VpiRel + fglau*[0, -1, 0]'*dot( inducedVeln1',[0, -1, 0] ));
aux2  = aux1*4*pi*rho*radio*F;
        
inducedQSVeln = [0; -3*norm(lift)*sin(rad2deg(betaRel))/aux2; -3*norm(lift)*cos(rad2deg(betaRel))/aux2 ] ;

if DWMbool
    [inducedVeln, inducedIntVeln] = uBEMdynamicInflow(inducedVeln1, inducedQSVeln, ...
                            inducedQSVeln1, inducedIntVeln1, V0, radio, Rrot, a, deltaT) ;
    
    inducedVeln    = ( (Rroof'*Rr')'*inducedVeln  )'   ; % Convert back to global coordinates
    inducedIntVeln = ( (Rroof'*Rr')'*inducedIntVeln )'; % Convert back to global coordinates
    inducedQSVeln  = ( (Rroof'*Rr')'*inducedQSVeln )' ; % Convert back to global coordinates
else
    inducedVeln    = ( (Rroof'*Rr')'*inducedQSVeln )';
    inducedIntVeln = [0, 0, 0] ;
end


function [inducedVeln, inducedIntVeln, inducedQSVeln] = uBEMinducedVelocity(lift, vFlowG, VpiRel, inducedVeln1, ...
                            inducedQSVeln1, inducedIntVeln1, radio, Rrot, Rhub, betaRel, rho, ...
                            deltaT, Rroof, Rr, R0, Rrblade, Rcone, L2, DWMbool) 

% Project induced wake componentes into deformed coordiantes
inducedVeln1    = L2*Rroof'*Rr'*inducedVeln1 ;
inducedQSVeln1  = L2*Rroof'*Rr'*inducedQSVeln1 ;
inducedIntVeln1 = L2*Rroof'*Rr'*inducedIntVeln1 ;

% Normal vector to rotor plane projected into deformed coordiantes
n  = (Rroof'*Rr')*Rrblade*Rcone'*[ 0, 0, -1 ]';

% Project wind velocity into deformed coordinates
V0     = L2*Rroof'*Rr'*vFlowG;
Vprime = norm( V0 + n*dot( inducedVeln1, n ));

a  = ( norm(V0) - Vprime ) / norm(V0);
if(a<0.333)
    fglau=1;
else
    fglau=0.25*(5-3*a);
end

% Apply prandtl tip and hub correction factor
F     = uBEMprandtlFactor(radio, Rrot, Rhub, betaRel);
aux1  = norm( V0 + fglau*n*dot( inducedVeln1, n ));
aux2  = aux1*4*pi*rho*radio*F;
        
inducedQSVeln = -3*lift/aux2 ;

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
end

function [ inducedVel, inducedIntVel, inducedQSVel, nodIdx1, nodIdx2 ] = uBEMupdateInducedVelocity(analysisSettings, elemCoords, ...
                                                                                inducedVeln1, inducedIntVeln1, inducedQSVeln1,  ...
                                                                                Ue, Udote, uBEMdataCoords, ...
                                                                                nextTime, DWMbool)

% HAWT global parameters
global polars ; aeroData   = polars ; 
global radius ; radiusList = radius ;
global chord  ; chordList  = chord  ;
global twist  ; twistList  = twist  ;
global Rrot   ; global Rhub;

%% compute corotational rotation matrices
% Elem reference coordinates:
xs = elemCoords(:);

%% extract fluid properties
densityFluid   = analysisSettings.fluidProps{1,1} ;
userFlowVel    = analysisSettings.fluidProps{3,1} ;

%% Import blade data including, radius, twist, chord and thick
[ aoast, liftCoef, dragCoef, momCoef ] = uBEMAeroProps( aeroData ) ;

aeroCoefs{1} = aoast    ;
aeroCoefs{2} = liftCoef ;
aeroCoefs{3} = dragCoef ;
aeroCoefs{4} = momCoef  ;

%% Define frame geom properties
[elemTwist, elemChord, elemIDsection, elemRadius, nodIdx1, nodIdx2] = uBEMframeProps(uBEMdataCoords, xs);
% uBEM system HAWT configuration matrix
global Rbem;   global R12;     global R34;
if nodIdx2 <= (3 + (length(radiusList) + 1)*1)
    Rb = Rbem(:,:,1) ;
elseif nodIdx2 > (3 + (length(radiusList) + 1)*1) && nodIdx2 <= (3 + (length(radiusList) + 1)*2)
    Rb = Rbem(:,:,2) ;
elseif nodIdx2 > (3 + (length(radiusList) + 1)*2) && nodIdx2 <= (3 + (length(radiusList) + 1)*3)
    Rb = Rbem(:,:,3) ;
end

%% Compute co-rotational transformation matrix
[R0, Rr, Rg1, Rg2, ~, ~] = corotRotMatrices( Ue, elemCoords ) ;

baseBool = false ;
% local rotations
if ~baseBool ;
    Rroof1 = Rr' * Rg1 * R0 ;% Eq(30) J.-M. Battini 2002
    Rroof2 = Rr' * Rg2 * R0 ;% Eq(30) J.-M. Battini 2002
else
    Rroof1 = Rr' * R0 * Rg1 ;
    Rroof2 = Rr' * R0 * Rg2 ;
end

tl1 = logar( Rroof1 ) ; % Eq(31) J.-M. Battini 2002
tl2 = logar( Rroof2 ) ;% Eq(31) J.-M. Battini 2002

% auxiliary matrix created for uBEM method to compute twist angle
L1 = [ 1 0 0  ;
       0 0 0  ;
       0 0 0] ;

% auxiliary matrix created to project transversal velocity
L2 = [ 0 0 0 ;
       0 1 0   ;
       0 0 1 ] ;

% 90 degrees rotation matrix
L3 = expon( [pi/2 0 0] ) ;

%% -------------------------------------------------------------------------------- 
% Convert twist angle and chord from uBEM local coordinate system to
% global coordiante system
global C0;
C0    = R0'*Rb;
elemTwist         = [ (elemTwist(1:3)')'     (elemTwist(4:6)')'    ] ; % twist angle expresed in global HAWT coordinates
chordVectorG      = [ (R0*elemChord(1:3)')'  (R0*elemChord(4:6)')' ] ; % chord local coords expresed in global HAWT coordinates
dimCharacteristic = [ norm( elemChord(1:3) ), norm( elemChord(4:6) ) ]' ;

% Node 1
% Velocity params
% Node 1 frame velocity
udotFrame1 = Udote( 1:2:6 ) ;
% Node 1 flow velocity
udotFlowNode1 = feval( userFlowVel, elemCoords(1:3)' + Ue(1:2:6), nextTime )  ;
% Node 1 induced velocity of previous time step
inducedVelNod1n1    = inducedVeln1(1:3)';
inducedQSVelNod1n1  = inducedQSVeln1(1:3)';
inducedIntVelNod1n1 = inducedIntVeln1(1:3)';

% Compute relative velocity of node 1 
[VpiRelNode1, VpiRelperpNode1, VrelGnode1] = uBEMcomputeVpiRels( udotFlowNode1, udotFrame1, inducedVelNod1n1, ...
                                                                Rroof1, Rr, L2, L3 ) ;

% Compute angle of attack of node 1
defChordVecNod1 = expon( deg2rad( elemTwist(1:3) ) )*( Rroof1'*Rr'*chordVectorG(1:3)') ;
tch1            = ( defChordVecNod1 / norm( defChordVecNod1 )) ;
if( norm( VpiRelNode1 ) == 0 )
    td1 = tchNod1 ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
else % the drag direction at a generic cross section in deformed coordinates is:
    td1 = VpiRelNode1 / norm( VpiRelNode1 ) ;
end

cosBetanod1 = dot( tch1, td1 ) / ( norm(td1) * norm(tch1) ) ;
sinBetanod1 = dot( cross(tch1, td1), [1 0 0] ) / ( norm( td1 ) * norm( tch1 ) ) ;
betaRelnod1 = sign( sinBetanod1 ) * acos( cosBetanod1 ) ;

% Compute node 1 aero params
global clstat1; global cdstat1; global cmstat1;
[clstat1, cdstat1, cmstat1] = uBEMinterpAeroParams(aeroCoefs, elemIDsection(1), betaRelnod1);

% Compute node 1 induced velocity
fll1  =  1/2 * densityFluid * clstat1 / 2 * dimCharacteristic(1) * norm( VpiRelNode1 ) * VpiRelperpNode1 ;

[inducedVelNod1, inducedIntVelNod1, inducedQSVelNod1] = uBEMinducedVelocity(fll1, udotFlowNode1, VpiRelNode1, inducedVelNod1n1, inducedQSVelNod1n1, ...
                                                     inducedIntVelNod1n1, elemRadius(1), Rrot, Rhub, betaRelnod1, densityFluid, ...
                                                     analysisSettings.deltaT, Rroof1, Rr, R0, Rb, R34, DWMbool) ;

% ----------------------------------------------------------------------------------
% Node 2
% Velocity params

% Node 2 frame velocity
udotFrame2 = Udote( 7:2:12 ) ;
% Node 2 flow velocity
udotFlowNode2 = feval( userFlowVel, elemCoords(4:6)' + Ue(7:2:12), nextTime )  ;
% Node 2 induced velocity of previous time step
inducedVelNod2n1    = inducedVeln1(4:6)';
inducedQSVelNod2n1  = inducedQSVeln1(4:6)';
inducedIntVelNod2n1 = inducedIntVeln1(4:6)';

% Compute relative velocity of node 2 
[VpiRelNode2, VpiRelperpNode2, VrelGnode2] = uBEMcomputeVpiRels( udotFlowNode2, ...
                            udotFrame2, inducedVelNod2n1, Rroof2, Rr, L2, L3 );

% Compute angle of attack of node 2
chordVecNod2 = expon( deg2rad( elemTwist(4:6) ) )*( Rroof2'*Rr'*chordVectorG(4:6)') ;
tch2         = ( chordVecNod2 / norm( chordVecNod2 )) ;
if( norm( VpiRelNode2 ) == 0 )
    td2 = tchNod2 ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
else % the drag direction at a generic cross section in deformed coordinates is:
    td2 = VpiRelNode2 / norm( VpiRelNode2 ) ;
end

cosBetanod2 = dot( tch2, td2 ) / ( norm(td2) * norm(tch2) ) ;
sinBetanod2 = dot( cross(tch2, td2), [1 0 0] ) / ( norm( td2 ) * norm( tch2 ) ) ;
betaRelnod2 = sign( sinBetanod2 ) * acos( cosBetanod2 ) ;
global betaTest; 
betaTest = [betaTest, rad2deg(betaRelnod2)];

% Compute node 1 aero params
[clstat2, cdstat2, cmstat2] = uBEMinterpAeroParams(aeroCoefs, elemIDsection(2), betaRelnod2);
global cl;
cl = [ cl, clstat2 ];

% Compute node 1 induced velocity
fll2  =  1/2 * densityFluid * clstat2 / 2 * dimCharacteristic(2) * norm( VpiRelNode2 ) * VpiRelperpNode2 ;

[inducedVelNod2, inducedIntVelNod2, inducedQSVelNod2] = uBEMinducedVelocity(fll2, udotFlowNode2, VpiRelNode2, inducedVelNod2n1, inducedQSVelNod2n1, ...
                                                     inducedIntVelNod2n1, elemRadius(2), Rrot, Rhub, betaRelnod2, densityFluid, ...
                                                     analysisSettings.deltaT, Rroof1, Rr, R0, Rb, R34, DWMbool) ; 

inducedVel     = [ ( (Rroof1'*Rr')'*inducedVelNod1'  )'    , ( (Rroof2'*Rr')'*inducedVelNod2'    )' ];
inducedIntVel  = [ ( (Rroof1'*Rr')'*inducedIntVelNod1' )'  , ( (Rroof2'*Rr')'*inducedIntVelNod2' )' ]; 
inducedQSVel   = [ ( (Rroof1'*Rr')'*inducedQSVelNod1'  )'  , ( (Rroof2'*Rr')'*inducedQSVelNod2'  )' ];

end

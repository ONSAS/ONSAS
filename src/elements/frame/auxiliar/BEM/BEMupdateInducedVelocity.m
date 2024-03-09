function [ waket1p, wakeINTt1p, wakeQSt1p ] = BEMupdateInducedVelocity( analysisSettings, elemCoords, ...
                                                                        elemWakeOld, elemWakeINTOld, elemWakeQSOld,      ...
                                                                        Ue, Udote, BEMparams, airFoilPolars, ...
                                                                        dynStallParams, nextTime, DWMbool)
global Rb; global Rcone; global Rrot; global Rhub;

% HAWT global parameters
%% compute corotational rotation matrices
% Elem reference coordinates:
xs = elemCoords(:);

%% extract fluid properties
densityFluid   = analysisSettings.fluidProps{1,1} ;
userFlowVel    = analysisSettings.fluidProps{3,1} ;

%% Define frame geom properties
[ nodID, nodRadio, nodChords, nodTwist, nodes ] = BEMframeProps(BEMparams, xs);

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

% Node 1
% Velocity params
% Node 1 frame velocity
udotFrame1 = Udote( 1:2:6 ) ;
% Node 1 flow velocity
udotFlowNode1 = feval( userFlowVel, elemCoords(1:3)' + Ue(1:2:6), nextTime )  ;
% Node 1 induced velocity of previous time step
inducedVelNod1n1 = elemWakeOld(1:2:6 )';

% Compute relative velocity of node 1 
[VpiRelNode1, VpiRelperpNode1, VrelGnode1] = BEMcomputeVpiRels( udotFlowNode1, udotFrame1, inducedVelNod1n1, ...
                                                                Rroof1, Rr, L2, L3 ) ;

% Compute angle of attack of node 1
dimCharacteristic1 = norm( nodChords(:,1) );
tchRef1 = expon(-nodTwist(:,1))*nodChords(:,1)  ; 
tch1    = tchRef1 / norm( tchRef1 ) ;
if( norm( VpiRelNode1 ) == 0 )
    td1 = tchNod1 ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
else % the drag direction at a generic cross section in deformed coordinates is:
    td1 = VpiRelNode1 / norm( VpiRelNode1 ) ;
end

cosBetanod1 = dot( td1, tch1 ) / ( norm(td1) * norm(tch1) ) ;
sinBetanod1 = dot( cross(td1, tch1), [1 0 0] ) / ( norm( td1 ) * norm( tch1 ) ) ;
betaRelnod1 = sign( sinBetanod1 ) * acos( cosBetanod1 );
% global betaTest1; 
% betaTest1 = [betaTest1, rad2deg(betaRelnod1)];

% Compute node 1 aero params
[clstat(1), ~, ~] = BEMinterpAeroParams( airFoilPolars, nodID(1), betaRelnod1);

% Compute node 1 induced velocity
fll1  =  -1/2 * densityFluid * clstat(1) / 2 * dimCharacteristic1 * norm( VpiRelNode1 ) * VpiRelperpNode1 ;
% Static and Inter value of wake model
inducedQSVelNod1n1  = elemWakeQSOld(1:2:6 )'  ;
inducedIntVelNod1n1 = elemWakeINTOld(1:2:6 )' ;

[inducedVelNod1, inducedINTVelNod1, inducedQSVelNod1] = BEMinducedVelocity(fll1, udotFlowNode1, VpiRelNode1, inducedVelNod1n1, inducedQSVelNod1n1, ...
                                                        inducedIntVelNod1n1, nodRadio(1), Rrot, Rhub, betaRelnod1, densityFluid, ...
                                                        analysisSettings.deltaT, Rroof1, Rr, R0, Rb, Rcone, L2, DWMbool) ;

% ----------------------------------------------------------------------------------
% Node 2
% Velocity params

% Node 2 frame velocity
udotFrame2 = Udote( 7:2:12 ) ;
% Node 2 flow velocity
udotFlowNode2 = feval( userFlowVel, elemCoords(4:6)' + Ue(7:2:12), nextTime )  ;
% Node 2 induced velocity of previous time step
inducedVelNod2n1 = elemWakeOld(7:2:12 )';

% Compute relative velocity of node 2 
[VpiRelNode2, VpiRelperpNode2, VrelGnode2] = BEMcomputeVpiRels( udotFlowNode2, ...
                            udotFrame2, inducedVelNod2n1, Rroof2, Rr, L2, L3 );

% Compute angle of attack of node 2
dimCharacteristic2 = norm( nodChords(:,2) );
tchRef2 = expon(-nodTwist(:,2))*nodChords(:,2)  ; 
tch2    = tchRef2 / norm( tchRef2 ) ;
if( norm( VpiRelNode2 ) == 0 )
    td2 = tchNod2 ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
else % the drag direction at a generic cross section in deformed coordinates is:
    td2 = VpiRelNode2 / norm( VpiRelNode2 ) ;
end

cosBetanod2 = dot( td2, tch2 ) / ( norm(td2) * norm(tch2) ) ;
sinBetanod2 = dot( cross(td2, tch2), [1 0 0] ) / ( norm( td2 ) * norm( tch2 ) ) ;
betaRelnod2 = sign( sinBetanod2 ) * acos( cosBetanod2 );     

% Compute node 1 aero params
[clstat(2), ~, ~] = BEMinterpAeroParams( airFoilPolars, nodID(2), betaRelnod2);

% Compute node 1 induced velocity
fll2  =  -1/2 * densityFluid * clstat(2) / 2 * dimCharacteristic2 * norm( VpiRelNode2 ) * VpiRelperpNode2 ;

% Static and Inter value of wake model
inducedQSVelNod2n1  = elemWakeQSOld(7:2:12 )'  ;
inducedIntVelNod2n1 = elemWakeINTOld(7:2:12 )' ;

[inducedVelNod2, inducedINTVelNod2, inducedQSVelNod2] = BEMinducedVelocity(fll2, udotFlowNode2, VpiRelNode2, inducedVelNod2n1, inducedQSVelNod2n1, ...
                                                        inducedIntVelNod2n1,  nodRadio(2), Rrot, Rhub, betaRelnod2, densityFluid, ...
                                                        analysisSettings.deltaT, Rroof2, Rr, R0, Rb, Rcone, L2, DWMbool) ; 

% global check; global param;

waket1p     = [ inducedVelNod1'    , inducedVelNod2'    ];
wakeINTt1p  = [ inducedINTVelNod1' , inducedINTVelNod2' ]; 
wakeQSt1p   = [ inducedQSVelNod1'  , inducedQSVelNod2'  ];

end

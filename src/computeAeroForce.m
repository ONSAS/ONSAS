function fag = computeAreoFroce(  elemCoords, elemCrossSecParams, Ue, Udote, Udotdote,...
                                  elemAeroParams, userWindVel, numGaussPoints)
    % Nodal Winds:
    udotWindNode1 = feval( userWindVel, elemCoords(1) ); 
    udotWindNode2 = feval( userWindVel, elemCoords(4) ); 
    udotWindElem  = [udotWindNode1; udotWindNode2];
    % Elem reference coordinates:
    xs = elemCoords(:);

    % Read aerodinamic profile
    vecChordUndef = elemAeroParams(1:3)';
    dimCaracteristic = elemAeroParams(4); 

    % Material and cross section props:
    [Area, J, Iyy, Izz, ~ ] = crossSectionProps ( elemCrossSecParams, 0 );

    % Change indexes according to battini's nomenclature
    permutIndxs = [1:2:5 2:2:6 ([1:2:5]+6) ([2:2:6]+6) ];
    dg       = Ue      ( permutIndxs );
    ddotg    = Udote   ( permutIndxs );
    ddotdotg = Udotdote( permutIndxs );   
    
    % Compute rotations matrixes:
    % rotation global matrices
    tg1 = dg(  4:6  );
    tg2 = dg( 10:12 );
    Rg1 = expon( tg1 );
    Rg2 = expon( tg2 );

    % rotation matrix to reference configuration
    x21 = xs(4:6) - xs(1:3);
    d21 = dg(7:9) - dg(1:3);
    lo = sqrt( ( x21       )' * ( x21       ) ); %
    l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ); %
    Ro = beamRefConfRotMat( x21 );

    % rigid rotation matrix:
    % deformed x axis
    e1 = ( x21 + d21 ) / l;
    q1 = Rg1 * Ro * [0 1 0]';
    q2 = Rg2 * Ro * [0 1 0]';
    q  = ( q1 + q2 ) / 2; 
    % deformed z local axis
    e3 = cross (e1, q);
    e3 = e3 / norm(e3); % normalization
    % deformed y local axis
    e2 = cross (e3, e1);
    Rr = [ e1 e2 e3 ];
    
    % Compute nus eneries in reference configuration
    q  = Rr' *  q;
    q1 = Rr' * q1;
    nu = q(1)/q(2);
    nu11 = q1(1)/q(2);
    nu12 = q1(2)/q(2);
    nu21 = 2*nu-nu11;
    nu22 = 2-nu12;

    % local rotations
    Re1 = Rr' * Rg1 * Ro;
    Re2 = Rr' * Rg2 * Ro;
    tl1 = logar( Re1 );
    tl2 = logar( Re2 );

    %auxiliar matrix
    I3 = eye(3);
    O3 = zeros(3);
    O1 = zeros(1,3);
    
    II=[O3 I3 O3 O3
        O3 O3 O3 I3];

    G=[ 0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
        0   0    1/l     0        0     0  0  0    -1/l     0        0     0
        0  -1/l  0       0        0     0  0  1/l   0       0        0     0]';    

    P = II - [G'; G'];
    %tensor to rotat magnitudes to rigid configuration
    EE=[Rr O3 O3 O3
        O3 Rr O3 O3
        O3 O3 Rr O3
        O3 O3 O3 Rr];
    
    % auxilair created to proyect transversal velocity
    L2 = [0 0 0 
          0 1 0
          0 0 1];
    L3 = expon([-pi/2 0 0]);         

    %angular velocity from rigid component
    wdoter = G' * EE' * ddotg;% Eq. 65

    %Extract points and wieghts for numGausspoints selected
    [xIntPoints, wIntPoints] = GaussPointsAndWeights( numGaussPoints );


    % Compute fag local by using RR to rotate the local aerodynamic force vector
    boolRigidRotation = false;
    fagElem = zeros(12,1);
    for ind = 1 : length( xIntPoints )
        xGauss = lo/2 * (xIntPoints( ind ) + 1); 
        fagElem =  fagElem + lo/2 * wIntPoints(ind) * integAeroForce(xGauss, ddotg, udotWindElem, lo, l, nu, nu11, nu12, nu21, nu22, tl1, tl2, Rr, Ro, vecChordUndef, dimCaracteristic, I3, O3, P, G, EE, L2, L3, boolRigidRotation ) ;
    end
    % express aerodinamic force in ONSAS nomencalture  [force1 moment1 force2 moment2  ...];
    fag = Cambio_Base(fagElem);
end

function integAeroForce = integAeroForce(x, ddotg, udotWindElem, lo, l, nu, nu11, nu12, nu21, nu22, tl1, tl2, Rr, Ro, vecChordUndef, dimCaracteristic, I3, O3, P, G, EE, L2, L3, boolRigidRotation )
    % Compute udot(x) and velWind(x):
    % Shape functions:
    % linear
    N1 = 1 -x/lo            ;
    N2 = x/lo               ;
    % cubic
    N3 = x*(1-x/lo)^2;
    N4 = -(1-x/lo)*(x^2)/lo ;
    N5 = (1-3*x/lo)*(1-x/lo);
    N6 = (3*x/lo-2)*(x/lo)	;
    N7 = N3 +N4             ;
    N8 = N5 +N6 - 1         ;

    % Kinematc variables inside the element
    % auxiliar matrices
    P1 = [  0      0     0   0   0    0     ; ...
            0      0     N3  0   0    N4    ; ...
            0     -N3    0   0   -N4  0     ]; % Eq. 38

    P2 = [ N1       0       0       N2      0       0 ; ...
            0      N5       0       0       N6      0 ; ...
            0      0        N5      0       0       N6] ; % Eq. 39

    N  = [ N1*I3   O3   N2*I3    O3 ] ;

    ul = P1 * [ tl1; tl2 ]              ; % Eq. 38
    H1 = N + P1 * P - 1*skew( ul ) * G' ; % Eq. 59
    H2 = P2*P+G'                        ; % Eq. 72 
    
    % Element velocity inside the element:
    % udotG_loc =   P1(x) * P * EE' * ddotg; %Eq. A.9
    udotG = Rr * H1 * EE' * ddotg; %Eq. 61
    % Global Rotation RgG(x)inside the element:
    thethaRoof  = P2 * [tl1;tl2]      ;% Eq. 39
    Rex         = expon( thethaRoof ) ;
    RgGx        = Rr * Rex * Ro'      ;

    % Wind velocity inside te element
    udotWindG = udotWindElem(1:3) * N1 + udotWindElem(4:6) * N2;
    %Transverse wind velocity inside the element:
    % proyect velocity and chord vector into transverse plane
    VrelG       = udotWindG - udotG  ;
    VpiRelG     = L2 * RgGx' * VrelG ;
    VpiRelGperp = L3 * VpiRelG       ;
    % Calculate relative incidence angle
    if( norm(VpiRelG) == 0)
        fprintf('WARNING: Relative velocity is zero \n')
        td = [0, 0, 1];%define random vector to compute zero force
    else
        td = VpiRelG / norm(VpiRelG);
    end
    % rotate chord vector
    tch = RgGx * vecChordUndef;
    betaRelG =  acos ( dot(tch ,td ) ) ;
    % Aero coefficients:
    C_d =  dragCoefFunction   (betaRelG) ;
    C_l =  liftCoefFunction   (betaRelG) ;
    C_m =  momentCoefFunction (betaRelG) ;    

    % Aero forces
    rhoAire = 1.2;
    fdl =  1/2 * rhoAire * C_d * dimCaracteristic * norm(VpiRelG) * VpiRelG     ; 
    fll =  1/2 * rhoAire * C_l * dimCaracteristic * norm(VpiRelG) * VpiRelGperp ; 
    fal =  fdl + fll;
    ma =  1/2 * rhoAire * C_m * VpiRelG' * VpiRelG * dimCaracteristic * (RgGx*Ro*[1 0 0]'); 

    % Rotate with RG matrix to global rotation matrix:
    RG =   [ RgGx    O3      O3     O3
             O3      RgGx    O3     O3
             O3      O3      RgGx   O3
             O3      O3      O3     RgGx ];    


    % integralTermAeroForceLoc  =   H1' * fal + H2' * ma;  %Eq 78
    if boolRigidRotation
        integAeroForce  =  EE *( H1' * fal + H2' * ma );  %Eq 78
    else
        integAeroForce  =  RG *( H1' * fal + H2' * ma );  %Eq 78
    end
end

function [xIntPoints, wIntPoints] = GaussPointsAndWeights (numGaussPoints )
    if numGaussPoints == 1
        xIntPoints = 0;
        wIntPoints = 2;
    elseif numGaussPoints == 2
        xIntPoints = [ -sqrt(1/3) sqrt(1/3) ];
        wIntPoints = [     1          1     ];        
    elseif numGaussPoints == 3
        xIntPoints = [ -sqrt(3/5)     0  sqrt(3/5)      ];
        wIntPoints = [        5/9	  8/9        5/9    ];
    elseif numGaussPoints == 4
        xIntPoints = [ -sqrt( 3 - 2 * sqrt(6 / 5) ) / sqrt(7),  sqrt( 3 - 2 * sqrt(6 / 5) ) / sqrt(7) ...
                       -sqrt( 3 + 2 * sqrt(6 / 5) ) / sqrt(7),  sqrt( 3 + 2 * sqrt(6 / 5) ) / sqrt(7)   ];
        wIntPoints = [ ( 18 + sqrt(30) ) / 36                   ( 18 + sqrt(30) ) / 36      ... 
                       ( 18 - sqrt(30) ) / 36                   ( 18 - sqrt(30) ) / 36                  ];
    else
        error("The number of gauss cuadrature points introduced are not implemented, only 1 2 3 or 4")
    end
end
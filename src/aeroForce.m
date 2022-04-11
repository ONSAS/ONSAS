function fagElem = aeroForce( elemCoords, elemCrossSecParams,...
                              Ue, Udote, Udotdote, userDragCoef,... 
                              userLiftCoef, userMomentCoef, elemTypeAero,...
                              userFlowVel,geometricNonLinearAero, nextTime ) 
  %Implementation Booleans for internal test it hardcoded into baseBool, essentaly this will compute the aerodynamic force in the rigid confguration
 baseBool = false ;

  %Boolean to compute aerodinamic force with ut = 0 compte aerodynamic force in the reference configuration 
  if ~geometricNonLinearAero 
    Ue = zeros(12,1) ;
  end
  % Nodal Flows:
  if ~isempty(userFlowVel)
    udotFlowNode1 = feval( userFlowVel, elemCoords(1), nextTime ) ; 
    udotFlowNode2 = feval( userFlowVel, elemCoords(4), nextTime ) ; 
    udotFlowElem  = [udotFlowNode1; udotFlowNode2] ;
  else
    error('A userFlowVel field with the name of Flow velocty function file must be defined into analysiSettings struct')
  end
  % Elem reference coordinates:
  xs = elemCoords(:) ;

  % Read aerodinamic profile
  vecChordUndef    = elemTypeAero( 1:3 )'  ;
  dimCaracteristic = norm( vecChordUndef ) ; 

  % Material and cross section props:
  [Area, J, Iyy, Izz, ~ ] = crossSectionProps ( elemCrossSecParams, 0 ) ;

  % Change indexes according to battini's nomenclature
  permutIndxs = [ 1:2:5 2:2:6 ([1:2:5]+6) ([2:2:6]+6) ];
  dg       = Ue      ( permutIndxs ) ;
  ddotg    = Udote   ( permutIndxs ) ;
  ddotdotg = Udotdote( permutIndxs ) ;   
  
  % The rotations matrixes according to Le, Battini 2014 and Battini 2002:
  % rotation global matrices
  tg1 = dg(  4:6  ) ;
  tg2 = dg( 10:12 ) ;
  Rg1 = expon( tg1 ) ;
  Rg2 = expon( tg2 ) ;

  % rotation matrix to reference configuration
  x21 = xs( 4:6 ) - xs( 1:3 ) ;
  d21 = dg( 7:9 ) - dg( 1:3 ) ;
  lo = sqrt( ( x21       )' * ( x21       ) ) ; %
  l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ) ; %
  R0 = beamRefConfRotMat( x21 ) ;

  % rigid rotation matrix:
  % deformed x axis
  e1 = ( x21 + d21 ) / l  ;
  q1 = Rg1 * R0 * [0 1 0]';
  q2 = Rg2 * R0 * [0 1 0]';
  q  = ( q1 + q2 ) / 2; 
  % deformed z local axis
  e3 = cross( e1, q )     ;
  e3 = e3 / norm( e3 )    ; % normalization
  % deformed y local axis
  e2 = cross ( e3, e1 )   ;
  Rr = [ e1 e2 e3 ]       ;
  
  % Compute nus enteries in reference configuration
  q  = Rr' *  q           ;
  q1 = Rr' * q1           ;
  nu = q( 1 ) / q( 2 )    ;
  nu11 = q1( 1 ) / q( 2 ) ;
  nu12 = q1( 2 ) / q( 2 ) ;
  nu21 = 2 * nu - nu11    ;
  nu22 = 2 - nu12         ;

  % local rotations
  if ~baseBool  ;
    Re1 = Rr' * Rg1 * R0 ;
    Re2 = Rr' * Rg2 * R0 ;
    tl1 = logar( Re1 ) ;
    tl2 = logar( Re2 ) ;
  elseif baseBool ;
    Re1 = Rr' * R0 * Rg1 ;
    Re2 = Rr' * R0 * Rg2 ;
    tl1 = logar( Re1 ) ;
    tl2 = logar( Re2 ) ;
  end

  %auxiliar matrix
  I3 = eye(3)     ;
  O3 = zeros(3)   ;
  O1 = zeros(1,3) ;
  
  II=[ O3 I3 O3 O3
       O3 O3 O3 I3 ];

  G=[ 0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
      0   0    1/l     0        0     0  0  0    -1/l     0        0     0
      0  -1/l  0       0        0     0  0  1/l   0       0        0     0 ]' ;    

  P = II - [G'; G'] ;
  %tensor to rotat magnitudes from rigid to global configuration
  EE=[ Rr O3 O3 O3
       O3 Rr O3 O3
       O3 O3 Rr O3
       O3 O3 O3 Rr ] ; 
  
  % auxilair matrix created to proyect transversal velocity
  L2 = [ 0 0 0  
         0 1 0 
         0 0 1 ] ;
  % auxilair matrix created to rotate 90 degrees (this will define the lift force once the drag is computed)
  L3 = expon( [pi/2 0 0] ) ;         

  %Extract points and wieghts for numGausspoints selected
  numGaussPoints = elemTypeAero(4);
  [xIntPoints, wIntPoints] = GaussPointsAndWeights( numGaussPoints ) ;

  %Compute the element fluid force by the equivalent virtual work theory
  fagElem = zeros(12,1) ;
  for ind = 1 : length( xIntPoints )
      %The Gauss integration coordinate is:
      xGauss = lo/2 * (xIntPoints( ind ) + 1) ;
      %Integrate for different cross section inner to the element   
      fagElem =  fagElem ...
                 +lo/2 * wIntPoints(ind) * integAeroForce( xGauss, ddotg, udotFlowElem,... 
                                                          lo, l, nu, nu11, nu12, nu21, nu22, tl1, tl2, Rr, R0,... 
                                                          vecChordUndef, dimCaracteristic,...
                                                          I3, O3, P, G, EE, L2, L3,...
                                                          userDragCoef, userLiftCoef, userMomentCoef ) ;
  end
  % express aerodinamic force in ONSAS nomencalture  [force1 moment1 force2 moment2  ...];
  fagElem = Cambio_Base(fagElem) ;
end

function integAeroForce = integAeroForce( x, ddotg, udotFlowElem,...
                                          lo, l, nu, nu11, nu12, nu21, nu22, tl1, tl2, Rr, R0,... 
                                          vecChordUndef, dimCaracteristic, I3, O3, P, G, EE, L2, L3,...
                                          userDragCoef, userLiftCoef, userMomentCoef )
  
  % Shape functions of Euler Bernoulli element to interpolate displacements and velocites for the cross section:
  % linear
  N1 = 1 -x / lo                           ;
  N2 = x / lo                              ;
  % cubic
  N3 = x * ( 1 - x / lo )^2                ;
  N4 = - ( 1 - x / lo ) * ( x^2 ) / lo     ;
  N5 = ( 1 - 3 * x / lo) * ( 1 - x / lo )  ;
  N6 = ( 3 * x / lo - 2 ) * ( x / lo )	   ;
  N7 = N3 + N4                             ;
  N8 = N5 + N6  -  1                       ;

  % Auxiliar shape function matrices variables for the cross section:
  P1 = [  0   0   0   0   0    0  ; ...
          0   0   N3  0   0    N4 ; ...
          0  -N3  0   0   -N4  0  ] ; % Eq. 38 Le, Battini 2014

  P2 = [  N1  0   0   N2  0    0  ; ...
           0  N5  0   0   N6   0  ; ...
           0  0   N5  0   0    N6 ] ; % Eq. 39 Le, Battini 2014

  N  = [ N1 * I3   O3   N2 * I3    O3 ] ;

  ul = P1 * [ tl1; tl2 ]                ; % Eq. 38 Local axial displacements of the centroid cross section
  H1 = N + P1 * P - 1 * skew( ul ) * G' ; % Eq. 59 auxiliar matrix for linear displacement interpolation (See Le, Battini 2014)
  H2 = P2 * P + G'                      ; % Eq. 72 auxiliar matrix for angular displacement interpolation (See Le, Battini 2014)
  
  % Rotations matrices to compute the flow velocity proyection of the for the cross section:
  % angular local rotation
  thethaRoof  = P2 * [tl1 ; tl2]    ; % Eq. 39 Le, Battini 2014
  % local Rroof rotation matrix is
  Rroofx      = expon( thethaRoof ) ; 
  
  % Kinematic velocities for the generic cross section
  % cross section centroid rigid velocity in global coordiantes:
  udotG = Rr * H1 * EE' * ddotg     ; %Eq. 61 Le, Battini 2014
  % cross section abssolute fluid flow velocity in global coordiantes interpolated with linear shape functions:
  udotFlowG = udotFlowElem(1:3) * N1 + udotFlowElem(4:6) * N2 ;
  
  %Transverse flow velocity of the cross section to compute drag lift and moment:
  %the relative flow velocity in global coordantes is
  VrelG       = udotFlowG - udotG   ;
  %then the proyection of the relative flow velocity on the bending plane is 
  %using deformed configuration coordiantes
  VpiRelG   = L2 * Rroofx' * Rr' * VrelG   ;

  %the perpendicular flow relative velocity proyection in the rigid confugration coordiantes is
  VpiRelGperp = L3 * VpiRelG       ;
  
  %Compute relative incidence angle
  %the chord vector orientation in the deformed confugration coordiantes to compute incidence flow angle is
  tch = ( vecChordUndef / norm( vecChordUndef ) ) ;
  
  % Calculate relative incidence angle in the rigid confugration
  if( norm( VpiRelG) == 0 )
      td = tch;%define tch equal to td if vRel is zero to compute force with zero angle of attack
  else %the drag direction at a generic cross section in rigid configuration coordinates is:
      td = VpiRelG / norm( VpiRelG ) ;
  end
  cosBeta = dot(tch, td) / ( norm(td) * norm(tch) ) ;
  sinBeta = dot( cross(td,tch), [1 0 0] ) / ( norm(td) * norm(tch) ) ;
  betaRelG =  sign( sinBeta ) * acos( cosBeta ) ;
  
  %Check aerodynamic coefficients existence and the load the dimensionless coefficent value if not set to 0:  
  if ~isempty( userDragCoef )
    C_d = feval( userDragCoef, betaRelG ) ;
  else
    C_d = 0 ;
  end
  if ~isempty( userLiftCoef )
    C_l = feval( userLiftCoef, betaRelG ) ; 
  else
    C_l = 0 ;
  end
  if ~isempty( userMomentCoef )
    C_m = feval( userMomentCoef, betaRelG  ) ; 
  else
    C_m = 0 ;
  end

  %The cross section aero forece in rigid configuration coordiantes is:
  % air denisty is hardcoed
  rhoFliud = 1.225 ;
  %drag cross section force vector in rigid or deformed configuration coordiantes
  fdl     =  1/2 * rhoFliud * C_d * dimCaracteristic * norm( VpiRelG) * VpiRelG     ; 
  %lift cross section force vector
  fll     =  1/2 * rhoFliud * C_l * dimCaracteristic * norm( VpiRelG) * VpiRelGperp ; 
  %drag +lift cross section force vector in rigid configuration coordiantes
  fal     =  fdl + fll ;
  ma      =  1/2 * rhoFliud * C_m * VpiRelG' * VpiRelG * dimCaracteristic * ( Rr * Rroofx * [1 0 0]' ) ;

  % Rotate to global coordiantes with EE matrix for rigid configuration formulation
  % compute the inegral term of the current cross section
    integralTermAeroForceRigid  =   H1' * Rroofx * fal + H2' * Rroofx * ma ;  %Similar to Eq 78 with the different that Rr' is not necessary beacous fal is rigid coordinates
    integAeroForce  =  EE *( integralTermAeroForceRigid ) ; %Rotate from rigid to global coordiantes

end

function [xIntPoints, wIntPoints] = GaussPointsAndWeights (numGaussPoints )
%Integration Gauss Points based on https://keisan.casio.com/exec/system/1329114617
  if numGaussPoints == 1
      
      xIntPoints = 0;
     
      wIntPoints = 2;
  
  elseif numGaussPoints == 2
      
      xIntPoints = [ -sqrt(1/3) sqrt(1/3) ];

      wIntPoints = [     1          1     ];        
  
  elseif numGaussPoints == 3
      
      xIntPoints = [ -sqrt(3/5)     0  sqrt(3/5) ];

      wIntPoints = [     5/9	    8/9     5/9    ];

  elseif numGaussPoints == 4
      
      xIntPoints = [ -sqrt( 3 - 2 * sqrt(6 / 5) ) / sqrt(7),  sqrt( 3 - 2 * sqrt(6 / 5) ) / sqrt(7) ...
                     -sqrt( 3 + 2 * sqrt(6 / 5) ) / sqrt(7),  sqrt( 3 + 2 * sqrt(6 / 5) ) / sqrt(7)   ];
      
      wIntPoints = [ ( 18 + sqrt(30) ) / 36                   ( 18 + sqrt(30) ) / 36      ... 
                     ( 18 - sqrt(30) ) / 36                   ( 18 - sqrt(30) ) / 36                  ];
  
  elseif numGaussPoints == 5
      
      xIntPoints = [ -0.9061798459386639927976, -0.5384693101056830910363, 0                       , ...
                      0.5384693101056830910363,  0.9061798459386639927976                                 ] ;
      
      wIntPoints = [  0.2369268850561890875143,  0.4786286704993664680413, 0.5688888888888888888889, ...
                      0.4786286704993664680413,  0.2369268850561890875143                                 ] ;
  
  elseif numGaussPoints == 6
      
      xIntPoints = [ -0.9324695142031520278123, -0.661209386466264513661, -0.2386191860831969086305, ...
                      0.238619186083196908631 ,  0.661209386466264513661, 0.9324695142031520278123        ] ;
      
      wIntPoints = [  0.1713244923791703450403, 0.3607615730481386075698, 0.4679139345726910473899,  ...
                      0.46791393457269104739,   0.3607615730481386075698, 0.1713244923791703450403        ] ;
  
  elseif numGaussPoints == 7
  
      xIntPoints = [ -0.9491079123427585245262, -0.7415311855993944398639, -0.4058451513773971669066, ...
                      0                       ,  0.4058451513773971669066,  0.7415311855993944398639, ...
                      0.9491079123427585245262                                                            ] ;
  
      wIntPoints = [ 0.1294849661688696932706 ,  0.2797053914892766679015, 0.38183005050511894495  , ...
                     0.417959183673469387755  ,  0.38183005050511894495  ,  0.279705391489276667901, ...
                     0.129484966168869693271                                                             ] ;
  
  elseif numGaussPoints == 8

      xIntPoints = [ -0.9602898564975362316836, -0.7966664774136267395916, -0.5255324099163289858177, ...
                     -0.1834346424956498049395, 0.1834346424956498049395, 0.5255324099163289858177, ...
                      0.7966664774136267395916, 0.9602898564975362316836                                  ] ;

      wIntPoints = [  0.1012285362903762591525, 0.2223810344533744705444, 0.313706645877887287338, ...
                      0.3626837833783619829652, 0.3626837833783619829652, 0.313706645877887287338, ...
                      0.222381034453374470544, 0.1012285362903762591525                                   ] ;
  
  elseif numGaussPoints == 9
      
      xIntPoints = [ -0.9681602395076260898356, -0.8360311073266357942994, -0.6133714327005903973087, ...
                     -0.3242534234038089290385,  0                       ,  0.3242534234038089290385, ...
                      0.6133714327005903973087,  0.8360311073266357942994,  0.9681602395076260898356      ] ;
      
      wIntPoints = [  0.0812743883615744119719,  0.1806481606948574040585,  0.2606106964029354623187, ...
                      0.312347077040002840069,   0.330239355001259763165 ,  0.312347077040002840069 , ...
                      0.260610696402935462319,   0.1806481606948574040585,  0.081274388361574411972       ] ;
  
  elseif numGaussPoints == 10
      
      xIntPoints = [ -0.973906528517171720078,  -0.8650633666889845107321, -0.6794095682990244062343, ...
                     -0.4333953941292471907993, -0.1488743389816312108848,  0.1488743389816312108848, ...
                      0.4333953941292471907993,  0.6794095682990244062343,  0.8650633666889845107321, ...
                      0.973906528517171720078                                                             ] ;
      
      wIntPoints = [ 0.0666713443086881375936,   0.149451349150580593146,   0.219086362515982043996, ...
                     0.2692667193099963550912,   0.2955242247147528701739,  0.295524224714752870174, ...
                     0.269266719309996355091,    0.2190863625159820439955,  0.1494513491505805931458, ...
                     0.0666713443086881375936                                                             ] ;  
  else
      error('The number of gauss cuadrature points introduced are not implemented, only between 1 and 10')
  end
end
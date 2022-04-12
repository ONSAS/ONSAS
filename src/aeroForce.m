% Copyright (C) 2022, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

% This function computes the fluid loads within the quasi-steady theory for co-rotational dynamic frame elements proposed by Lee, Battini 2014
function fagElem = aeroForce( elemCoords, elemCrossSecParams,...
                              Ue, Udote, Udotdote, userDragCoef,... 
                              userLiftCoef, userMomentCoef, elemTypeAero,...
                              userFlowVel,geometricNonLinearAero, nextTime ) 

  %Implementation Booleans for internal test, baseBool changes the local angles computation
  baseBool = true ;

  % Boolean to compute fluid force with ut = 0, this should be used if the fluid loads are computed in the reference 
  % configuration and thereafter nonlinear effects are not considered. 
  if ~geometricNonLinearAero 
    Ue = zeros(12,1) ;
  end
  % fluid velocity at the nodes of the element evaluated in deformed configuration (spatial points):
  if ~isempty(userFlowVel)
    udotFlowNode1 = feval( userFlowVel, elemCoords(1) + Ue(1:2:6), nextTime ) ; 
    udotFlowNode2 = feval( userFlowVel, elemCoords(4) + Ue(7:2:12), nextTime ) ;
    % compact them into a single vector for the element 
    udotFlowElem  = [udotFlowNode1; udotFlowNode2] ;
  else
    error('A userFlowVel field with the name of Flow velocity function file must be defined into analysiSettings struct')
  end
  
  % Elem reference coordinates:
  xs = elemCoords(:) ;

  % Load element properties to fluid loads 
  % chord vector
  vecChordUndef    = elemTypeAero( 1:3 )'  ;
  % length of the chord vector 
  dimCaracteristic = norm( vecChordUndef ) ; 

  % Change indexes according to Battini's nomenclature ( ONSAS uses for each node [u1 theta1 u2 theta2 u3 theta3] and Battini [u1 u2 u3 theta1 theta2 theta3] )
  permutIndxs = [ 1:2:5 2:2:6 ([1:2:5]+6) ([2:2:6]+6) ];
  dg       = Ue      ( permutIndxs ) ;
  ddotg    = Udote   ( permutIndxs ) ;
  ddotdotg = Udotdote( permutIndxs ) ;   
  
  % The rotations matrixes according to Le, Battini 2002 and Le 2014:
  % rotation global matrices
  % select thetas
  tg1 = dg(  4:6  )  ;
  tg2 = dg( 10:12 )  ;
  % compute matrices with exp Rodrigue's formula
  Rg1 = expon( tg1 ) ; % Eq(5) J.-M. Battini 2002
  Rg2 = expon( tg2 ) ; % Eq(5) J.-M. Battini 2002

  % rotation reference configuration matrix 
  % element director vector in the reference configuration
  x21 = xs( 4:6 ) - xs( 1:3 ) ;
  % element director vector of the displacements of displacements vector of the element
  d21 = dg( 7:9 ) - dg( 1:3 ) ;
  % reference length of the element
  lo = sqrt( ( x21       )' * ( x21       ) ) ;
  % deformed length of the element
  l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ) ; % Eq(26) J.-M. Battini 2002
  % the function computes RO with e3 such that E3ref maximize the scalar product with E3glob (unless E1ref is E3glob in that case set E1ref = E3glob)
  R0 = beamRefConfRotMat( x21 ) ;

  % rigid rotation matrix:
  % deformed x axis
  e1 = ( x21 + d21 ) / l  ;% Eq(25) J.-M. Battini 2002
  q1 = Rg1 * R0 * [0 1 0]';% Eq(27) J.-M. Battini 2002
  q2 = Rg2 * R0 * [0 1 0]';% Eq(27) J.-M. Battini 2002
  q  = ( q1 + q2 ) / 2    ;% Eq(27) J.-M. Battini 2002 
  % deformed z local axis
  e3 = cross( e1, q )     ;% Eq(28) J.-M. Battini 2002
  e3 = e3 / norm( e3 )    ;% Eq(28) J.-M. Battini 2002 
  % deformed y local axis
  e2 = cross ( e3, e1 )   ;% Eq(28) J.-M. Battini 2002
  Rr = [ e1 e2 e3 ]       ;% Eq(24) J.-M. Battini 2002
  
  % local rotations
  if ~baseBool  ;
    Re1 = Rr' * Rg1 * R0 ;% Eq(30) J.-M. Battini 2002
    Re2 = Rr' * Rg2 * R0 ;% Eq(30) J.-M. Battini 2002
  elseif baseBool ;
    Re1 = Rr' * R0 * Rg1 ;
    Re2 = Rr' * R0 * Rg2 ;
  end
  tl1 = logar( Re1 )   ;% Eq(31) J.-M. Battini 2002
  tl2 = logar( Re2 )   ;% Eq(31) J.-M. Battini 2002

  % Auxiliary matrices computation 
  % compute auxiliary q vectors and nu in rigid coordinates
  q  = Rr' *  q           ;
  q1 = Rr' * q1           ;
  nu = q( 1 ) / q( 2 )    ;% Eq(58) J.-M. Battini 2002
  nu11 = q1( 1 ) / q( 2 ) ;% Eq(58) J.-M. Battini 2002
  nu12 = q1( 2 ) / q( 2 ) ;% Eq(58) J.-M. Battini 2002
  nu21 = 2 * nu - nu11    ;% Eq(58) J.-M. Battini 2002
  nu22 = 2 - nu12         ;% Eq(58) J.-M. Battini 2002

  % identity and null auxiliary matrices
  I3 = eye(3)     ;
  O3 = zeros(3)   ;
  O1 = zeros(1,3) ;
  
  II=[ O3 I3 O3 O3
       O3 O3 O3 I3 ];

  G=[ 0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
      0   0    1/l     0        0     0  0  0    -1/l     0        0     0
      0  -1/l  0       0        0     0  0  1/l   0       0        0     0 ]' ;% Eq(58) J.-M. Battini 2002    

  P = II - [G'; G'] ;% Eq(55) J.-M. Battini 2002

  % matrix E to rotate magnitudes from rigid to global coordinates
  EE=[ Rr O3 O3 O3
       O3 Rr O3 O3
       O3 O3 Rr O3
       O3 O3 O3 Rr ] ; 
  
  % matrix created to project transversal velocity
  L2 = [ 0 0 0  
         0 1 0 
         0 0 1 ] ;
  % auxilair matrix created to rotate 90 degrees (this will define the lift force once the drag is computed)
  L3 = expon( [pi/2 0 0] ) ;         

  % Extract points and weights for numGausspoints selected
  numGaussPoints = elemTypeAero(4);
  [xIntPoints, wIntPoints] = GaussPointsAndWeights( numGaussPoints ) ;

  % Compute the element fluid force by the equivalent virtual work theory
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
  % Transform element fluid forces to ONSAS nomenclature  [force1 moment1 force2 moment2  ...];
  fagElem = Cambio_Base(fagElem) ;
end

function integAeroForce = integAeroForce( x, ddotg, udotFlowElem,...
                                          lo, l, nu, nu11, nu12, nu21, nu22, tl1, tl2, Rr, R0,... 
                                          vecChordUndef, dimCaracteristic, I3, O3, P, G, EE, L2, L3,...
                                          userDragCoef, userLiftCoef, userMomentCoef )
  
  % Shape functions of Euler Bernoulli element to interpolate displacements for a generic cross section:
  % linear
  N1 = 1 -x / lo                           ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N2 = x / lo                              ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  % cubic
  N3 = x * ( 1 - x / lo )^2                ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N4 = - ( 1 - x / lo ) * ( x^2 ) / lo     ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N5 = ( 1 - 3 * x / lo) * ( 1 - x / lo )  ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N6 = ( 3 * x / lo - 2 ) * ( x / lo )	   ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N7 = N3 + N4                             ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N8 = N5 + N6  -  1                       ; % Eq.(37)  T-N Le J.-M. Battini et al 2014

  % Auxiliary shape function matrices variables for the cross section:
  P1 = [  0   0   0   0   0    0  ; ...
          0   0   N3  0   0    N4 ; ...
          0  -N3  0   0   -N4  0  ] ; % Eq.(38)  T-N Le J.-M. Battini et al 2014
  ul = P1 * [ tl1; tl2 ]           ; %% Eq.(38)  T-N Le J.-M. Battini et al 2014

  P2 = [  N1  0   0   N2  0    0  ; ...
           0  N5  0   0   N6   0  ; ...
           0  0   N5  0   0    N6 ] ; % Eq.(39)  T-N Le J.-M. Battini et al 2014
  thethaRoof  = P2 * [tl1 ; tl2]    ; % Eq.(39)  T-N Le J.-M. Battini et al 2014

  N  = [ N1 * I3   O3   N2 * I3    O3 ] ;% Eq.(51)  T-N Le J.-M. Battini et al 2014

  H1 = N + P1 * P - 1 * skew( ul ) * G' ; % Eq.(59) auxiliary matrix for linear displacement interpolation (See Le, Battini 2014)
  H2 = P2 * P + G'                      ; % Eq.(72) auxiliary matrix for angular displacement interpolation (See Le, Battini 2014)
  
  % Rotations matrices to compute the flow velocity projection of the for the cross section:
  % local Rroof rotation matrix is
  Rroofx = expon( thethaRoof ) ; 
  
  % Kinematic velocities for the generic cross section
  % cross section centroid rigid velocity in global coordinates:
  udotG = Rr * H1 * EE' * ddotg     ; % Eq.(61)  T-N Le J.-M. Battini et al 2014
  % cross section absolute fluid flow velocity in global coordinates interpolated with linear shape functions:
  udotFlowG = udotFlowElem(1:3) * N1 + udotFlowElem(4:6) * N2 ;
  
  % Transverse flow velocity of the cross section to compute drag lift and moment:
  % the relative flow velocity in global coordinates is
  VrelG = udotFlowG - udotG   ;
  % then the projection (in t2,t3 plane) of the relative flow velocity in the deformed coordinates is:
  VpiRelG = L2 * Rroofx' * Rr' * VrelG   ;
  % the perpendicular flow relative velocity projection in deformed coordinates is:
  VpiRelGperp = L3 * VpiRelG       ;
  
  % Compute relative incidence angle
  % the chord vector orientation in the deformed coordinates to compute incidence flow angle is:
  tch = ( vecChordUndef / norm( vecChordUndef ) ) ;
  
  % Calculate relative incidence angle in the rigid configuration
  if( norm( VpiRelG) == 0 )
      td = tch;% define tch equal to td if vRel is zero to compute force with zero angle of attack
  else % the drag direction at a generic cross section in deformed coordinates is:
      td = VpiRelG / norm( VpiRelG ) ;
  end
  cosBeta = dot(tch, td) / ( norm(td) * norm(tch) ) ;
  sinBeta = dot( cross(td,tch), [1 0 0] ) / ( norm(td) * norm(tch) ) ;
  betaRelG =  sign( sinBeta ) * acos( cosBeta ) ;
  
  % Check fluid coefficients existence and the load it values:  
  if ~isempty( userDragCoef )
    c_d = feval( userDragCoef, betaRelG ) ;
  else
    c_d = 0 ;
  end
  if ~isempty( userLiftCoef )
    c_l = feval( userLiftCoef, betaRelG ) ; 
  else
    c_l = 0 ;
  end
  if ~isempty( userMomentCoef )
    c_m = feval( userMomentCoef, betaRelG  ) ; 
  else
    c_m = 0 ;
  end

  % The cross section fluid forces in deformed coordinates is:
  % at the moment air density is hardcoed
  rhoFliud = 1.225 ;
  %
  %%%%%%%%%%%%%%FIXME
  %
  % drag cross section force vector in deformed coordinates
  fdl = 1/2 * rhoFliud * c_d * dimCaracteristic * norm( VpiRelG) * VpiRelG     ; 
  % lift cross section force vector in deformed coordinates
  fll = 1/2 * rhoFliud * c_l * dimCaracteristic * norm( VpiRelG) * VpiRelGperp ; 
  % drag + lift cross section force vector in deformed coordinates
  fal = fdl + fll ;
  % torsional moment fluid load in deformed coordinates
  ma = 1/2 * rhoFliud * c_m * VpiRelG' * VpiRelG * dimCaracteristic * ( [1 0 0]' ) ;

  % Compute the element fluid load forces vector in global coordinates
  % compute the integral term of the current cross section in rigid coordinates
    integralTermAeroForceRigid  =   H1' * Rroofx * fal + H2' * Rroofx * ma ;  
  % rotate to global coordinates with EE matrix for rigid configuration formulation
    integAeroForce  =  EE *( integralTermAeroForceRigid ) ; 
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
      error('The number of gauss quadrature points introduced are not implemented, only between 1 and 10')
  end
end
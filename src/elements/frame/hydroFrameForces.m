% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, J. Bruno Bazzano,
% Joaquin Viera, Marcelo Forets, Jean-Marc Battini. 
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
function fagElem = hydroFrameForces( elemCoords,... 
                                     Ue, Udote, Udotdote,... 
                                     aeroCoefs, elemTypeAero, analysisSettings,...
                                     nextTime ) 
  % Check all required parameters are defined
  if isempty(analysisSettings.fluidProps) 
    error(' define correctly row cell analysisSettings.fluidProps = {fluidDensity; fluidViscosity; fluidVelocityFunction } ')
  end
  if isempty(elemTypeAero) 
    error(' define correctly elements.elemTypeAero = [chordVec1 chordVec2 chordVec3 numGauss  ] ')
  end
  if isempty( aeroCoefs ) 
    error(' define correctly row cell of strings elements.aeroCoefs = {dragFunc; liftFunc; momentFunc } ')
  end

  % Implementation Booleans for internal test, baseBool changes the local angles computation
  baseBool = false ;
  % extract fluid properties
  densityFluid       = analysisSettings.fluidProps{1,1} ;
  viscosityFluid = analysisSettings.fluidProps{2,1} ;
  userFlowVel    = analysisSettings.fluidProps{3,1} ;
  % extract nonLinearity in aero force boolean
  geometricNonLinearAero = analysisSettings.geometricNonLinearAero ;
  
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
  dimCharacteristic = norm( vecChordUndef ) ; 

  % Change indexes according to Battini's nomenclature ( ONSAS uses for each node [u1 theta1 u2 theta2 u3 theta3] and Battini [u1 u2 u3 theta1 theta2 theta3] )
  permutIndxs = [ 1:2:5 2:2:6 ([1:2:5]+6) ([2:2:6]+6) ];
  dg       = Ue      ( permutIndxs ) ;
  ddotg    = Udote   ( permutIndxs ) ;
  ddotdotg = Udote   ( permutIndxs ) ;
  
  % The rotations matrixes according to Le, Battini 2002 and Le 2014:
  % rotation global matrices
  % select thetas
  tg1 = dg(  4:6  ) ;
  tg2 = dg( 10:12 ) ;
 % compute matrices with exp Rodrigue's formula
  Rg1 = expon( tg1 ) ;% Eq(5) J.-M. Battini 2002
  Rg2 = expon( tg2 ) ;% Eq(5) J.-M. Battini 2002

  % rotation reference configuration matrix 
  % element director vector in the reference configuration
  x21 = xs( 4:6 ) - xs( 1:3 ) ;
  % element director vector of the displacements of displacements vector of the element
  d21 = dg( 7:9 ) - dg( 1:3 ) ;
  % reference length of the element
  lo = sqrt( ( x21       )' * ( x21       ) ) ; 
  % deformed length of the element
  l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ) ; 
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
  
  % Auxiliary matrices computation 
  % compute auxiliary q vectors and nu in rigid coordinates
  q  = Rr' *  q           ;
  q1 = Rr' * q1           ;
  nu = q( 1 ) / q( 2 )    ;% Eq(58) J.-M. Battini 2002
  nu11 = q1( 1 ) / q( 2 ) ;% Eq(58) J.-M. Battini 2002
  nu12 = q1( 2 ) / q( 2 ) ;% Eq(58) J.-M. Battini 2002
  nu21 = 2 * nu - nu11    ;% Eq(58) J.-M. Battini 2002
  nu22 = 2 - nu12         ;% Eq(58) J.-M. Battini 2002

  % local rotations
  if ~baseBool ;
    Re1 = Rr' * Rg1 * R0 ;% Eq(30) J.-M. Battini 2002
    Re2 = Rr' * Rg2 * R0 ;% Eq(30) J.-M. Battini 2002
  elseif baseBool ;
    Re1 = Rr' * R0 * Rg1 ;
    Re2 = Rr' * R0 * Rg2 ;
  end
  tl1 = logar( Re1 ) ; % Eq(31) J.-M. Battini 2002
  tl2 = logar( Re2 ) ;% Eq(31) J.-M. Battini 2002
  
  % identity and null auxiliary matrices
  I3 = eye(3)     ;
  O3 = zeros(3)   ;
  
  II=[ O3 I3 O3 O3
       O3 O3 O3 I3 ];

  G=[ 0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
      0   0    1/l     0        0     0  0  0    -1/l     0        0     0
      0  -1/l  0       0        0     0  0  1/l   0       0        0     0 ]' ;% Eq(58) J.-M. Battini 2002     

  P = II - [G'; G'] ; % Eq(55) J.-M. Battini 2002
  % tensor to rotate magnitudes from rigid to global configuration
  EE=[ Rr O3 O3 O3
       O3 Rr O3 O3
       O3 O3 Rr O3
       O3 O3 O3 Rr ] ;% Eq.(30)  T-N Le J.-M. Battini et al 2014
  
  % auxiliary matrix created to project transversal velocity
  L2 = [ 0 0 0  
         0 1 0 
         0 0 1 ] ;
  % auxillary matrix created to rotate 90 degrees (this will define the lift force once the drag is computed)
  L3 = expon( [pi/2 0 0] ) ;         

  % Extract points and weights for numGausspoints selected
  numGaussPoints = elemTypeAero(4);
  [xIntPoints, wIntPoints] = GaussPointsAndWeights( numGaussPoints ) ;
  % WOM computation call for cases with VIVbool equal to true
  if exist('VIVBool')~=0
    if VIVBool == 'true'
      % extract the accelerations and the velocities of nodes
      % node 1
      udotFrame1    = Udote(1:2:6)      ;
      udotdotFrame1 = Udotdote(1:2:6)   ;
      % node 2
      udotFrame2    = Udote(7:2:end)    ;
      udotdotFrame2 = Udotdote(7:2:end) ;

      % projected velocities at nodes 1 and 2 
      %node 1 
      [VpiRel1, VpiRelPerp1, Vrel1] = computeVpiRels( udotFlowNode1, udotFrame1, Re1, Rr, L2, L3 ) ;
      %node 2 
      [VpiRel2, VpiRelPerp2, Vrel2] = computeVpiRels( udotFlowNode2, udotFrame2, Re2, Rr, L2, L3 ) ;
      % directions of lift forces
      tl1 = VpiRelPerp1 / norm( VpiRelPerp1 ) ;
      tl2 = VpiRelPerp2 / norm( VpiRelPerp2 ) ;
      % evaluate cl0(at the moment cl0 cannot depend on Reynolds and relative incidence angle)
      userLiftCoef   = aeroCoefs{2} ;
      if ~isempty( userLiftCoef )
        c_l0 = feval( userLiftCoef, 0, 0  ) ;
      else
        error('for the VIV force you need to define a userLiftCoef returning c_l0')
      end
      % compute
      q = WOMV1(VpiRel1, VpiRel2, udotdotFrame1, udotdotFrame2, tl1, tl2, dimCharacteristic, c_l0, nextTime, analysiSettings.deltaT ) ; 
    end
    else
      q = 2 ;
  end
  % Compute the element fluid force by the equivalent virtual work theory
  fagElem = zeros(12,1) ;
  for ind = 1 : length( xIntPoints )
      %The Gauss integration coordinate is:
      xGauss = lo/2 * (xIntPoints( ind ) + 1) ;
      %Integrate for different cross section inner to the element   
      fagElem =  fagElem ...
                 +lo/2 * wIntPoints(ind) * integAeroForce( xGauss, ddotg, udotFlowElem,... 
                                                          lo, tl1, tl2, Rr, ... 
                                                          vecChordUndef, dimCharacteristic,...
                                                          I3, O3, P, G, EE, L2, L3,...
                                                          aeroCoefs, densityFluid, viscosityFluid, q ) ;
  end
  % express aerodynamic force in ONSAS nomenclature  [force1 moment1 force2 moment2  ...];
  fagElem = Cambio_Base(fagElem) ;
end

function integAeroForce = integAeroForce( x, ddotg, udotFlowElem,...
                                          lo, tl1, tl2, Rr, ... 
                                          vecChordUndef, dimCharacteristic, I3, O3, P, G, EE, L2, L3,...
                                          aeroCoefs, densityFluid, viscosityFluid, q )
  
  % Shape functions of Euler Bernoulli element to interpolate displacements and velocites for the cross section:
  % linear
  N1 = 1 -x / lo                           ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N2 = x / lo                              ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  % cubic
  N3 = x * ( 1 - x / lo )^2                ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N4 = - ( 1 - x / lo ) * ( x^2 ) / lo     ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N5 = ( 1 - 3 * x / lo) * ( 1 - x / lo )  ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N6 = ( 3 * x / lo - 2 ) * ( x / lo )	   ; % Eq.(37)  T-N Le J.-M. Battini et al 2014

  % Auxiliary shape function matrices variables for the cross section:
  P1 = [  0   0   0   0   0    0  ; ...
          0   0   N3  0   0    N4 ; ...
          0  -N3  0   0   -N4  0  ] ; % Eq.(38)  T-N Le J.-M. Battini et al 2014

  P2 = [  N1  0   0   N2  0    0  ; ...
           0  N5  0   0   N6   0  ; ...
           0  0   N5  0   0    N6 ] ; % Eq.(39)  T-N Le J.-M. Battini et al 2014

  N  = [ N1 * I3   O3   N2 * I3    O3 ] ; % Eq.(51)  T-N Le J.-M. Battini et al 2014

  ul = P1 * [ tl1; tl2 ]                ; % Eq.(38)  T-N Le J.-M. Battini et al 2014
  H1 = N + P1 * P - 1 * skew( ul ) * G' ; % Eq.(59)  T-N Le J.-M. Battini et al 2014
  H2 = P2 * P + G'                      ; % Eq.(72)  T-N Le J.-M. Battini et al 2014
  
  % Rotations matrices to compute the flow velocity projection of the for the cross section:
  % angular local rotation
  thethaRoof  = P2 * [tl1 ; tl2] ; % Eq. 39 Le, Battini 2014
  % local Rroof rotation matrix is
  Rroofx      = expon( thethaRoof ) ; 
 
  % Kinematic velocities for the generic cross section
  % cross section centroid rigid velocity in global coordinates:
  udotG = Rr * H1 * EE' * ddotg ; % Eq.(61)  T-N Le J.-M. Battini et al 2014
  % cross section absolute fluid flow velocity in global coordinates interpolated with linear shape functions:
  udotFlowG = udotFlowElem(1:3) * N1 + udotFlowElem(4:6) * N2 ;
  
  % Relative, perpendicular and projected  flow velocity of the cross section to compute drag lift and moment:
  [VpiRelG, VpiRelGperp, VrelG] = computeVpiRels(udotFlowG, udotG, Rroofx, Rr, L2, L3 )  ;
  
  % Compute relative incidence angle
  % the chord vector orientation in the deformed coordinates to compute incidence flow angle is:
  tch = (vecChordUndef / norm( vecChordUndef )) ;

  % Calculate relative incidence angle in the deformed configuration
  if( norm( VpiRelG) == 0 )
      td = tch ;%define tch equal to td if vRel is zero to compute force with zero angle of attack
  else % the drag direction at a generic cross section in deformed coordinates is:
      td = VpiRelG / norm( VpiRelG ) ;
  end
  cosBeta  = dot(tch, td) / ( norm(td) * norm(tch) ) ;
  sinBeta  = dot( cross(td,tch), [1 0 0] ) / ( norm(td) * norm(tch) ) ;
  betaRelG = sign( sinBeta ) * acos( cosBeta ) ;
  
  % Delete spaces
  userDragCoef   = aeroCoefs{1} ;
  userLiftCoef   = aeroCoefs{2} ;
  userMomentCoef = aeroCoefs{3} ;

  % Computation of Renynolds number
  Re = norm(udotFlowG) * dimCharacteristic / viscosityFluid ;
  
  % Check fluid coefficients existence and the load it values if not set 0:  
  if ~isempty( userDragCoef )
    c_d = feval( userDragCoef, betaRelG, Re  ) ;
  else
    c_d = 0 ;
  end
  if ~isempty( userLiftCoef )
    c_l = feval( userLiftCoef, betaRelG, Re  ) ; 
  else
    c_l = 0 ;
  end
  if ~isempty( userMomentCoef )
    c_m = feval( userMomentCoef, betaRelG, Re) ; 
  else
    c_m = 0 ;
  end
  % The cross section fluid forces in deformed coordinates is:
  % drag cross section force vector in deformed coordinates
  fdl =  1/2 * densityFluid * c_d * dimCharacteristic * norm( VpiRelG) * VpiRelG     ; 
  % lift cross section force vector in deformed coordinates
  fll =  1/2 * densityFluid * c_l * dimCharacteristic * norm( VpiRelG) * VpiRelGperp * q / 2 ; %note that if there is no VIV effect q is 2
  % drag + lift cross section force vector in deformed coordinates
  fal =  fdl + fll ;
 % torsional moment fluid load in deformed coordinates
  ma =  1/2 * densityFluid * c_m * VpiRelG' * VpiRelG * dimCharacteristic * ( [1 0 0]' ) ;

  % Compute the element fluid load forces vector in global coordinates
  % compute the integral term of the current cross section in rigid coordinates
  integralTermAeroForceRigid  =   H1' * Rroofx * fal + H2' * Rroofx * ma ;  
  % rotate to global coordinates with EE matrix for rigid configuration formulation
  integAeroForce  =  EE *( integralTermAeroForceRigid ) ; %Rotate from rigid to global coordinates

end
% This function return the relative projected velocity
function [VpiRel, VpiRelPerp, VrelG] = computeVpiRels( udotFlow, udotFrame, Rroof, Rr, L2, L3 )
  % the relative velocity is:
  VrelG = udotFlow - udotFrame ;
  % then the projection (in t2,t3 plane) of the relative flow velocity in the deformed coordinates is:
  VpiRel = L2 * Rroof' * Rr' * VrelG ;
  % the perpendicular flow relative velocity projection in deformed coordinates is:
  VpiRelPerp = L3 * VpiRel ;
end
% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
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
function  [ fs, ks ]= frame_inertial_force( elemCoords, ...
                        elemCrossSecParams, elemConstitutiveParams, ...
                        Ue, Udote, Udotdote, elemrho, massMatType ) ;
  global massratio;
  % element coordinates
  xs = elemCoords(:) ;

  % ----- material and geometric params ------
  E   = elemConstitutiveParams(2) ;
  nu  = elemConstitutiveParams(3) ;
  G   = E/(2*(1+nu)) ;
  rho = elemrho ;
  if ~isempty( massratio ) % massratio = rho_structure/rho_fluid
    rho = rho * (1+(1/massratio)) ;
  end
  % ----- extract cross section properties ---
  [Area, J, Iyy, Izz, Jrho] = crossSectionProps ( elemCrossSecParams, rho ) ;
  % ------------------------------------------

  % compute corotational matrices rotation 
  [R0, Rr, Rg1, Rg2, Rroof1, Rroof2] = corotRotMatrices( Ue, elemCoords ) ;

  % --- global kinematics ---
  % permut indexes according to Battini's nomenclature
  dg = switchToBattiniNom( Ue ) ;
  % global thetas
  tg1 = dg(  4:6  ) ;
  tg2 = dg( 10:12 ) ;
  % global acel and vels
  ddotg    = switchToBattiniNom( Udote ) ;
  ddotdotg = switchToBattiniNom( Udotdote ) ;
  % -------------------------------

  % length and coords of the element  
  [x21, d21, l, l0] = corotLenCoords(xs ,dg) ;

  % --- local displacements ---
  % axial displacement
  u   = l - l0 ;
  % local rotations
  tl1 = logar( Rroof1 ) ;
  tl2 = logar( Rroof2 ) ;  
  locDisp = [ u tl1' tl2' ] ;
  % -------------------------------

  % --- auxiliary vector and matrices  ---
  % aux zero matrices 
  [I3, O3, O1, II] = corotZeros() ;

  % auxiliary q base 
  q1g = Rg1 * R0 * [0 1 0]' ;
  q2g = Rg2 * R0 * [0 1 0]' ;
  qg  = ( q1g + q2g ) / 2   ;

  [nu, nu11, nu12, nu21, nu22, e1, e2, e3, r, Gaux, P, EE ] = corotVecMatAuxStatic(...
                                                                  R0, Rr, Rg1, Rg2, l, II, O3, O1);
                

if strcmp( massMatType, 'consistent' )
    sumInterForce  = zeros (12, 1 ) ;
    sumGyro        = zeros (12    ) ;
    sumMass        = zeros (12    ) ;

    % Compute interial force by 3 quadrature gauss points
    [xIntPoints, wIntPoints] = gaussPointsAndWeights (3) ;

    % Compute internalForce
    for ind = 1 : length( xIntPoints )

    xGauss = l0/2 * (xIntPoints( ind ) + 1) ;

    [interTermInertialForce, interTermMassMatrix, interTermGyroMatrix ] = interElementBeamForces ( ... 
        xGauss, l0, l, tl1, tl2, ddotg, ddotdotg, r, P, EE, I3, O3, O1, Rr, R0, Jrho, rho, Area, Gaux) ;

    sumInterForce = sumInterForce ...
        + l0/2 * wIntPoints( ind ) * interTermInertialForce ;
    %
    sumGyro = sumGyro ...
        + l0/2 * wIntPoints( ind ) * interTermGyroMatrix  ;
    %
    sumMass = sumMass ...
        + l0/2 * wIntPoints( ind ) * interTermMassMatrix ;
    end

    Fine       = EE * sumInterForce ;
    GyroMatrix = EE * sumGyro * EE' ;
    MassMatrix = EE * sumMass * EE' ;

    % % Add Bt Matrix for expon angular update
    % Bt=[I3   O3       O3      O3
    %     O3 inv(Dg1)'    O3      O3
    %     O3     O3      I3      O3
    %     O3     O3      O3      inv(Dg2)' ];
    % MassMatrix = MassMatrix * Bt ;
    % GyroMatrix = GyroMatrix * Bt ;

    Fine       = swtichToONSASBase(Fine); % Format [f1 m1 ...];
    GyroMatrix = swtichToONSASBase(GyroMatrix); % Format [u1 theta1 u2 theta2 u3 theta3];
    MassMatrix = swtichToONSASBase(MassMatrix); % Format [u1 theta1 u2 theta2 u3 theta3];

    %~ Fine
    fs{3} = Fine ;

    ks{2} = GyroMatrix ;
    ks{3} = MassMatrix ;

elseif strcmp( massMatType, 'lumped' )
    Me = sparse(12,12)                                      ;
    Me (1:2:end, 1:2:end) = rho * Area * l0 * 0.5 * eye(6)  ;
    Fine = Me * Udotdote                                    ;

    fs{3} = Fine      ;
    ks{2} = zeros(12) ;
    ks{3} = Me        ;
else
    error('the massMatType field into the elements struct must be or consistent or lumped' )
end %endIfConsistentBoolean

end%endFunction

function [IntegrandoForce, IntegrandoMassMatrix, IntegrandoGyroMatrix ] = interElementBeamForces (...
x, l0, l, tl1, tl2, ddotg, ddotdotg, r, P, EE, I3, O3, O1, Rr, Ro, Jrho, rho, Area, Gaux )
    % Bernoulli weight function
    [N1, N2, N3, N4, N5, N6, N7, N8] = bernoulliInterpolWeights(x, l0) ;

    % Auxiliary matrices 
    [P1, P2, N, N1, N2] = corotVecMatAuxDyn( N1, N2, N3, N4, N5, N6, N7, N8, tl1, tl2, Gaux, I3, O3, P ) ;

    % --- local displacements of a generic cross section ---
    ul = P1 * [ tl1; tl2 ]                   ; % Eq.(38)  T-N Le J.-M. Battini et al 2014
    % angular local rotation
    thethaRoof  = P2 * [tl1 ; tl2] ; % Eq. 39 Le, Battini 2014
    % local Rroof rotation matrix is
    Rroofx      = expon( thethaRoof ) ; 
    % rigid velocity
    wdoter  = Gaux' * EE' * ddotg ;% Eq. 65

    H1  = N + P1 * P - 1 * skew( ul ) * Gaux' ;

    A1  = [ O1           O1   O1         O1 ;
            0   -1  0    O1   0  1  0	 O1 ;
            0   0 	-1   O1	  0  0  1    O1 ] ; %Eq. A.4

    udotl   = P1 * P * EE' * ddotg ; %Ec A.9

    % Auxiliar matrix to veclocites
    H1dot   = N7 / (l^2) * A1 * (r' * ddotg) - skew( udotl ) * Gaux' ; %Ec A.8

    % compute skew only one time
    skewWdoter = skew(wdoter);

    ET      = [ skewWdoter   O3         O3 		      O3	          ;
                O3		       skewWdoter O3   		    O3	          ;
                O3			     O3 			  skewWdoter  O3      ;
                O3			     O3  			  O3          skewWdoter   ];

    C1      = skewWdoter * H1 + H1dot - H1 * ET; % Ec  66

    % Linear and angular veclocites and acelerations
    udot    = Rr * H1 * EE' * ddotg; %Ec 61
    udotdot = Rr * H1 * EE' * ddotdotg + Rr * C1 * EE' * ddotg; % Ec 67

    H2      = P2 * P + Gaux'; %Ec 72 se puede usar para comprobar con ec A.10
    wdot  = Rr * H2 * EE' * ddotg;%Ec74

    A2    = [  O1                 O1         O1            O1;
                0      0      1   O1         0 0 -1        O1;
                0      -1     0   O1         0 1  0        O1] ;%Ec A.12
    H2dot  = N8 / l^2 * A2 * (r'*ddotg) ;%Ec A.14

    C2     = skewWdoter * H2 + H2dot - H2 * ET ;%Ec 76

    wdotdot= Rr * H2 * EE' * ddotdotg  + Rr * C2 * EE' * ddotg ;%Ec 77
    
    % Compute global rotation Rg(x)
    thethaRoof  = P2 * [tl1;tl2]    ;% Ec 39
    Rroofx      = expon(thethaRoof) ; %Ec 19 elevado en ambos lados
    Rgx  	    = Rr * Rroofx * Ro' ;

    % Compute dyadic tensor
    Irho  = Rgx * Ro * Jrho * (Rgx*Ro)' ; %Ec 45
    Irhoe = Rr' * Irho * Rr             ; %Ec 80
    
    % Calculate integral Force
    IntegrandoForce  =      H1'* Rr' * Area * rho * udotdot ...
                            + H2' * Rr' * ( Irho*wdotdot + skew(wdot) * Irho * wdot ) ;  %Eq 78

    IntegrandoMassMatrix  = H1' * Area *rho * H1 + H2' * Irhoe * H2 ;

    % Compute C3 and C4
    h1 = H1 * ddotg ; %Eq B6
    h2 = H2 * ddotg ;

    rElem = [ [-1 0 0]   O1  [1 0 0] O1]; %Ec B10

    F1    = [skew(ddotg(1:3))' skew(ddotg(4:6))' skew(ddotg(7:9))' skew(ddotg(10:12))']' ;

    C3    = -skew(h1) * Gaux'  + (N7 / l^2) * A1 *(ddotg * rElem)...
                  +skewWdoter * P1 * P + H1 * F1 * Gaux'; % B13

    C4  = -skew(h2)*Gaux' + ( N8 / l^2 ) * A2 * ddotg * rElem + H2 * F1 * Gaux'; %B14

    % Compute Gyroscopic Matrix
    IntegrandoGyroMatrix  = H2' * ( ( skewWdoter * Irhoe ) - skew( Irhoe * wdoter) ) * H2 ...
                            + H1' * Area * rho*(C1 + C3)  + H2'*Irhoe*(C2+C4) ; %Ec88

end

% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% J. Bruno Bazzano, Joaquin Viera, Marcelo Forets, Jean-Marc Battini. 
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

function  [ fs, ks ]= elementFrameMassForce( ...
  elemCoords, elemCrossSecParams, elemConstitutiveParams, Ue, Udote, Udotdote, elemrho, massMatType ) ;

  % element coordinates
  xs = elemCoords(:) ;
  % --- material constitutive  params ---
  rho = elemrho ;
  % -------------------------------

  % -------------------------------
  [Area, J, Iyy, Izz, Jrho ] = crossSectionProps ( elemCrossSecParams, elemrho ) ;

  % auxiliar matrices
  I3 = eye(3)     ;
  O3 = zeros(3)   ;
  O1 = zeros(1,3) ;

  permutIndxs = [ 1:2:5 2:2:6 ([1:2:5] + 6) ([2:2:6] + 6) ] ;

  dg       = Ue      ( permutIndxs ) ;
  ddotg    = Udote   ( permutIndxs ) ;
  ddotdotg = Udotdote( permutIndxs ) ;

  % global thetas
  tg1 = dg(  4:6  ) ;
  tg2 = dg( 10:12 ) ;

  % rotation matrices
  Rg1 = expon( tg1 ) ;
  Rg2 = expon( tg2 ) ;

  x21 = xs( 4:6 ) - xs( 1:3 ) ;
  d21 = dg( 7:9 ) - dg( 1:3 ) ;

  lo = sqrt( ( x21       )' * ( x21       ) ) ; %
  l  = sqrt( ( x21 + d21 )' * ( x21 + d21 ) ) ; %

  % rotation matrix to reference configuration
  Ro = beamRefConfRotMat( x21 ) ;


  % --- rigid rotation ---
  % deformed x axis
  e1 = ( x21 + d21 ) / l   ;

  q1 = Rg1 * Ro * [0 1 0]' ;
  q2 = Rg2 * Ro * [0 1 0]' ;
  q  = ( q1 + q2 ) / 2     ;

  % deformed z local axis
  e3 = cross(e1, q)    ;
  e3 = e3 / norm( e3 ) ; % normalization

  % deformed y local axis
  e2 = cross (e3, e1);

  % rotation matrix
  Rr = [ e1 e2 e3 ] ;
  % -------------------

  % --- local displacements ---
  % axial displacement
  u  = l - lo;

  % local rotations
  % Rr * Re1 * u = Rg1 * R0 * u
  Re1 = Rr' * Rg1 * Ro;
  Re2 = Rr' * Rg2 * Ro;

  tl1 = logar( Re1 ) ;
  tl2 = logar( Re2 ) ;
  % -----------------------

  % --- auxiliary matrices ---

  q  = Rr' *  q ;
  q1 = Rr' * q1 ;

  nu    = q(1)  / q(2) ;
  nu11  = q1(1) / q(2) ;
  nu12  = q1(2) / q(2) ;
  nu21  = 2*nu - nu11  ;
  nu22  = 2 - nu12     ;

  % transformation to the global coordinates an axial vector
  r = [ -e1' O1  e1' O1 ]' ;


  G=[0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
     0   0    1/l     0        0     0  0  0    -1/l     0        0     0
     0  -1/l  0       0        0     0  0  1/l   0       0        0     0]';

  II=[O3 I3 O3 O3
      O3 O3 O3 I3];

  P = II - [G'; G'] ;

  % rotation to global coords tensor
  EE=[Rr O3 O3 O3
      O3 Rr O3 O3
      O3 O3 Rr O3
      O3 O3 O3 Rr];


if strcmp( massMatType, 'consistent' )
    sumInterForce  = zeros (12, 1 ) ;
    sumGyro        = zeros (12    ) ;
    sumMass        = zeros (12    ) ;

    % Compute interial force by quadrature
    xIntPoints = [ -sqrt(3/5)     0  sqrt(3/5)  ] ;
    wIntPoints = [        5/9	  8/9        5/9  ] ;

    % Compute internalForce
    for ind = 1 : length( xIntPoints )

    xGauss = lo/2 * (xIntPoints( ind ) + 1) ;

    [interTermInertialForce, interTermMassMatrix, interTermGyroMatrix ] = interElementBeamForces ( ... 
        xGauss, lo, l, tl1, tl2, ddotg, ddotdotg, r, P, EE, I3, O3, O1, Rr, Ro, Jrho, rho, Area, G) ;

    sumInterForce = sumInterForce ...
        + lo/2 * wIntPoints( ind ) * interTermInertialForce ;
    %
    sumGyro = sumGyro ...
        + lo/2 * wIntPoints( ind ) * interTermGyroMatrix  ;
    %
    sumMass = sumMass ...
        + lo/2 * wIntPoints( ind ) * interTermMassMatrix ;
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

    Fine       = Cambio_Base(Fine); % En formato [f1 m1 ...];
    GyroMatrix = Cambio_Base(GyroMatrix); % En formato [u1 theta1 u2 theta2 u3 theta3];
    MassMatrix = Cambio_Base(MassMatrix); % En formato [u1 theta1 u2 theta2 u3 theta3];

    %~ Fine
    fs{3} = Fine ;

    ks{2} = GyroMatrix ;
    ks{3} = MassMatrix ;

elseif strcmp( massMatType, 'lumped' )
    Me = sparse(12,12)                                      ;
    Me (1:2:end, 1:2:end) = rho * Area * lo * 0.5 * eye(6)  ;
    Fine = Me * Udotdote                                    ;

    fs{3} = Fine      ;
    ks{2} = zeros(12) ;
    ks{3} = Me        ;
else
    error('the massMatType field into the elements struct must be or consistent or lumped' )
end %endIfConsistentBoolean

end%endFunction

function [IntegrandoForce, IntegrandoMassMatrix, IntegrandoGyroMatrix ] = interElementBeamForces ( x, lo, l, tl1, tl2, ddotg, ddotdotg, r, P, EE, I3, O3, O1, Rr, Ro, Jrho, rho, Area, G )

    % Compute Bernoulli linear interpolation fucntions:
    % linear
    N1 = 1 - x / lo	                     ;
    N2 = x/lo	                            ;
    % cubic
    N3 = x * ( 1 - x / lo )^2	       ;
    N4 = -(1 - x / lo ) * ( x^2 ) / lo    ;
    N5 = ( 1- 3* x / lo ) * ( 1 - x / lo );
    N6 = ( 3 * x/lo - 2 )*( x/lo )        ;
    N7 = N3 + N4 		              ;
    N8 = N5 + N6 -1		              ;

    % Compute auxiliar matrixes:
    P1  =[ 0    0       0       0       0        0 ; ...
           0    0       N3      0       0       N4 ; ...
           0    -N3     0       0       -N4      0 ] ;

    P2  =[ N1   0       0       N2      0       0 ; ...
           0    N5      0       0       N6      0 ; ...
           0    0       N5      0       0       N6 ] ;

    N   =[ N1*I3    O3  N2*I3   O3 ];

    % local displacements
    ul  = P1 * [ tl1; tl2 ] ; % Eq. 38 local disp
    wdoter  = G' * EE' * ddotg ;% Eq. 65

    H1  = N + P1 * P - 1 * skew( ul ) * G' ;

    A1  =[  O1           O1   O1         O1 ;
            0   -1  0    O1 	 0  1  0	O1 ;
            0   0 	-1  O1	 0  0  1    O1 ] ; %Eq. A.4

    udotl   = P1 * P * EE' * ddotg ; %Ec A.9

    % Auxiliar matrix to veclocites
    H1dot   = N7 / (l^2) * A1 * (r' * ddotg) - skew( udotl ) * G' ; %Ec A.8

    % compute skew only one time
    skewWdoter = skew(wdoter);
    ET      = [ skewWdoter   O3         O3 		      O3	          ;
                O3		       skewWdoter O3   		    O3	          ;
                O3			     O3 			  skewWdoter  O3	          ;
                O3			     O3  			  O3          skewWdoter   ];

    C1      = skewWdoter*H1 + H1dot - H1*ET; % Ec  66

    % Linear and angular veclocites and acelerations
    udot    = Rr * H1 * EE' * ddotg; %Ec 61
    udotdot = Rr * H1 * EE' * ddotdotg + Rr * C1 * EE' * ddotg; % Ec 67

    H2      = P2 * P + G'; %Ec 72 se puede usar para comprobar con ec A.10
    wdot  = Rr * H2 * EE' * ddotg;%Ec74

    A2    = [  O1                 O1         O1            O1;
                0      0      1   O1         0 0 -1        O1;
                0      -1     0   O1         0 1  0        O1] ;%Ec A.12
    H2dot  = N8 / l^2 * A2 * (r'*ddotg) ;%Ec A.14

    C2     = skewWdoter*H2 + H2dot - H2*ET ;%Ec 76

    wdotdot= Rr*H2*EE'*ddotdotg  + Rr*C2*EE'*ddotg ;%Ec 77
    % Compute global rotation Rg(x)
    thethaRoof  = P2*[tl1;tl2];% Ec 39
    Rex         = expon(thethaRoof); %Ec 19 elevado en ambos lados
    Rgx  	      = Rr*Rex*Ro';

    % Compute dyadic tensor
    Irho	  = Rgx * Ro * Jrho * (Rgx*Ro)'; %Ec 45
    Irhoe       = Rr' * Irho * Rr;   	 		%Ec 80

     % Calculate integral Force
    IntegrandoForce  =      H1'* Rr' * Area * rho * udotdot ...
                            + H2' * Rr' * ( Irho*wdotdot + skew(wdot) * Irho * wdot ) ;  %Eq 78

    IntegrandoMassMatrix  = 1*H1'*Area*rho*H1 + 1*H2'*Irhoe*H2;

    % Compute C3 and C4
    h1 = H1 * ddotg ; %Eq B6
    h2 = H2 * ddotg ;

    rElem = [ [-1 0 0]   O1  [1 0 0] O1]; %Ec B10

    % Could be worng
    F1    = [skew(ddotg(1:3))' skew(ddotg(4:6))' skew(ddotg(7:9))' skew(ddotg(10:12))']'; %Chequear con los nodales

    C3    = -skew(h1) * G'  + (N7 / l^2) * A1 *(ddotg * rElem)...
                  +skewWdoter * P1 * P + H1 * F1 * G'; % B13

    C4  = -skew(h2)*G' + ( N8 / l^2 )*A2*ddotg*rElem + H2 * F1 * G'; %B14

    % Compute Gyroscopic Matrix
    IntegrandoGyroMatrix  =    H2' * ( ( skewWdoter * Irhoe ) - skew( Irhoe * wdoter) ) * H2 ...
                                + H1' * Area*rho*(C1 + C3)  + H2'*Irhoe*(C2+C4) ; %Ec88

end

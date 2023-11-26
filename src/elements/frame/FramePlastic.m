% Copyright 2023, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini,
% Sergio A. Merlino.
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
 
% =========================================================================

% Euler-Bernoulli element with embeded discontinuity
% Numerical modeling of softening hinges in thin Euler–Bernoulli beams
% Francisco Armero, David Ehrlich / University of California, Berkeley

% =========================================================================

function [ fs, ks, finteLocalCoor ] = FramePlastic(elemCoords, elemCrossSecParams, massMatType, density, hyperElasModel, hyperElasParams, Ut, Udotdotte, intBool, matFintBool, elem)
  
  ndofpnode = 6 ;
  
  % --- material constit params ---
	E   = hyperElasParams(1) ;
	nu  = hyperElasParams(2) ;
	G   = E/(2*(1+nu)) ;
	
	[A, J, Iy, Iz] = crossSectionProps ( elemCrossSecParams, density ) ;
	
  % --- elem lengths and rotation matrix
	[ local2globalMats, l ] = beamParameters( elemCoords ) ;
	R = RotationMatrix(ndofpnode, local2globalMats) ;
	
  % temporary
  % ------------------------
  elemReleases = [0 0 0 0] ;
  % ------------------------

  % --- set the local degrees of freedom corresponding to each behavior
  LocAxialdofs  = [ 1 7 ] ;
  LocTorsndofs  = [ 2 8 ] ;
  LocBendXYdofs = [ 3 6 9 12 ] ;
  LocBendXZdofs = [ 5 4 11 10 ] ;

  KL = zeros ( 2*ndofpnode, 2*ndofpnode ) ;

  Kaxial = E*A/l * [ 1 -1 ; ...
                    -1  1 ] ;
  KL( LocAxialdofs , LocAxialdofs ) = Kaxial ;

  % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

  uvector = [u1; u2] ;
  vvector = [v1; v2] ;
  thetavector = [theta1; theta2] ;

  dvector = [uvector', vvector', thetavector']' ;
  
  Bu = [-1/l 1/l] ;

  Kfd = zeros(6,6) ;
  Kfalfa = zeros(6,6) ;
  Khd = zeros(6,6) ;
  Khalfa = zeros(6,6) ;

  % Gauss-Lobatto Quadrature with 3 integration points [a (a+b)/2 b]

  npi = 3 ;
  xpi = [a (a+b)/2 b] ;
  wpi = [1/3 4/3 1/3] ;

  % Integration

  for j = 1:npi
  [Kfdj, Kfalfaj, Khdj, Khalfaj] = integrand(xpi(j),wpi(j)) ;
  Kfd = Kfd + Kfdj ;
  Kfalfa = Kfalfa + Kfalfaj ;
  Khd = Khd + Khdj ;
  Khalfa = khalfa + Khalfaj ;
  end

  function [Kfd, Kfalfa, Khd, Khalfa] = integrand(xpi)

  N = bendingInterFuns (xpi, le, 2) ;
  Bv = [N(1) N(3)] ;
  Btheta = [N(2) N(4)] ;

  Bd = [Bu 0 0; 0 Bv Btheta] ;

  Cep
  Ghat

  Kfd = Bd'*E*A*Bd ;
  Kfalfa = Bd'*E*A*Bd ;
  Khd = Bd'*E*A*Bd ;
  Khalfa = Bd'*E*A*Bd ;
  
  end

  % /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\

  % bending XY
  if     elemReleases(3) == 0 && elemReleases(4) == 0
    KbendXY = E * Iz / l^3 * kBendNoRelease ;
  elseif elemReleases(3) == 1 && elemReleases(4) == 0
    KbendXY = E * Iz / l^3 * kBendReleaseLef ;
  elseif elemReleases(3) == 0 && elemReleases(4) == 1
    KbendXY = E * Iz / l^3 * kBendReleaseRig ;
  else
    KbendXY = zeros(4,4) ;
  end
  
	% bending XZ
	RXYXZ = eye(4) ; RXYXZ(2,2) = -1; RXYXZ(4,4) = -1;
	if     elemReleases(1) == 0 && elemReleases(2) == 0
		KbendXZ = E * Iy / l^3 * RXYXZ * kBendNoRelease * RXYXZ ;
	elseif elemReleases(1) == 1 && elemReleases(2) == 0
		KbendXZ = E * Iy / l^3 * RXYXZ * kBendReleaseLef * RXYXZ ;
	elseif elemReleases(1) == 0 && elemReleases(2) == 1
		KbendXZ = E * Iy / l^3 * RXYXZ * kBendReleaseRig * RXYXZ ;
	else
		KbendXZ = zeros(4,4) ;
	end
	
  Ktorsn = G*J/l * [  1 -1  ; ...
                     -1  1  ] ;
	
	
  KL( LocBendXYdofs , LocBendXYdofs ) = KbendXY ;
  KL( LocBendXZdofs , LocBendXZdofs ) = KbendXZ ;
  KL( LocTorsndofs  , LocTorsndofs  ) = Ktorsn ;

  KGelem = R * KL * R' ;
  Finte = KGelem * Ut ;

  finteLocalCoor =  R' * Finte ;

  fs{1} = Finte  ;
  ks{1} = KGelem ;

  if density > 0
    Xe = elemCoords(:) ;
    localAxisRef = Xe(4:6) - Xe(1:3) ;
    lini = sqrt( sum( localAxisRef.^2 ) ) ;
    Me = sparse( 12, 12 ) ;

    if strcmp(massMatType, 'consistent')
    
          MeBending = density * A *  l / 420 *       [156     22*l    54     -13*l   ;...
                                                22*l    4*l^2   13*l   -3*l^2  ;...
                                                54      13*l    156    -22*l   ;...
                                                -13*l   -3*l^2  -22*l  4*l^2  ] ;

          MeAxial   = density * A * l / 6 *     [ 2  1;...
                                                  1  2];
                                                  
          Me(LocBendXYdofs,LocBendXYdofs) = MeBending ;
          Me(LocBendXZdofs,LocBendXZdofs) = RXYXZ * MeBending * RXYXZ ;
          Me(LocAxialdofs, LocAxialdofs)  = MeAxial;
          
    elseif strcmp(massMatType, 'lumped')
    
          Me (1:2:end, 1:2:end) = density * A * lini * 0.5 * eye(6) ;
          
    else
      error('the massMatType field into the elements struct must be or consistent or lumped' )
    end
    
    Fmasse = Me * Udotdotte ;

    fs{3} = Fmasse  ;
    ks{3} = Me      ;
  elseif density == 0
    fs{3} = zeros(12,1) ;
    ks{2} = zeros(12)   ;
    ks{3} = zeros(12)   ;
  end

end

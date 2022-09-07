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
 
% --------------------------------------------------------------------------------------------------

% =============================================================================
function [ fs, ks ] = linearStiffMatBeam3D(elemCoords, elemCrossSecParams, massMatType, density, hyperElasModel, hyperElasParams, Ut, Udotdotte, intBool, boolMatNonLin, matFintBool, elem)
  
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
  %~ ------------------------
  elemReleases = [0 0 0 0] ;
  %~ ------------------------

  % --- set the local degrees of freedom corresponding to each behavior
  LocAxialdofs  = [ 1 7 ] ;
  LocTorsndofs  = [ 2 8 ] ;
  LocBendXYdofs = [ 3 6 9 12 ] ;
  LocBendXZdofs = [ 5 4 11 10 ] ;

  KL = zeros ( 2*ndofpnode, 2*ndofpnode ) ;

  Kaxial = E*A/l * [ 1 -1 ; ...
                    -1  1 ] ;
  KL( LocAxialdofs , LocAxialdofs ) = Kaxial ;

  kBendReleaseRig = [ 3    3*l   -3   0 ; ...
                      3*l  3*l^2 -3*l 0 ; ...
                     -3   -3*l    3   0 ; ...
                      0    0      0   0 ] ;

  kBendReleaseLef = [  3   0 -3   3*l   ; ...
                       0   0  0   0     ; ...
                      -3   0  3  -3*l   ; ...
                       3*l 0 -3*l 3*l^2 ] ;

  % K bending in local coordinates
  kBendNoRelease = [  12     6*l    -12     6*l   ; ...
                       6*l   4*l^2   -6*l   2*l^2 ; ...
                     -12    -6*l     12    -6*l   ; ...
                       6*l   2*l^2   -6*l   4*l^2 ] ;

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
  
	% ====================================================================  
	% Material non Linearity 
	% ====================================================================
	if boolMatNonLin == 0
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
		
	elseif boolMatNonLin == 1
	
		finteKTe_int		
	
	end % endif mat non linearity
	
  Ktorsn = G*J/l * [  1 -1  ; ...
                     -1  1  ] ;
	
	
  KL( LocBendXYdofs , LocBendXYdofs ) = KbendXY ;
  KL( LocBendXZdofs , LocBendXZdofs ) = KbendXZ ;
  KL( LocTorsndofs  , LocTorsndofs  ) = Ktorsn ;

  KGelem = R * KL * R' ;
  Finte = KGelem * Ut ;
  
  %~ RXYXZ = eye(4) ; RXYXZ(2,2) = -1; RXYXZ(4,4) = -1;
	%~ KbendXZ = E * Iy / l^3 * RXYXZ * kBendNoRelease * RXYXZ 
  
  if boolMatNonLin == 1
		Finte(LocBendXZdofs) = R(LocBendXZdofs,LocBendXZdofs)*finte ;
  end

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
% ==============================================================================
%~ % ==============================================================================
%~ function [KGelem] = linearStiffTimoMatBeam3D

  %~ G = E / ( 2*(1 + nu) ) ;
  %~ ky = 5/6 ;
  %~ kz = 5/6 ; % Shear factor

  %~ Da = E * A ;
  %~ Dby = E * Iy ;
  %~ Dbz = E * Iz ;
  %~ Dsy = G * ky * A ;
  %~ Dsz = G * kz * A ;
  %~ Dt = G * J ;

  %~ % Local matrices

  %~ Ke11 = [ Da/l  0       0       0     0                     0                   ; ...
           %~ 0     Dsy/l   0       0     0                     Dsy/2               ; ...
           %~ 0     0       Dsz/l   0    -Dsz/2                 0                   ; ...
           %~ 0     0       0       Dt/l  0                     0                   ; ...
           %~ 0     0      -Dsz/2   0     (Dsz*l/4 + Dby/l)     0                   ; ...
           %~ 0     Dsy/2   0       0     0                     (Dsy*l/4 + Dbz/l)   ] ;

  %~ Ke12 = [ -Da/l   0       0        0      0                   0                 ; ...
            %~ 0     -Dsy/l   0        0      0                   Dsy/2             ; ...
            %~ 0      0      -Dsz/l    0     -Dsz/2               0                 ; ...
            %~ 0      0       0       -Dt/l   0                   0                 ; ...
            %~ 0      0       Dsz/2    0      (Dsz*l/4 + Dby/l)   0                 ; ...
            %~ 0     -Dsy/2   0        0      0                   (Dsy*l/4 - Dbz/l) ] ;

  %~ Ke22 = [ Da/l   0       0       0       0                 0                 ; ...
           %~ 0      Dsy/l   0       0       0                -Dsy/2             ; ...
           %~ 0      0       Dsz/l   0       Dsz/2             0                 ; ...
           %~ 0      0       0       Dt/l    0                 0                 ; ...
           %~ 0      0       Dsz/2   0       (Dsz*l/4 + Dby)   0                 ; ...
           %~ 0     -Dsy/2   0       0       0                 (Dsy*l/4 + Dbz/l) ] ;

  %~ Ke21 = abs(Ke12)' ;

  %~ Kelem = [ Ke11 Ke12 ; ...
            %~ Ke21 Ke22 ] ;


%~ end
% ==============================================================================

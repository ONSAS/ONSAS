<<<<<<< HEAD:src/elements/frame/linearStiffMatBeam3D.m
% Copyright 2023, Jorge M. Perez Zerpa, Joaquin Viera, Mauricio Vanzulli.
=======
% Copyright 2023, ONSAS Authors (see documentation)
>>>>>>> master:src/elements/frame/elementFrameLinear.m
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
<<<<<<< HEAD:src/elements/frame/linearStiffMatBeam3D.m
 
function [fs, ks, finteLocalCoor] = linearStiffMatBeam3D(elemCoords, elemCrossSecParams, ...
                                                         massMatType, density, modelName, modelParams, ...
                                                         Ut, Udotdotte)
  ndofpnode = 6 ;
  
  % --- material constitutive params ---
	if strcmp(modelName,'linearElastic')
    E   = modelParams(1) ;
	  nu  = modelParams(2) ;
	  G   = E/(2*(1+nu)) ;
	elseif strcmp(modelName,'isotropicHardening')

  end

=======
% 
% --------------------------------------------------------------------------------------------------

% =============================================================================
function [ fs, ks, finteLocalCoor ] = elementFrameLinear(elemCoords, elemCrossSecParams, massMatType, density, modelName, modelParams, Ut, Udotdotte)
  
  ndofpnode = 6 ;
  
  % --- material constit params ---
	E   = modelParams(1) ;
	nu  = modelParams(2) ;
	G   = E/(2*(1+nu)) ;
	
>>>>>>> master:src/elements/frame/elementFrameLinear.m
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

if strcmp(modelName,'isotropicHardening')

global ne
global secFinteVec
global elemFuncs

KbendXZ = zeros(4,4) ;
finte 	= zeros(4,1) ;

elemCrossSecParamsVec = elemCrossSecParams{2} 

RXYXZ = eye(4) ; RXYXZ(2,2) = -1; RXYXZ(4,4) = -1;

% Elem Gauss points
[xge, we] = gaussPointsAndWeights(ne) ;
pgeVec = ( l/2  * xge' + l/2 ) ;	

for j = 1:length(we)
	secFint 	= 0 ;
	secKTe 		= 0 ;
	pge 			= pgeVec(j) ; 
	
	% Bending intern functions second derivative
	B = bendingInterFuns(pge, l, 2)*RXYXZ ;
	
	if intBool == 1
		% Tangent stiffness matrix
		secKTe = quadv('secKT', -elemCrossSecParamsVec(2)/2, elemCrossSecParamsVec(2)/2, [], [], ...
														 elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), B, R(LocBendXZdofs,LocBendXZdofs), ...
														 Ut(LocBendXZdofs), hyperElasParams, hyperElasModel) ;
		KbendXZ = l/2*( B'*secKTe*B*we(j) ) + KbendXZ ;	
	end
	
	secFint = quadv('secFint', -elemCrossSecParamsVec(2)/2, elemCrossSecParamsVec(2)/2, [], [], ...
															elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), B, R(LocBendXZdofs,LocBendXZdofs), ...
															Ut(LocBendXZdofs), hyperElasParams, hyperElasModel) ;
	finte = l/2*secFint*we(j) + finte ;	
	
	% ==============================
	if matFintBool == 1 && ( elem == elemFuncs || elem == elemFuncs+1 )
		if elem == elemFuncs
			secFinteVec(1,end+1) = secFint(4) ;
		else
			secFinteVec(2, (size(secFinteVec,2)-ne)+j ) = secFint(4) ;
		end	
	end
		
		
	% ==============================
													
end % endfor we







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
    
    fs{3} = Me * Udotdotte  ;
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



function FintInt = secFint( z, ty, tz, B, R, Ut, hyperElasParams, hyperElasModel)
		
	epsk = epsVal(z, B, R, Ut) ;
	[sigma, ~] = constitutiveModel(hyperElasParams, hyperElasModel, epsk') ;
	FintInt = ty * -B' * (z .* sigma)' ;


function KTInt = secKT( z, ty, tz, B, R, Ut, hyperElasParams, hyperElasModel )
		
	epsk = epsVal(z, B, R, Ut) ;
	[~, dsigdeps] = constitutiveModel(hyperElasParams, hyperElasModel, epsk') ;
	KTInt = ty * dsigdeps * z.^2  ;

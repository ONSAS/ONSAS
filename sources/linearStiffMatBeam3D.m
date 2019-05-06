%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

% --------------------------------------------------------------------------------------------------

% ==============================================================================
function [KGelem] = linearStiffMatBeam3D(E, nu, A, Iy, Iz, J, l, elemReleases, R)

  ndofpnode = 6 ;  
  G  = E / ( 2*(1+nu) ) ;
  
  % set the locales degrees of freedom corresponding to each behavior
  LocAxialdofs  = [ 1  1+ndofpnode                        ] ;
  LocTorsndofs  = [ 2  2+ndofpnode                        ] ;
  LocBendXYdofs = [ 3  6           3+ndofpnode 6+ndofpnode] ;
  LocBendXZdofs = [ 5  4           5+ndofpnode 4+ndofpnode] ;
  
  KL = zeros ( ndofpnode*2, ndofpnode*2 ) ;
  
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
  
  Ktorsn = G*J/l * [  1 -1  ; ...
                     -1  1  ] ;
    
  KL( LocBendXYdofs , LocBendXYdofs ) = KbendXY ;
  KL( LocBendXZdofs , LocBendXZdofs ) = KbendXZ ;
  KL( LocTorsndofs  , LocTorsndofs  ) = Ktorsn ;
      
    
  KGelem = R * KL * R' ;
    
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

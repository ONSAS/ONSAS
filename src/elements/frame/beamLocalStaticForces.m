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
 
function [fl,kl, strain, stress] = beamLocalStaticForces (u, tl1, tl2, L, E, G, A, Iyy, Izz, J);

  ltx1 = tl1(1);
  lty1 = tl1(2);
  ltz1 = tl1(3);

  ltx2 = tl2(1);
  lty2 = tl2(2);
  ltz2 = tl2(3);

  % --- internal forces vector --
  fl = zeros(7,1);

  % Auxiliar epsilon x
  strain = u/L ;
  % Auxiliar sigma x
  stress = E*strain ;

  fl(1) = E*A*u/L;

  fl(2) = - G * J / L * (ltx2-ltx1) ;
  fl(5) = -fl(2);

  fl(3) = 2*E*Iyy/L * ( 2*lty1 +   lty2 ) ;
  fl(6) = 2*E*Iyy/L * (   lty1 + 2*lty2 ) ;

  fl(4) = 2*E*Izz/L * ( 2*ltz1 +   ltz2 ) ;
  fl(7) = 2*E*Izz/L * (   ltz1 + 2*ltz2 ) ;

  % stiffness matrix
  kl = zeros(7,7);

  kl(1,1) = E*A/L ;

  kl([2 5], [2 5] ) = G * J / L  * [ 1 -1 ; -1 1 ];

  kl([3 6], [3 6] ) = 2* E * Iyy / L * [ 2  1 ; 1 2 ];

  kl([4 7], [4 7] ) = 2* E * Izz / L * [ 2  1 ; 1 2 ];

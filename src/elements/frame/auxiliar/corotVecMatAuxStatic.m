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

function [nu, nu11, nu12, nu21, nu22, e1, e2, e3, r, Gaux, P, EE ] = corotVecMatAuxStatic(R0, Rr, Rg1, Rg2, l, II, O3, O1)

  % global coords
  q1g = Rg1 * R0 * [0 1 0]' ;
  q2g = Rg2 * R0 * [0 1 0]' ;
  qg  = ( q1g + q2g ) / 2     ;

  % local coords
  q  = Rr' *  qg ;
  q1 = Rr' * q1g ;

  nu   = q(1) / q(2)   ;
  nu11 = q1(1) / q(2)  ;
  nu12 = q1(2) / q(2)  ;
  nu21 = 2 * nu - nu11 ;
  nu22 = 2 - nu12      ;

  Gaux = [0   0    nu/l  nu12/2  -nu11/2  0  0  0    -nu/l  nu22/2  -nu21/2  0
          0   0    1/l     0        0     0  0  0    -1/l     0        0     0
          0  -1/l  0       0        0     0  0  1/l   0       0        0     0]';

  P = II - [Gaux'; Gaux'] ;

  % rigid base 
  e1 = Rr(:, 1) ;
  e2 = Rr(:, 2) ;
  e3 = Rr(:, 3) ;
  %  transform local to global axial force  
  r = [ -e1' O1  e1' O1 ]' ;

  EE=[ Rr O3 O3 O3 
       O3 Rr O3 O3 
       O3 O3 Rr O3 
       O3 O3 O3 Rr ] ;
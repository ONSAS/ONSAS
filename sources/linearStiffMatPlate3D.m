% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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

% ==============================================================================
function [Kelem] = linearStiffMatPlate3D(E, nu, t, Lx, Ly)
  
  D  = E * t^3 / ( 12 * (1-nu^2) ) ;
  a  = Lx /  2  ;
  b  = Ly /  2  ;
  
  Ke1 = b / ( 6 * a^3 )  * ...
    [ 6         6*a     0    -6      6*a       0    -3       3*a     0     3      3*a       0  ;
      6*a       8*a^2   0   -6*a     4*a^2     0   -3*a     2*a^2    0    3*a     4*a^2     0  ;
      0          0      0     0       0        0     0        0      0     0       0        0  ;
      -6       -6*a     0     6     -6*a       0     3      -3*a     0    -3     -3*a       0  ;
      6*a      4*a^2    0   -6*a     8*a^2     0   -3*a     4*a^2    0    3*a     2*a^2     0  ;
      0         0       0     0       0        0     0        0      0     0       0        0  ;
      -3      -3*a      0     3     -3*a       0     6     -6*a      0    -6     -6*a       0  ;
      3*a      2*a^2    0   -3*a     4*a^2     0   -6*a     8*a^2    0    6*a    4*a^2      0  ;
      0         0       0     0       0        0     0        0      0     0       0        0  ;
      3        3*a      0    -3      3*a       0    -6       6*a     0     6      6*a       0  ;
      3*a     4*a^2     0   -3*a     2*a^2     0   -6*a      4*a^2   0    6*a     8*a^2     0  ;
      0         0       0     0       0        0     0        0      0     0       0        0  ] ;
  
  Ke2 = a / (6*b^3) * ...
    [ 6       0      6*b       3    0    3*b       -3       0      3*b     -6      0     6*b   ;
      0       0       0        0    0     0         0       0       0       0      0      0    ;
      6*b     0      8*b^2    3*b   0    4*b^2    -3*b      0      2*b^2  -6*b     0     4*b^2 ;
      3       0      3*b       6    0    6*b       -6       0      6*b     -3      0     3*b   ;
      0       0       0        0    0     0         0       0       0       0      0      0    ;
      3*b     0      4*b^2    6*b   0    8*b^2    -6*b      0      4*b^2  -3*b     0     2*b^2 ;
      -3      0     -3*b      -6    0   -6*b        6       0     -6*b      3      0    -3*b   ;
      0       0       0        0    0     0         0       0       0       0      0      0    ;
      3*b     0      2*b^2    6*b   0    4*b^2    -6*b      0      8*b^2  -3*b     0     4*b^2 ;
      -6      0     -6*b      -3    0   -3*b        3       0     -3*b      6      0    -6*b   ;
      0       0       0        0    0     0         0       0       0       0      0      0    ;
      6*b     0      4*b^2    3*b   0    2*b^2    -3*b      0      4*b^2  -6*b     0     8*b^2 ] ;
  
  Ke3 = nu / (2*a*b) * ...
    [  1         a          b        -1      0         -b        1     0       0        -1      -a        0    ;
       a         0         2*a*b      0      0          0        0     0       0        -a       0        0    ;
       b        2*a*b       0        -b      0          0        0     0       0         0       0        0    ;
      -1         0         -b         1     -a          b       -1     a       0         1       0        0    ;
       0         0         0         -a      0        -2*a*b     a     0       0         0       0        0    ;
      -b         0         0          b    -2*a*b       0        0     0       0         0       0        0    ;
       1         0         0         -1      a          0        1    -a      -b        -1       0        b    ;
       0         0         0          a      0          0       -a     0      2*a*b      0       0        0    ;
       0         0         0          0      0          0       -b    2*a*b    0         b       0        0    ;
      -1        -a         0          1      0          0       -1     0       b         1       a       -b    ;
      -a         0         0          0      0          0       0      0       0         a       0      -2*a*b ; 
       0         0         0          0      0          0       b      0       0        -b     -2*a*b     0    ] ;
  
  
  Ke4 = (1-nu)/ (30*a*b) * ...
    [ 21       3*a        3*b      -21        3*a        -3*b         21       -3*a         -3*b        -21       -3*a        3*b   ;
      3*a      8*a^2       0       -3*a      -2*a^2        0          3*a       2*a^2         0         -3*a      -8*a^2       0    ;
      3*b       0         8*b^2    -3*b        0         -8*b^2       3*b        0	         2*b^2      -3*b        0        -2*b^2 ;
      -21     -3*a       -3*b	      21       -3*a         3*b        -21        3*a          3*b         21        3*a       -3*b   ;
      3*a     -2*a^2       0       -3*a       8*a^2        0          3*a      -8*a^2         0         -3*a       2*a^2       0    ;
      -3*b      0        -8*b^2     3*b        0          8*b^2      -3*b        0          -2*b^2       3*b        0         2*b^2 ;
      21       3*a        3*b      -21         3*a       -3*b         21       -3*a         -3*b        -21       -3*a        3*b   ;
      -3*a     2*a^2       0        3*a      -8*a^2        0         -3*a       8*a^2         0          3*a      -2*a^2       0    ;
      -3*b      0         2*b^2     3*b        0         -2*b^2	     -3*b        0           8*b^2       3*b        0        -8*b^2 ;
      -21     -3*a       -3*b       21       -3*a         3*b        -21        3*a          3*b         21        3*a       -3*b   ;
      -3*a    -8*a^2       0        3*a       2*a^2        0         -3*a      -2*a^2         0          3*a       8*a^2       0    ;
      3*b       0        -2*b^2    -3*b        0          2*b^2       3*b        0          -8*b^2      -3*b        0         8*b^2 ] ;
  
  
  Kelem = D * (Ke1 + Ke2 + Ke3 + Ke4 ) ;

end
% ==============================================================================

% ==============================================================================

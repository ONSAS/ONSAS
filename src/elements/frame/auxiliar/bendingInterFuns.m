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

function N = bendingInterFuns (x , l, derivdeg )

assert( iscolumn(x)   ,'x must be scalar or column' )
assert( length(l)==1                                )

switch derivdeg

case 0
  N1 =  (2*x.^3 - 3*l*x.^2 + l.^3 )     / l^3 ;
  N2 =  (x.^3 -2 * l * x.^2 + l.^2*x  ) / l^2 ;
  N3 = -( 2*x.^3 - 3*l*x.^2 )           / l^3 ;
  N4 =  (x.^3 - l*x.^2 )                / l^2 ;

case 1
  N1 =  (6*x.^2 - 6*x*l )     / l^3 ;
  N2 =  (3*x.^2 -4*l*x +l.^2) / l^2 ;
  N3 = -(6*x.^2 - 6*x*l)      / l^3 ;
  N4 =  (3*x.^2 - 2*l*x )     / l^2 ;

case 2
  N1 =  (12*x-6*l )  / l^3 ;
  N2 =  (6*x-4*l )   / l^2 ;
  N3 = -(12*x - 6*l) / l^3 ;
  N4 =  (6*x - 2*l ) / l^2 ;

case 3
  N1 =  (12)  / l^3 ;
  N2 =  (6 )  / l^2 ;
  N3 = -(12)  / l^3 ;
  N4 =  (6 )  / l^2 ;

end

N= [N1 N2 N3 N4] ;

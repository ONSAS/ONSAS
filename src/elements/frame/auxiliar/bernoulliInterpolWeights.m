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

function [ N1, N2, N3, N4, N5, N6, N7, N8 ] = bernoulliInterpolWeights(x, l0) 

  % Shape functions of Euler Bernoulli element to interpolate displacements and velocites for the cross section:
  % linear
  N1 = 1 -x / l0                           ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N2 = x / l0                              ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  % cubic
  N3 = x * ( 1 - x / l0 )^2                ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N4 = - ( 1 - x / l0 ) * ( x^2 ) / l0     ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N5 = ( 1 - 3 * x / l0) * ( 1 - x / l0 )  ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N6 = ( 3 * x / l0 - 2 ) * ( x / l0 )	   ; % Eq.(37)  T-N Le J.-M. Battini et al 2014
  N7 = N3 + N4 		                         ;
  N8 = N5 + N6 -1		                       ; 

end
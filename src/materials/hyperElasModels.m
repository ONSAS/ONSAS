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
 
% function for hyperelastic models input.

function [sigma, dsigdeps ] = hyperElasModels (epsilon, paramsmodel )

nummodel    = paramsmodel(1)     ;
paramsmodel = paramsmodel(2:end) ;

switch nummodel
case 1
  % Linear material model (saint-venant-kirchhoff )
  E = paramsmodel(1) ;
  sigma = E * epsilon ;
  dsigdeps = E ;
case 2
  % Linear material model (saint-venant-kirchhoff )
  E = paramsmodel(1) ;
  sigma = E * epsilon ;
  dsigdeps = E ;
case 3
  % Linear material model (saint-venant-kirchhoff )
  E = paramsmodel(1) ;
  sigma = E * epsilon ;
  dsigdeps = E ;

case 7
  % Bi-modulus material model: param1: Tension young modulus, param2: Compression modulus
  ET = paramsmodel(1) ;
  EC = paramsmodel(2) ;

  if epsilon >=(-1e-15)
    sigma    = ET * epsilon ;
    dsigdeps = ET ;
  else
    sigma    = EC * epsilon ;
    dsigdeps = EC ;
  end

case 8
  % Bi-modulus material model with pre-strain: param1: Tension young modulus,
  % param2: Compression modulus, param3: epszero (pre-strain)
  ET = paramsmodel(1) ;
  EC = paramsmodel(2) ;
  epszero = paramsmodel(3) ;

  if ( epsilon - epszero ) >=(-1e-15)
    sigma    = ET * ( epsilon - epszero ) ;
    dsigdeps = ET ;
  else
    sigma    = EC * ( epsilon - epszero ) ;
    dsigdeps = EC ;
  end

case 9
  % Example 5.11 from Reddy nonlinear sqrt stress-strain constitutive relation
  E = paramsmodel(1) ;
  sigma    = sign(epsilon) * E * sqrt( abs( epsilon) )  ;
  if epsilon == 0
    dsigdeps = E^2 ;
  else
    dsigdeps = 0.5 * E / sqrt( abs(epsilon) ) ;
  end

end

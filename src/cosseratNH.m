% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro  
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


function [S, ConsMat] = cosseratNH( consParams, Egreen, consMatFlag)

young   = consParams(1) ;
nu      = consParams(2) ;

lambda  = young * nu / ( (1 + nu) * (1 - 2*nu) ) ;
shear   = young      / ( 2 * (1 + nu) )          ;

C       = 2*Egreen + eye(3);  % Egreen = 1/2 (C - I)
invC    = inv(C);
detC    = det(C); % TODO use analyDet ?
J       = sqrt(detC);
S       = shear * (eye(3) - invC) + lambda * log(J) * invC;

if consMatFlag == 0 % only stress computed
  ConsMat = [] ;

elseif consMatFlag == 1 % complex-step computation expression

  ConsMat = zeros(6,6);
  ConsMat = complexStepConsMat( 'cosseratNH', consParams, Egreen ) ;

else
  error("the analytical expression for the Neo-Hookean constitutive law is not available")

end

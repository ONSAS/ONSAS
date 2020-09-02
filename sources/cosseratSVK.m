% Copyright (C) 2020, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera, 
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


function [S, ConsMat] = cosseratSVK( consParams, Egreen, consMatFlag )

young   = consParams(1) ;
nu      = consParams(2) ;

lambda  = young * nu / ( (1 + nu) * (1 - 2*nu) ) ;
shear   = young      / ( 2 * (1 + nu) )          ;

S       = lambda * trace(Egreen) * eye(3)  +  2 * shear * Egreen ;

if consMatFlag == 0 % only stress computed
  ConsMat = [] ;
  
elseif consMatFlag == 1 % complex-step computation expression

  ConsMat = zeros(6,6);
  ConsMat = complexStepConsMat( 'cosseratSVK', consParams, Egreen ) ;

elseif consMatFlag == 2  % analytical expression

  ConsMat = zeros(6,6);
  ConsMat (1,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ 1-nu , nu   , nu   ] ; 
  ConsMat (2,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ nu   , 1-nu , nu   ] ;
  ConsMat (3,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ nu   , nu   , 1-nu ] ;
  ConsMat (4,4  ) = shear ;
  ConsMat (5,5  ) = shear ;
  ConsMat (6,6  ) = shear ;

end

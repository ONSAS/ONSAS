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
 

function [S, ConsMat] = cosseratSVK( consParams, Egreen, consMatFlag )

lambda   = consParams(1) ;
shear    = consParams(2) ;

S       = lambda * trace(Egreen) * eye(3)  +  2 * shear * Egreen ;

if consMatFlag == 0 % only stress computed
  ConsMat = [] ;

elseif consMatFlag == 1 % analytical expression
  ConsMat              = zeros(6,6);
  ConsMat ( 1:3, 1:3 ) = 2*shear * eye(3) + lambda * ones(3) ;
  ConsMat ( 4:6, 4:6 ) =   shear * eye(3)                   ;

elseif consMatFlag == 2  % complex-step computation expression
  ConsMat = zeros(6,6);
  ConsMat = complexStepConsMat( 'cosseratSVK', consParams, Egreen ) ;

end

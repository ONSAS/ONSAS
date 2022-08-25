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
 

function [S, ConsMat] = cosseratNHC( consParams, Egreen, consMatFlag)

shear    = consParams(1) ;
bulk     = consParams(2) ;

C       = 2*Egreen + eye(3);  % Egreen = 1/2 (C - I)
invC    = inv(C);
detC    = det(C); % TODO use analyDet ?
J       = sqrt(detC);
S       = shear * ( eye(3) - invC ) + bulk * ( J * (J-1)* invC) ;

if consMatFlag == 0 % only stress computed
  ConsMat = [] ;

elseif consMatFlag == 1 % analytic expression
  error("the analytical expression for the Neo-Hookean constitutive law is not available")

elseif consMatFlag == 2 % complex-step computation expression

  ConsMat = zeros(6,6);
  ConsMat = complexStepConsMat( 'cosseratNHC', consParams, Egreen ) ;
end

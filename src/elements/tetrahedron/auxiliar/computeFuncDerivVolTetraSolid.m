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
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

function [ funder, jacobianmat, vol, tetCoordMat ] = computeFuncDerivVolTetraSolid( elemCoords )

tetCoordMat        = reshape( elemCoords', 3, 4 ) ;

eleCoordMat = tetCoordMat               ;

xi = 0.25 ;  wi = 1/6  ;

% matriz de derivadas de fun forma respecto a coordenadas isoparametricas
deriv = shapeFuns( xi, xi , xi , 1 ) ;
 
% jacobiano que relaciona coordenadas materiales con isoparametricas
jacobianmat = eleCoordMat * deriv'  ;

vol = analyDet( jacobianmat ) / 6.0 ;

if vol<0,  vol, error('Element with negative volume, check connectivity.'), end

funder = inv(jacobianmat)' * deriv ;

% Copyright 2023, ONSAS Authors (see documentation)
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
%
function [ funder, jacobianmat, vol, tetCoordMat ] = computeFuncDerivVolTetraSolid( elemCoords )

tetCoordMat = reshape( elemCoords', 3, 4 ) ;

% Compute the derivatives of the shape functions in a matrix at x = xi
% point to evaluate 
xi = 0.25 ;  wi = 1/6  ;
deriv = shapeFuns( xi, xi , xi , 1 ) ;
 
% Jacobian linking material and isoparametric coordinates
jacobianmat = tetCoordMat * deriv'  ;

% tetrahedron volume
vol = analyDet( jacobianmat ) / 6.0 ;

if vol<0,  vol, error('Element with negative volume, check connectivity.'), end

funder = inv(jacobianmat)' * deriv ;

%~ Copyright (C) 2019, Jorge M. Pérez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquín Viera, Mauricio Vanzulli  

%~ This file is part of ONSAS.

%~ ONSAS is free software: you can redistribute it and/or modify
%~ it under the terms of the GNU General Public License as published by
%~ the Free Software Foundation, either version 3 of the License, or
%~ (at your option) any later version.

%~ ONSAS is distributed in the hope that it will be useful,
%~ but WITHOUT ANY WARRANTY; without even the implied warranty of
%~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%~ GNU General Public License for more details.

%~ You should have received a copy of the GNU General Public License
%~ along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.

% function that calculate the visual size for forces

function [maxNorm2F, visualloadfactor] =visualLoadFac(strucsize, variableFext, constantFext, nnodes)
	
	ndofpnode = 6 ;
	if ( norm( variableFext ) + norm( constantFext ) ) == 0
		visualloadfactor = 0 ;
    maxNorm2F = 1 ;
	else
		maxNorm2F = max( [ sqrt( max( ( [ reshape(variableFext, nnodes, ndofpnode)' ...
                                    reshape(constantFext, nnodes, ndofpnode)' ].^2 )' ) ) ] ) ;
		visualloadfactor = 0.1 * strucsize / maxNorm2F ;
	end

end

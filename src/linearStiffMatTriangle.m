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

function [Kelem] = linearStiffMatTriangle(E, nu, Nodes, t, state)

	xcoord = Nodes(:,1) ;
	ycoord = Nodes(:,2) ;

	A = det( [ ones(1,3) ; xcoord' ; ycoord' ] ) / 2; 

	B = 1 / (2*A) * [ y(2)-y(3)  				 0  y(3)-y(1)  				 0  y(1)-y(2)  				 0  ; ...
														0  x(3)-x(2) 				  0  x(1)-x(3)  				0  x(2)-x(1)  ; ...
										x(3)-x(2)  y(2)-y(3)  x(1)-x(3)  y(3)-y(1)  x(2)-x(1)  y(1)-y(2)  ] ;
										
	if state == 1 % Plane stress 

		D = E / (1-nu^2) * [ 	1    nu  				  0   ;
												 nu  		1   				0   ;
													0   	0   (1-nu)/2 ] ;

	elseif state == 2 % Plane strain									
			
		D = E * (1-nu) / ( (1+nu)*(1-2*nu) ) = [ 				 1 	nu/(1-nu) 											0 ; ... 
																						 nu/(1-nu)	 			  1 										 	0 ; ...
																										 0 				  0 		(1-2*nu)/(2*(1-nu)) ] ;
										
	end									
										
	Kelem = B' * D * B * A * t ;

end

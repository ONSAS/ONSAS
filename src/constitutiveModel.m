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

function [sigma, dsigdeps] = constitutiveModel(hyperElasParams, hyperElasModel, epsk)

	E = hyperElasParams(1) ;
			
	% Linear elastic 
	if strcmp(hyperElasModel, 'linearElastic') 
		dsigdeps = E ;
		sigma = E * epsk ;
	% Bi Linear	
	elseif strcmp(hyperElasModel, 'biLinear')
		sigmaY = hyperElasParams(3) ;
		sigma_tr = abs( E*epsk ) ;
		if sigma_tr >= sigmaY
			K = hyperElasParams(4) ;
			epsY = sigmaY/E ;
			%~ sigma = sigmaY*sign(epsk) + K * ( epsk - epsY*sign(epsk) ) ; 
			dsigdeps = E*K / (E+K) ; 
			sigma = sigmaY*sign(epsk) + dsigdeps * ( epsk - epsY*sign(epsk) ) ; 
		else
			dsigdeps = E ;
			sigma = sigma_tr * sign(epsk) ;
		end
	elseif strcmp(hyperElasModel, 'userFunc')
		[sigma, dsigdeps] = userConsModel(hyperElasParams, epsk) ;
	end
	
end

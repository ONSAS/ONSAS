% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, J. Bruno Bazzano,
% Joaquin Viera, Marcelo Forets, Jean-Marc Battini. 
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

function [sigma, dsigdeps] = constitutiveModel(hyperElasParams, hyperElasModel, epsk, matFintBool, elem, elemAux, xgej,maxxge)

	% user function
	global userFuncBool

	E = hyperElasParams(1) ;
			
	% Linear elastic 
	if strcmp(hyperElasModel, 'linearElastic') 
		sigma = E * epsk ;
		dsigdeps = E ;
	% Elastoplastic perfect
	elseif strcmp(hyperElasModel, 'elastoPlasticPerfect') 
		sigmaY = hyperElasParams(3) ;
		sigma_tr = abs( E*epsk ) ;
		if sigma_tr >= sigmaY
			sigma = sigmaY * sign(epsk) ;
			dsigdeps = 0 ;
		else
			sigma = sigma_tr * sign(epsk) ;
			dsigdeps = E ;
		end	
	% Linear hardening	
	elseif strcmp(hyperElasModel, 'linearHardening')
		sigmaY = hyperElasParams(3) ;
		sigma_tr = abs( E*epsk ) ;
		if sigma_tr >= sigmaY
			K = hyperElasParams(4) ;
			epsY = sigmaY/E ;
			sigma = sigmaY*sign(epsk) + K * ( epsk - epsY*sign(epsk) ) ;
			dsigdeps = E*K / (E+K) ; 
		else
			sigma = sigma_tr * sign(epsk) ;
			dsigdeps = E ;
		end
	elseif strcmp(hyperElasModel, 'userFunc')
		[sigma, dsigdeps] = userConsModel(hyperElasParams, epsk, matFintBool, elem, elemAux, xgej,maxxge) ;
	end
	
end

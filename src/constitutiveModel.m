function [sigma, dsigdeps] = constitutiveModel(hyperElasParams, hyperElasModel, epsk)

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
	elseif ~isempty(userFuncBool)
		[sigma, dsigdeps] = userConsModel(hyperElasParams, epsk) ;
	end
	
	
end

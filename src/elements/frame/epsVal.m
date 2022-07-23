function epsk = epsVal(z, B, R, Ut)

	%~ epsk = -z'.* B * R * Ut ;


	eps_sup = -0.00251990926724482 ;
	tz = 0.3 ;
	
	epsk = -2*eps_sup / tz * z' + eps_sup ;
	
	
end

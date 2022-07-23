

function FintInt = secFint( z, ty, tz, B, R, Ut, hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, xge, maxXge )
		

	% Strain
	epsk = epsVal(z, B, R, Ut) ;
	% Section width
	%~ t = secWidth(z', 2, ty, tz ) ;
	t = ty ;
	% stress & slope
	[sigma, dsigdeps] = constitutiveModel(hyperElasParams, hyperElasModel, epsk', matFintBool, elem, elemAux, xge, maxXge ) ;
	
	FintInt = t * -B' * (z .* sigma)' ;

end

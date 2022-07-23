function KTInt = secKT( z, ty, tz, B, R, Ut, hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, xge, maxXge )
		
	% Strain
	epsk = epsVal(z, B, R, Ut) ;
	% Section width
	%~ t = secWidth(z', 2, ty, tz ) ;
	t = ty ;
	% stress & slope
	[~, dsigdeps] = constitutiveModel(hyperElasParams, hyperElasModel, epsk', matFintBool, elem, elemAux, xge, maxXge ) ;
	
	KTInt = t * dsigdeps * z.^2  ;

end

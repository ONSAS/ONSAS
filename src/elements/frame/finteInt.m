

function finte = finteInt( x, ty, tz, R, Ut, hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, xge, maxXge, l )
		
		
	B = bendingInterFuns(x, l, 2) ;
	
	aux = quadv('secFint', -tz/2, tz/2, [], [], ty, tz, B, R, Ut, hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, xge, maxXge ) ;
		
	finte = aux ;


end

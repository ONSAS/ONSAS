

function KTe = KTint( x, ty, tz, R, Ut, hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, l )
	
	
	B = bendingInterFuns(x, l, 2) ;

	secKTe = quadv('secKT', -tz/2, tz/2, [], [], ty, tz, B, R, Ut, hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, [], []) ;	
	
	KTe = B' * secKTe * B ;

end

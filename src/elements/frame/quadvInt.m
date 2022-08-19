% --- geometry params ---
		elemCrossSecParamsVec = elemCrossSecParams{2} ;
		
		global sigmaVec
		global elemAux
		
		KTe = zeros(4,4) ;
		finte = zeros(4,1) ;
		
		if intBool == 1
			% Tangent stiffness matrix	
			KTe = quadv('KTint', 0, l, [], [], elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), R(LocBendXZdofs,LocBendXZdofs), Ut(LocBendXZdofs), hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, l ) ;
		end
	
		% Internal force
		finte = quadv('finteInt', 0, l, [], [], elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), R(LocBendXZdofs,LocBendXZdofs), Ut(LocBendXZdofs), hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, [], [], l) ;
	
	
	KbendXZ = KTe ;


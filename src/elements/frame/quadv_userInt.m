% --- geometry params ---

		global ne
		%~ global ns
		global sigmaVec
		global elemAux
		
		KTe = zeros(4,4) ;
		finte = zeros(4,1) ;
		
		% Elem Gauss points
		[xge, we] = GaussPointsAndWeights (ne) ;
		%~ [xgs, ws] = GaussPointsAndWeights(ns) ; 

		a = 0 ;
		b = l ;
		p1 = (b-a)/2 ;
		p2 = (b+a)/2 ;
		
		%~ as = -elemCrossSecParamsVec(2)/2 ;
		%~ bs = elemCrossSecParamsVec(2)/2 ;
		%~ ps1 = (bs-as)/2 ;
		%~ ps2 = (bs+as)/2 ;
		
		pgeVec = (p1  * xge' + p2 ) ;
		%~ pgsVec = (ps1 * xgs' + ps2) ;
		
		
		
		
		for j = 1:length(we)
			secFinte = 0 ;
			secKTe = 0 ;
			
			pge = pgeVec(j) ;
			
			%%%%%%%%%%%%%%%%
			%~ n = 0 ;
			sigmaAuxVec = [] ;
			%%%%%%%%%%%%%%%%
			
			% Bending intern functions second derivative

			B = bendingInterFuns(pge, l, 2) ;
			%~ epskVec = -pgsVec*B*R(LocBendXZdofs,LocBendXZdofs)*Ut(LocBendXZdofs) ;
			
			% Section thk
			%~ tVec = secWidth(pgsVec, 2, elemCrossSecParamsVec(1), elemCrossSecParamsVec(2)) ;
			
			%~ for m = 1:length(ws)
     
				%~ pgs = pgsVec(m) ;
				%~ % Elem strain
				%~ epsk = epskVec(m) ;
				
				%~ % Elem stress
				%~ [sigma, dsigdeps] = constitutiveModel(hyperElasParams, hyperElasModel, epsk, matFintBool, elem, elemAux, xge(j), max(xge) ) ;      
				%~ t = tVec(m) ;

				%~ % Integration in section
				%~ % --------------------------------------------------------------------
				
				%~ % to compute finte
				%~ secFinte = ps1 * ( t * (-B') * pgs * sigma * ws(m) ) + secFinte ;
				
				%~ if intBool == 1
					%~ % to compute KT
					%~ secKTe = ps1 * ( t * dsigdeps * pgs^2 * ws(m) ) + secKTe ;
				%~ end
			
			%~ %%%%%%%%%%%%%%%%
			%~ if matFintBool == 1
				%~ if elem == elemAux	
					%~ if xge(j) == max(xge)
						%~ sigmaAuxVec = [sigmaAuxVec, sigma] ;
					%~ end
				%~ end
			%~ end
			%~ %%%%%%%%%%%%%%%%
		
			%~ end % endif ws
			
			
			%~ secFinte = quadv('secFint', -elemCrossSecParamsVec(2)/2, elemCrossSecParamsVec(2)/2, [], [], elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), B, R(LocBendXZdofs,LocBendXZdofs), Ut(LocBendXZdofs), hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, xge(j), max(xge)) ;
			
			if intBool == 1
				% Tangent stiffness matrix
				secKTe = quadv('secKT', -elemCrossSecParamsVec(2)/2, elemCrossSecParamsVec(2)/2, [], [], elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), B, R(LocBendXZdofs,LocBendXZdofs), Ut(LocBendXZdofs), hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, xge(j), max(xge)) ;
				KTe = p1*( B'*secKTe*B*we(j) ) + KTe ;	
			end
		
		
			% Internal force
			%~ finte = p1 * we(j) * secFinte + finte ;		
		
		%%%%%%%%%%%%%%%%
		%~ if matFintBool == 1
			%~ if elem == elemAux	
				%~ if xge(j) == max(xge)
					%~ sigmaVec{end+1} = sigmaAuxVec ;
				%~ end
			%~ end
		%~ end
		%%%%%%%%%%%%%%%%

		end % endif we
		
		finte = quadv('finteInt', 0, l, [], [], elemCrossSecParamsVec(1), elemCrossSecParamsVec(2), R(LocBendXZdofs,LocBendXZdofs), Ut(LocBendXZdofs), hyperElasParams, hyperElasModel, matFintBool, elem, elemAux, [], [], l) ;
	
		
		KbendXZ = KTe ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 			Matriz Kg del elemento	      	 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Kg_e] = elementMatrizKglobal (Dte,Ddote,Ddotdote,params,Jrho,xelem,BooleanComplexModes,NewmarkParams,numericalMethodParams)
	
	 beta   = NewmarkParams(1);
     alpha  = NewmarkParams(2);
	 DeltaT = numericalMethodParams(4);

	%loads Bt
	[Bt_e] = elementMatrizBt (Dte);
	
	%loads Mass and Gyroscopic Matrix
	if BooleanComplexModes ==1
		[MassMatrix_e,GyroMatrix_e] = MatrizMasaCompleja (xelem, Dte, Ddote, Ddotdote, params,Jrho );
	   else
		[Fine_e,MassMatrix_e,GyroMatrix_e] = elementFuerzaInercial(xelem, Dte, Ddote, Ddotdote, params,Jrho );
	end
	
	%Loads Ktang
	[ ~, KT_e] = elementBeam3DInternLoads(xelem, Dte, params ) ;
	
	%Compute Kg
	%~ MassMatrix_e
	%~ GyroMatrix_e
	Kg_e = + KT_e + 1/(beta*DeltaT^2) * MassMatrix_e * Bt_e + alpha/(beta*DeltaT) * GyroMatrix_e * Bt_e; %Ec 103 y 102
   
    %~ KG_din= + 1/(beta*DeltaT^2) * MassMatrix_e * Bt_e + alpha/(beta*DeltaT) * GyroMatrix_e * Bt_e
end

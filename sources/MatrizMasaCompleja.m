%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% 	Matrices de masa con paso complejo   	%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MassMatrixComp,GyroMatrixComp] = MatrizMasaCompleja (xelem, Ue, Udote, Udotdote, params,Jrho )
	
	%-------------- Increment ----------------------%
	h  = 1e-10 ;
	
	%-------------- Mass Matrix Comp --------------%
	MassMatrixComp = zeros (12,12);
	VecIncrement = zeros (12,1)   ;
	
	for aux = 1:12
		VecIncrement(aux) = i*h;
		[FineMatrix] = elementFuerzaInercial (xelem, Ue, Udote, Udotdote+VecIncrement, params,Jrho);
		FineMatrix = 1/h * imag(FineMatrix);
		MassMatrixComp (:,aux) = FineMatrix;
		VecIncrement = zeros (12,1);
	end


	%-------------- Gyro Matrix Comp --------------%
    GyroMatrixComp = zeros (12,12);
	VecIncrement = zeros (12,1)   ;
	for aux = 1:12
		VecIncrement(aux) = i*h;
		[FineMatrix] = elementFuerzaInercial (xelem, Ue, Udote+VecIncrement, Udotdote, params,Jrho);
		FineMatrix = 1/h * imag(FineMatrix);
		GyroMatrixComp (:,aux) = FineMatrix;
		VecIncrement = zeros (12,1);
	end

end


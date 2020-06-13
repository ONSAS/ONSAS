function B = Cambio_Base(A)

%---------------- Cambio de base matrix  -------------------


Pch = sparse (6,6);
Pch (1,1) = 1;
Pch (3,2) = 1;
Pch (5,3) = 1;
Pch (2,4) = 1;
Pch (4,5) = 1;
Pch (6,6) = 1;

P = [Pch sparse(6,6);sparse(6,6) Pch];
	
    if size(A,2)>1
		B = P'*A*P;
	 else
		B = P'*A;
	end
end

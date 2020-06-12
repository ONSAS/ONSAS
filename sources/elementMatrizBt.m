%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% 			Matriz Bt del elemento	      	 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Bt_e] = elementMatrizBt (Dte)

Bt 		 = zeros(size(Dte,1),size(Dte,1));

%Angulos en el paso k	
Wtk = Dte(2:2:end);
%Matriz Bt
for aux= 1:6:size(Dte,1)
	Bt(aux:aux+2,aux:aux+2) = eye(3,3);
	
	%calculo de la Ts
	indexANG1 = aux/2 + 1/2 ;
	indexANG2	 = indexANG1 +2;
	Wti = Wtk(indexANG1:indexANG2);
	inversaTs=invTs(Wti);
	Bt(aux+3:aux+5,aux+3:aux+5) = inversaTs'; 
end

Bt_e = Cambio_Base(Bt);

 
	function [t]=logar(R); %Funcion de matriz logarítmica
	
		u=[R(3,2)-R(2,3)
		   R(1,3)-R(3,1)
		   R(2,1)-R(1,2)];
		
		nu=norm(u);
		
		if nu==0
		  t=[0;0;0];
		  else
		  t=asin(nu/2)/nu*u;
		end
	end

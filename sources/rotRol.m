function Ro =rotRol(x);  %construye los vectores E1 t E2 y E3 y arbitrariamebnte elige y z local
		r1 = x / norm(x) ;%versor x
		
		q = [0;0;1];
		
		r2 = cross (q,r1); %Vector perpendicular a la barra en su config inde y e3
		
		if norm(r2)~=0
		  r2=r2/norm(r2);
		  [r3]=cross(r1,r2);
		else
		  q=[0;1;0];
		  [r3]=cross(r1,q);
		  r3=r3/norm(r3);
		  [r2]=cross(r3,r1);
		end
		Ro=[r1 r2 r3];
end

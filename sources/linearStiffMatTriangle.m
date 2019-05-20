function [Kelem] = linearStiffMatTriangle(E, nu, Nodes, t, state)

	xcoord = Nodes(:,1) ;
	ycoord = Nodes(:,2) ;

	A = det( [ ones(1,3) ; xcoord' ; ycoord' ] ) / 2; 

	B = 1 / (2*A) * [ y(2)-y(3)  				 0  y(3)-y(1)  				 0  y(1)-y(2)  				 0  ; ...
														0  x(3)-x(2) 				  0  x(1)-x(3)  				0  x(2)-x(1)  ; ...
										x(3)-x(2)  y(2)-y(3)  x(1)-x(3)  y(3)-y(1)  x(2)-x(1)  y(1)-y(2)  ] ;
										
	if state == 1 % Plane stress 

		D = E / (1-nu^2) * [ 	1    nu  				  0   ;
												 nu  		1   				0   ;
													0   	0   (1-nu)/2 ] ;

	elseif state == 2 % Plane strain									
			
		D = E * (1-nu) / ( (1+nu)*(1-2*nu) ) = [ 				 1 	nu/(1-nu) 											0 ; ... 
																						 nu/(1-nu)	 			  1 										 	0 ; ...
																										 0 				  0 		(1-2*nu)/(2*(1-nu)) ] ;
										
	end									
										
	Kelem = B' * D * B * A * t ;

end

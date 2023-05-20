% Copyright 2022, Jorge M. Perez Zerpa, Mauricio Vanzulli, Alexandre Villié,
% Joaquin Viera, J. Bruno Bazzano, Marcelo Forets, Jean-Marc Battini.
%
% This file is part of ONSAS.
%
% ONSAS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% ONSAS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with ONSAS.  If not, see <https://www.gnu.org/licenses/>.
 
%md Function for computation of nodal forces and tangent stiffness matrix for 2D 3 nodes triangle element with linearElastic and linear isotropic hardening with von Mises flow rule models. The dofs are assumed to be on x-y.
%md

function [ fs, ks, stress, strain, acum_plas_strain ] = elementTriangSolid( ...
  elemCoords, elemDisps, hyperElasModel, elemConstitutiveParams, paramOut, t, planeStateFlag, dotdotdispsElem, density, previous_state )

	x = elemCoords(1:3:end)' ;  y = elemCoords(2:3:end)' ;

	A = 0.5 * det( [ ones(1,3) ; x' ; y' ] ) ; % element area

	assert( A>=0, 'Element with negative area, check connectivity.')

	B = 1 / (2*A) * [ y(2)-y(3)  0         0   y(3)-y(1)  0         0   y(1)-y(2)  0         0 ;
										0          x(3)-x(2) 0   0          x(1)-x(3) 0   0          x(2)-x(1) 0 ;
										x(3)-x(2)  y(2)-y(3) 0   x(1)-x(3)  y(3)-y(1) 0   x(2)-x(1)  y(1)-y(2) 0 ] ;

	% Compute strain 
	strain = B * elemDisps ;
	
	% Plastic strain
	acum_plas_strain =  previous_state{3} ;

	%md### compute internal loads and stiffness matrix

	E  = elemConstitutiveParams(2) ;
	nu = elemConstitutiveParams(3) ;
	
	if planeStateFlag == 1
		De = E / (1-nu^2) * [ 1   nu  0           ; ...
													 nu  1   0           ; ...
													 0   0   (1-nu )/2 ] ;

	elseif planeStateFlag == 2
		De = E * (1-nu) / ( (1+nu)*(1-2*nu) ) * ...
											 [ 1          nu/(1-nu)  0                   ; ...
												 nu/(1-nu)  1          0                   ; ...
												 0          0          (1-2*nu)/(2*(1-nu)) ] ;
	end

	if strcmp( hyperElasModel, 'linearElastic' )

		dstressdeps = De ;
		stress = dstressdeps * strain ;
		stress_n           = previous_state{1,:}(1:3)  ;
		
	elseif strcmp( hyperElasModel, 'isotropicHardening') 
		
		G = E/(2*(1+nu)) 	; % Shear modulus
		K = E/(3*(1-2*nu)) 	; % Bulk modulus
		
		% Previous vals
		stress_n           = previous_state{1,:}(1:3)  ;
		strain_n           = previous_state{2,:}(1:3)  ;
		acum_plas_strain_n =  previous_state{3} 	   ;
		
		% Isotropic Hardening variables
		H       = elemConstitutiveParams(4) ;
		sigma_Y_0   = elemConstitutiveParams(5) ;
		
		% Trial 
		stress_tr = stress_n + (De * (strain - strain_n'))' ; % Eq 7.81 
		sigma_Y_tr = sigma_Y_0 + H*acum_plas_strain_n ; % Eq 7.83
		
		% Hydrostatic stress
		p = 1/3*sum(stress_tr(1:2)) ;
		
		% Deviatoric stress
		s_tr = stress_tr - p*[1,1,0] ;
		
		% J2 invariant
		norm_str = sqrt(s_tr(1)*s_tr(1) + s_tr(2)*s_tr(2) + 2*s_tr(3)*s_tr(3) ) ;
		J2_tr = 1/2 * norm_str^2 ; % Invariant J2
		
		% Yielding function
		q_tr = sqrt(3*J2_tr) 		;
		phi_tr = q_tr - sigma_Y_tr 	; % Trial
		
		% Trial strain
		strain_e_tr = De\stress_n' + (strain-strain_n') ; % Eq 7.92 
		
		% Volumetric and deviatoric split 
		strain_e_v_tr = 1/3*sum(strain_e_tr(1:2))*[1;1;0] 	; % volumetric component
		strain_e_d_tr = strain_e_tr - strain_e_v_tr			; % deviatoric component
		strain_e_d_tr(3) = strain_e_d_tr(3)/2 				;
		
		% Volumetric and deviatoric stress
		st_hyd = sum(strain_e_v_tr) * K * [1;1;0] 	;
		st_dev = 2*G*strain_e_d_tr 					;
		
		if phi_tr <= 0 % elastic behavior

			stress           	= stress_tr' 			;
			dstressdeps       = De 						;
			acum_plas_strain 	= acum_plas_strain_n 	; 
			
		else % elasto-plastic behavior
			
			delta_gamma = phi_tr / ( 3*G + H ) ; % Plastic multiplier			
			Nhat = strain_e_d_tr / norm(strain_e_d_tr) ;
			
			Id = eye(3)*2/3 ;
			Id(1,2) = -1/3 ;
			Id(2,1) = -1/3 ;
			Id(3,3) =  0.5 ;

			IxI = eye(3) ;
			IxI(3,3) = 0 ;
			IxI(1,2) = 1 ;
			IxI(2,1) = 1 ;
			
			stress           = ( De - delta_gamma*6*G^2 / q_tr * Id ) * strain_e_tr ; % eq 7.93
			dstressdeps      = 2*G*( 1 - delta_gamma*3*G/q_tr ) * Id + 6*G^2*( delta_gamma/q_tr - 1/(3*G+H) ) * (s_tr'*s_tr) / (norm_str^2) + K * IxI ; % Eq 7.120
			acum_plas_strain = acum_plas_strain_n + delta_gamma ; 

		end
		
	else	
		error("not implemented yet.")
	end
	
	% Compute Tangent Stiffness Matrix
	KTe = B' * dstressdeps * B * A * t ;
	% Internal force
	Finte = B' * stress * A * t ;

	%md### compute inertial loads and mass matrix
	Mmase = density*A/3.0 * speye(9,9);
	Fmase = Mmase * dotdotdispsElem ;

	fs = { Finte, [], Fmase } ;
	ks = { KTe,   [], Mmase } ;

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
 
%md Function for computation of nodal forces and tangent stiffness matrix for 2D 3 nodes triangle element with linearElastic behavior. The dofs are assumed to be on x-y.
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
		
		stress_n           = previous_state{1,:}  ;
		%~ disp("stressSize")
		%~ size(stress_n)
		
	elseif strcmp( hyperElasModel, 'isotropicHardening') 
		% Previous vals
		stress_n           = previous_state{1,:}  ;
		strain_n           = previous_state{2,:}  ;
		acum_plas_strain_n =  previous_state{3} ;
		
		% Isotropic Hardening variables
		%~ H       = elemConstitutiveParams(4) ;
		H = 0 ;
		sigma_Y_0   = elemConstitutiveParams(5) ;
		
		% Trial 
		stress_tr = stress_n + (De * (strain - strain_n'))' ; % Eq 7.81 
		sigma_Y_tr = sigma_Y_0 + H*acum_plas_strain_n ; % Eq 7.83
		
		%~ s_tr = stress_tr - 1/3*trace(stress_tr)*eye(3) ; % deviator tensor
		s_tr = stress_tr - 1/2*sum(stress_tr(1:2))*[1,1,0] ; % deviator tensor as a vector
		
		J2_tr = 1/2 * norm(s_tr)^2 ; % Invariant J2
		
		phi_tr = sqrt(3*J2_tr) - sigma_Y_tr ; % Trial
		
		
		%%%%
		G = E/(2*(1+nu)) ; % Shear modulus
		K = E/(3*(1-2*nu)) + G/3 ; % Bulk modulus
		
		strain_e_tr = De\stress_tr' ;
		%~ strain_e_tr2 = De\stress'  
		strain;
		
		%~ De\ stress_tr' 
		%~ De\ De*strain 
		%~ strain
		
		strain_e_v_tr = 1/2*sum(strain_e_tr(1:2))*[1;1;0]; % volumetric component
		strain_e_d_tr = strain_e_tr - strain_e_v_tr; % deviatoric component
		
		
		st_hyd = sum(strain_e_v_tr) * K * [1;1;0];
		strain_e_d_tr(3) = strain_e_d_tr(3)/2 ;
		
		st_dev = 2*G*strain_e_d_tr;
		
		stress_tr;
		st_hyd+st_dev;
		
		s_tr2 = 2*G*strain_e_d_tr ;
		
		s_tr;
		s_tr2;
		
		%~ phi_tr2 = sqrt(3/2) * 2 * G * norm ( strain_e_d_tr ) - sigma_Y_tr
		%~ phi_tr
		
		q_tr = sqrt(3/2) * 2 * G * norm ( strain_e_d_tr ) ;
		
		if phi_tr < 0 % elastic behavior

			stress           	= stress_tr ;
			dstressdeps       					= De ;
			acum_plas_strain 	= acum_plas_strain_n ; 
			
		else % elasto-plastic behavior phi_tr == 0

			%~ strain_e_d_tr = stress_tr \ De ; %%%% no del todo seguro
			%~ s_tr = 2*G*strain_e_d_tr ;
			
			delta_gamma = phi_tr / ( 3*G + H ) ; % 
			
			
			Nhat = strain_e_d_tr / norm(strain_e_d_tr) ;
			
			Is = eye(3) ;
			Is(3) = 1/2 ;
			
			Id = Is - 1/3*eye(3) ;
			
			stress           = ( De - delta_gamma*6*G^2 / q_tr * Id ) * strain ;
			
			dstressdeps      = 2*G*( 1 - delta_gamma*3*G/q_tr ) * Id + 6*G^2*( delta_gamma/q_tr - 1/(3*G+H) ) * kron( Nhat, Nhat') + K * eye(3) ;
			
			acum_plas_strain = acum_plas_strain_n + delta_gamma ; 

		end
		
	else	
		error("not implemented yet.")
	end
	
	% Compute Tangent Stiffness Matrix
	KTe = B' * dstressdeps * B * A * t ;
	% Internal force
	Finte = KTe * elemDisps ;

	%md### compute inertial loads and mass matrix
	Mmase = density*A/3.0 * speye(9,9);
	Fmase = Mmase * dotdotdispsElem ;

	fs = { Finte, [], Fmase } ;
	ks = { KTe,   [], Mmase } ;

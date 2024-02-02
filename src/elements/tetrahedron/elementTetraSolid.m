% Copyright 2023, ONSAS Authors (see documentation)
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
% 
% function for computation of nodal forces and tangent stiffness matrix for 3D 4 nodes 
%tetraedron element with different constitutive behaviors. 
%The hyperelastic behavior cased is based on equation 4.9.25 from belytschko 2nd edition.
function [ fs, ks, stress ] = elementTetraSolid( ...
  elemCoords, elemDisps, dotdotdispsElem, elemConstitutiveParams, density, paramOut, consMatFlag )

  % Internal flag to compute the stiffness matrix of the element using an analytic expression
  booleanKTAnalytic = 1 ;

  [ funder, ~, vol, tetCoordMat ] = computeFuncDerivVolTetraSolid( elemCoords ) ;

  % Displacements and coordinates element matrix (3dofs (ux,uy,uz) and 4 nodes)
  eleDispsMat = reshape( elemDisps, 3, 4) ;
  eleCoordSpa = tetCoordMat + eleDispsMat ;

  % Computes the gradients of the deformation function 
  H = eleDispsMat * funder' ;
  F = H + eye(3) ;

  % Compute the Green-Lagrange strain tensor
  Egreen = 0.5 * ( H + transpose( H ) + transpose( H ) * H ) ;

  if elemConstitutiveParams(1) == 2 % Saint-Venant-Kirchhoff compressible solid

    [ S, ConsMat ] = cosseratSVK( elemConstitutiveParams(2:3), Egreen, consMatFlag ) ;

  elseif elemConstitutiveParams(1) == 3 % Neo-Hookean Compressible

    [ S, ConsMat ] = cosseratNHC( elemConstitutiveParams(2:3), Egreen, consMatFlag ) ;
  end

  % Computes internal forces of the elemnt
  matBgrande = BgrandeMats ( funder , F ) ;
  Svoigt = mat2voigt( S, 1 ) ;
  Finte    = transpose(matBgrande) * Svoigt * vol ;

  % Computes element strain and stresses
  strain = zeros(6,1);
  stress = Svoigt ;

  KTe = zeros(12,12) ;
  % If paramOut == 2 the stifness matrix is returned
  if paramOut == 2

    Kml        = matBgrande' * ConsMat * matBgrande * vol ;

    matauxgeom = funder' * S * funder  * vol ;
    Kgl        = zeros(12,12) ;
    for i=1:4
      for j=1:4
        Kgl( (i-1)*3+1 , (j-1)*3+1 ) = matauxgeom(i,j);
        Kgl( (i-1)*3+2 , (j-1)*3+2 ) = matauxgeom(i,j);
        Kgl( (i-1)*3+3 , (j-1)*3+3 ) = matauxgeom(i,j);
      end
    end


    KTe = Kml + Kgl ;

  end % if param out


  Mmase = density*vol/4.0 * speye(12,12);
	Fmase = Mmase * dotdotdispsElem ;

  fs = { Finte, [], Fmase } ;
	ks = { KTe,   [], Mmase } ;


% ======================================================================
% Auxiliar functions
% ======================================================================
function matBgrande = BgrandeMats ( deriv , F )

  matBgrande = zeros(6, 12) ;

  matBgrande(1:3,:) = [ diag( deriv(:,1) )*F' diag( deriv(:,2) )*F' diag( deriv(:,3) )*F' diag( deriv(:,4) )*F' ];

  for k=1:4
      matBgrande ( 4:6 , (k-1)*3 + (1:3) ) = [ deriv(2,k)*F(:,3)'+deriv(3,k)*F(:,2)' ; ...
                                               deriv(1,k)*F(:,3)'+deriv(3,k)*F(:,1)' ; ...
                                               deriv(1,k)*F(:,2)'+deriv(2,k)*F(:,1)' ] ;
  end

function BMat = BMats ( deriv )

  BMat = zeros(6,12) ;

  for k = 1:4

    for i = 1:3
      BMat ( i , (k-1)*3 + i  ) = deriv(i,k) ;
    end

    BMat ( 4 , (k-1)*3 + 2   ) = deriv(3,k) ;
    BMat ( 4 , (k-1)*3 + 3   ) = deriv(2,k) ;

    BMat ( 5 , (k-1)*3 + 1   ) = deriv(3,k) ;
    BMat ( 5 , (k-1)*3 + 3   ) = deriv(1,k) ;

    BMat ( 6 , (k-1)*3 + 1   ) = deriv(2,k) ;
    BMat ( 6 , (k-1)*3 + 2   ) = deriv(1,k) ;

  end
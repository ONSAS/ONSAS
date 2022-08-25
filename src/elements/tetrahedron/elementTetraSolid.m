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
 
% function for computation of nodal forces and tangent stiffness matrix for 3D 4 nodes tetraedron element with different constitutive behaviors. The hyperelastic behavior cased is based on equation 4.9.25 from belytschko 2nd edition.
%

function [ Finte, KTe, stress ] = elementTetraSolid( ...
  elemCoords, elemDisps, elemConstitutiveParams, paramOut, consMatFlag )

  Finte = zeros(12,1) ;
  booleanKTAnalytic = 1 ;

  [ funder, jacobianmat, vol, tetCoordMat ] = computeFuncDerivVolTetraSolid( elemCoords ) ;

  eleDispsMat = reshape( elemDisps, 3, 4) ;
  eleCoordSpa = tetCoordMat + eleDispsMat ;

  H = eleDispsMat * funder' ;

  F = H + eye(3) ;

  Egreen = 0.5 * ( H + transpose( H ) + transpose( H ) * H ) ;

  if elemConstitutiveParams(1) == 2 % Saint-Venant-Kirchhoff compressible solid

    [ S, ConsMat ] = cosseratSVK( elemConstitutiveParams(2:3), Egreen, consMatFlag ) ;

  elseif elemConstitutiveParams(1) == 3 % Neo-Hookean Compressible

    [ S, ConsMat ] = cosseratNHC( elemConstitutiveParams(2:3), Egreen, consMatFlag ) ;
  end

  matBgrande = BgrandeMats ( funder , F ) ;

  Svoigt = mat2voigt( S, 1 ) ;

  Finte    = transpose(matBgrande) * Svoigt * vol ;

  strain = zeros(6,1);
  stress = Svoigt ;

  KTe = zeros(12,12) ;

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



  %~ mu    = E / (2.0 * ( 1.0 + nu ) ) ;
  %~ Bulk  = E / (3.0 * ( 1.0 - 2.0 * nu ) ) ;

  %~ eyetres      = eye(3)     ;
  %~ eyevoig      = zeros(6,1) ;
  %~ eyevoig(1:3) = 1.0        ;

  %~ ConsMat = zeros(6,6) ;
  %~ ConsMat (1:3,1:3) = Bulk  + 2* mu * ( eyetres - 1.0/3.0 ) ;
  %~ ConsMat (4:6,4:6) =            mu *   eyetres ;

  %~ Kml    = BMat' * ConsMat * BMat * tetVol ;

  %~ KTe([1:2:end], [1:2:end])  =  Kml ;

  %~ strain = BMat * Ue ;
  %~ stress = ConsMat * strain ;
  %~ Fint   = stress' * BMat * tetVol ;

  %~ Finte(1:2:end) = Fint ;




% ======================================================================
% ======================================================================

%~ S = lambda * trace(Egreen) * eye(3)  +  2 * shear * Egreen ;


% ======================================================================
% ======================================================================
function matBgrande = BgrandeMats ( deriv , F )

  matBgrande = zeros(6, 12) ;

  matBgrande(1:3,:) = [ diag( deriv(:,1) )*F' diag( deriv(:,2) )*F' diag( deriv(:,3) )*F' diag( deriv(:,4) )*F' ];

  for k=1:4

    %~ for i=1:3 % fila
      %~ for j=1:3 % columna
        %~ matBgrande ( i , (k-1)*3 + j  ) = deriv(i,k) * F(j,i) ;
      %~ end
    %~ end


    %~ for j=1:3
      %~ matBgrande ( 4:6 , (k-1)*3 + j ) = [ deriv(2,k)*F(j,3)+deriv(3,k)*F(j,2) ; ...
                                           %~ deriv(1,k)*F(j,3)+deriv(3,k)*F(j,1) ; ...
                                           %~ deriv(1,k)*F(j,2)+deriv(2,k)*F(j,1) ] ;
    %~ end
      matBgrande ( 4:6 , (k-1)*3 + (1:3) ) = [ deriv(2,k)*F(:,3)'+deriv(3,k)*F(:,2)' ; ...
                                               deriv(1,k)*F(:,3)'+deriv(3,k)*F(:,1)' ; ...
                                               deriv(1,k)*F(:,2)'+deriv(2,k)*F(:,1)' ] ;
  end


% ======================================================================
% ======================================================================





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

% Copyright (C) 2019, Jorge M. Perez Zerpa, J. Bruno Bazzano, Jean-Marc Battini, Joaquin Viera, Mauricio Vanzulli  
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

% function for computation of nodal forces and tangent stiffness matrix for 3D 4 nodes tetraedron element with hyperelastic behavior. based on equation 4.9.25 from belytschko 2nd edition.
%

function [ Finte, KTe, strain, stress ] = elementTetraSVKSolidInternLoadsTangMat( tetcoordmat, Ue, params, paramOut )
  
  Finte = zeros(24,1) ;
  
  booleanKTAnalytic = 1 ;
  
  elecoordmat = tetcoordmat ;
  eledispmat = reshape( Ue, 3,4) ;
  elecoordspa = tetcoordmat + eledispmat ;

  xi = 0.25 ;  wi = 1/6  ;
      
  % matriz de derivadas de fun forma respecto a coordenadas isoparametricas
  deriv = shapeFuns( xi, xi , xi , 1 )  ;

  % jacobiano que relaciona coordenadas materiales con isoparametricas
  jacobianmat = elecoordmat * deriv'  ;

  vol = det( jacobianmat ) / 6 ;
  if vol<0,  vol, error('Element with negative volume, check connectivity.'), end

  funder = inv(jacobianmat)' * deriv ;
 
  H = eledispmat * funder' ;

  F = H + eye(3) ;

  %~ Egreen = 0.5 * ( H + H' + H' * H ) ;
  Egreen = 0.5 * ( H + transpose(H) + transpose(H) * H ) ;

  % constitutive params
  young  = params(1) ;
  nu     = params(2) ;

  S = cosserat( young, nu, Egreen) ;

  matBgrande = BgrandeMats ( funder , F ) ;
     
  Svoigt = tranvoigSin2( S ) ;
  
  ConsMat = constTensor ( young, nu, Egreen ) ;
  
  Scons = ConsMat * tranvoigCon2(Egreen) ;

  %~ Finte(1:2:end)    = matBgrande' * Svoigt * vol ;
  Finte(1:2:end)    = transpose(matBgrande) * Svoigt * vol ;
    
  strain = zeros(6,1);
  stress = zeros(6,1);

  if paramOut == 2

    KTe = sparse(24,24) ;

    if booleanKTAnalytic
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
      
      KTe(1:2:end,1:2:end) = Kml + Kgl ;

    else
    
      step = 1e-8;
      
      for i=1:2:12
        ei = zeros(12,1);   ei(i) = step * j ;
        
        FinteComp = elementTetraSVKSolidInternLoadsTangMat( tetcoordmat, Ue + ei, params, 1 ) ; 
        KTe(:,i) = imag( FinteComp ) / step;
      end % for 
      
    end % if analytic calc
  end % if param out




% ======================================================================
% ======================================================================
function S = cosserat( young, nu, Egreen)

lambda  = young * nu / ( (1 + nu) * (1 - 2*nu) ) ;
shear   = young      / ( 2 * (1 + nu) )          ;

S = lambda * trace(Egreen) * eye(3)  +  2 * shear * Egreen ;


% ======================================================================
% ======================================================================
function matBgrande = BgrandeMats ( deriv , F )

  matBgrande = zeros(6, 12) ;
  
  for k=1:4

    for i=1:3 % fila
      for j=1:3 % columna
        matBgrande ( i , (k-1)*3 + j  ) = deriv(i,k) * F(j,i) ;
      end
    end          

    for j=1:3
      matBgrande ( 4 , (k-1)*3 + j ) = deriv(2,k) * F(j,3) + deriv(3,k) * F(j,2) ;
      matBgrande ( 5 , (k-1)*3 + j ) = deriv(1,k) * F(j,3) + deriv(3,k) * F(j,1) ;
      matBgrande ( 6 , (k-1)*3 + j ) = deriv(1,k) * F(j,2) + deriv(2,k) * F(j,1) ;
    end
  end


% ======================================================================
% ======================================================================
function v = tranvoigCon2(Tensor)
    
  if norm( Tensor - Tensor' ) > 1e-10
    Tensor
    norm( Tensor - Tensor' )
    %~ error('tensor not symmetric')
  end
  
  v = zeros(6,1) ;
  
  v(1) = Tensor(1,1) ;
  v(2) = Tensor(2,2) ;
  v(3) = Tensor(3,3) ;
  v(4) = Tensor(2,3)*2 ;
  v(5) = Tensor(1,3)*2 ;
  v(6) = Tensor(1,2)*2 ;

function v = tranvoigSin2(Tensor)
    
  if norm( Tensor - Tensor' ) > 1e-5
    %~ Tensor
    norm( Tensor - Tensor' )
    %~ error('tensor not symmetric')
  end
  
  v = zeros(6,1) ;
  
  v(1) = Tensor(1,1) ;
  v(2) = Tensor(2,2) ;
  v(3) = Tensor(3,3) ;
  v(4) = Tensor(2,3) ;
  v(5) = Tensor(1,3) ;
  v(6) = Tensor(1,2) ;


% ======================================================================
% ======================================================================
function ConsMat = constTensor ( young, nu, Egreen )

  booleanCAnalytic = 1 ;

  ConsMat= zeros(6,6);
  
  if booleanCAnalytic 

    shear   = young / ( 2 * (1+ nu) ) ;

    ConsMat (1,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ 1-nu , nu   , nu   ] ; 
    ConsMat (2,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ nu   , 1-nu , nu   ] ;
    ConsMat (3,1:3) = ( shear / (1 - 2 * nu) ) * 2 * [ nu   , nu   , 1-nu ] ;
    ConsMat (4,4  ) = shear ;
    ConsMat (5,5  ) = shear ;
    ConsMat (6,6  ) = shear ;

  else
  
    compStep = 1e-8 ;
  
    for i=1:6
  
      direction = zeros( 3,3) ;
      if i<4
        direction(i,i) = 1 ;
      elseif i==4
        direction(2,3) = 1 ;
        direction(3,2) = 1 ;
      elseif i==5
        direction(1,3) = 1 ;
        direction(3,1) = 1 ;
      elseif i==6
        direction(1,2) = 1 ;
        direction(2,1) = 1 ;
      end
  
      EgreenComp = Egreen  + compStep * direction * j ;
    
      Scomp = cosserat( young, nu, EgreenComp) ;
      
      dsde = tranvoigSin2 ( imag( Scomp ) / compStep ) ;
    
      %~ ConsMat(1:3,i) = dsde(1:3) ;
      %~ ConsMat(4:6,i) = dsde(4:6)*0.5 ;
      ConsMat(:,i) = dsde ;

    end % for compone
  end


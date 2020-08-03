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

% function for computation of the normal force and tanget matrix 
% of 3D truss elements using engineering strain.




% Co-rotational Truss with Engineering strain




function [Finte, KTe, stress, dstressdeps, strain ] = ...
  elementTrussInternForce( Xe, Ue, hyperelasparams, A , paramout, booleanCSTangs )

  % vector reference element
  localAxisRef = Xe( [(1:2:5)+6]) - Xe( [(1:2:5)]) ;
  
  % initial length
  lini = sqrt( sum( localAxisRef.^2 ) ) ;

  % normalized reference vector
  e1ref = localAxisRef / lini ;
  
  % vector deformed element
  Xedef = Xe + Ue ;
  
  localAxisDef = Xedef( [(1:2:5)+6]) - Xedef( [(1:2:5)]) ;

  % deformed length
  ldef = sqrt( sum( localAxisDef.^2 ) ) ;

  % normalized deformed co-rotational vector
  e1def = localAxisDef / ldef ;

  % --- strain ---
  strain = ( ldef^2 - lini^2 ) / ( lini * (lini + ldef) ) ;

  % --- stress and constitutive tensor ---
  [ stress, dstressdeps  ] = hyperElasModels ( strain, hyperelasparams ) ;
  
  TTcl              = zeros(12,1) ;
  TTcl(1:2:5)       = -e1def ;    TTcl([(1:2:5)+6]) =  e1def ;

  % internal forces computation
  Finte = stress * A * TTcl  ;

  if paramout == 2
  
    KTe = zeros(12,12) ;

    % tangent matrix computation  
    if booleanCSTangs == 0
    
      B = zeros( 12,3) ;
      B(1:2:5     , :) = -eye(3);   B([1:2:5]+6 , :) = +eye(3);
    
      KMe   = dstressdeps * A / lini * ( TTcl * TTcl' ) ;
      Ksige = A    * stress / (ldef) * ( B * B' - TTcl * (TTcl')  ) ;
    
      KTe = KMe + Ksige ;
    
    else
      %
      step = 1e-10;
      
      for i=1:2:12
        ei = zeros(12,1);   ei(i) = j ;
        
        FinteComp = elementTrussInternForce( Xe, Ue + ei*step , hyperelasparams, A, 1 );
        
        KTe(i,1:2:end) = imag( FinteComp ) / step;
      end
    end

  else
    KTe = [] ;
  end
  
  % reduce to keep only displacements dofs
  Finte = Finte(1:2:end);
  KTe   = KTe(1:2:end,1:2:end);


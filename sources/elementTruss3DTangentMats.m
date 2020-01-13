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

%

function [ KTe, KL0e ] = elementTruss3DTangentMats( ...
  Xe, Ue, hyperelasparams, A )

  Xe = (Xe(:))';
  Ue = (Ue(:))';
  
  booleanAnalytical = 1 ;
  
  if booleanAnalytical
    Xe ;
    localAxisRef = ( Xe( [(1:2:5)+6]) - Xe( [(1:2:5)]) )' ;
    lini = sqrt( sum( localAxisRef.^2 ) ) ;
  
    e1ref = localAxisRef / lini ;
    
    Xedef = Xe + Ue ;
    localAxisDef = ( Xedef( [(1:2:5)+6] ) - Xedef( [(1:2:5)]) )' ;
    ldef = sqrt( sum( localAxisDef.^2 ) ) ;
    
    e1def = localAxisDef / ldef ;
  
    strain = ( ldef^2 - lini^2 ) / ( lini * (lini + ldef) ) ;
  
    [ stress, dstressdeps  ] = hyperElasModels ( strain, hyperelasparams ) ;
  
    B = zeros( 12,3) ;
    B(1:2:5     , :) = -eye(3);
    B([1:2:5]+6 , :) = +eye(3);
  
    TTcl = B * e1def;
  
    KMe = dstressdeps * A / lini * ( TTcl * TTcl' ) ;
    
    Ksige = A * stress / (ldef) * ( B * B' - TTcl * (TTcl')  ) ;
  
    KTe = KMe + Ksige ;
  
    TrefTcl = B * e1ref;
    [ ~, dstressdeps0  ] = hyperElasModels ( +eps, hyperelasparams ) ;
    KL0e = dstressdeps0 * A / lini * ( TrefTcl * TrefTcl' ) ;





  else
    %
    KTe = zeros(12,12) ;

    %~ FinteU = elementTruss3DInternLoads( Xe, Ue, hyperelasparams, A );
    step = 1e-10;
    
    for i=1:2:12
      ei = zeros(12,1);
      ei(i) = j ;
      
      FinteComp = elementTruss3DInternLoads( Xe, Ue + ei*step , hyperelasparams, A );
      
      KTe(i,:) = imag( FinteComp ) / step;
      %~ pause
    end
    
    KL0e = [];
    %~ KTe
    %~ stop
  end

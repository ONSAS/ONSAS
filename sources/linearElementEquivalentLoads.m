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

% --------------------------------------------------------------------------------------------------

function unifDisNodal = linearElementEquivalentLoads( unifDisLoadG, ElemLengths, Conec,ndofpnode, Local2GlobalMats, 
  
  for i = 1:size(unifDisLoadG,1)

    nelem = unifDisLoadG(i,1) ;

    m = indexesElems( nelem ) ;

    if Conec(nelem,7) == 2
      l = ElemLengths(m) ;
      nodeselem = Conec(nelem,1:2)' ;
      aux = nodes2dofs( nodeselem, ndofpnode ) ;

      R   = RotationMatrix ( ndofpnode, Local2GlobalMats{m} ) ;
      exL = Local2GlobalMats{m}(:,1) ;
      eyL = Local2GlobalMats{m}(:,2) ;
      ezL = Local2GlobalMats{m}(:,3) ;

      if unifDisLoadG(i,2)
        qUnifDisxG = unifDisLoadG(i,5) ;
        ezL = Local2GlobalMats{m}(:,3) ;
        cosAlpha = ezL' * [1 0 0]' / norm(ezL) ; alpha = acos(cosAlpha) ; 
        senAlpha = cos(pi/2-alpha) ;
        q_perp = qUnifDisxG * cosAlpha * l * [ 0 0 0 -l/12 1/2 0 0 0 0 l/12 1/2 0 ]' ;
        q_long = qUnifDisxG * senAlpha * l * [ 1/2  0 0   0    0  0 1/2 0 0   0   0  0 ]' ;
        elemUnifDis{nelem,1} = R * q_perp + R * q_long ;
        unifDisNodal( aux ) = elemUnifDis{nelem,1} + unifDisNodal( aux ) ;
      end 
       
      if unifDisLoadG(i,3)
        qUnifDisyG = unifDisLoadG(i,6) ;
        eyL = Local2GlobalMats{m}(:,2) ;
        cosAlpha = eyL' * [0 1 0]' / norm(eyL) ; alpha = acos(cosAlpha) ; 
        senAlpha = cos(pi/2-alpha) ;
        q_perp = qUnifDisyG * cosAlpha * l * [ 0 0 1/2 0 0 l/12 0 0 1/2 0 0 -l/12 ]' ;
        q_long = qUnifDisyG * senAlpha * l * [ 1/2  0 0   0    0  0 1/2 0 0   0   0  0 ]' ;
        elemUnifDis{nelem,2} = R * q_perp + R * q_long ;
        unifDisNodal( aux ) = elemUnifDis{nelem,2} + unifDisNodal( aux ) ;
      end
      
      if unifDisLoadG(i,4) 
        qUnifDiszG = unifDisLoadG(i,7) ;
        ezL = Local2GlobalMats{m}(:,3) ;
        cosAlpha = ezL' * [0 0 1]' / norm(ezL) ; alpha = acos(cosAlpha) ; 
        senAlpha = cos(pi/2-alpha) ;
        q_perp = qUnifDiszG * cosAlpha * l * [ 0 0 0 -l/12 1/2 0 0 0 0 l/12 1/2 0 ]' ;
        q_long = qUnifDiszG * senAlpha * l * [ 1/2  0 0   0    0  0 1/2 0 0   0   0  0 ]' ;
        elemUnifDis{nelem,3} = R * q_perp + R * q_long ;
        unifDisNodal( aux ) = elemUnifDis{nelem,3} + unifDisNodal( aux ) ;  
      end 
   
    end
  end

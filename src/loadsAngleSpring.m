% Copyright (C) 2021, Jorge M. Perez Zerpa, J. Bruno Bazzano, Joaquin Viera,
%   Mauricio Vanzulli, Marcelo Forets, Jean-Marc Battini, Sebastian Toro  
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

% function for testing experimental feature.

function [ fextAngSpr, KTAngSprAppr] = loadsAngleSpring( Nodes, Conec, bendStiff )

n = size(Nodes,1);

fextAngSpr   = zeros ( 6*n,   1);
KTAngSprAppr = sparse( 6*n, 6*n);

for i = 1:n
  %

  if bendStiff(i) > 0
    %
    if Conec(i-1  ,2) ~= Conec(i,1);
      error('conectivity error for angle springs.');
    end

    EI = bendStiff (i);

    indxim1 = Conec(i-1,1) ;
    indxi   = Conec(i-1,2) ;
    indxip1 = Conec(i  ,2) ;

    ti   = ( Nodes(indxi  ,:) - Nodes(indxim1,:) )' ;
    tip1 = ( Nodes(indxip1,:) - Nodes(indxi  ,:) )' ;

    li   = norm( ti  ) ;
    lip1 = norm( tip1) ;

    dtds    = ( tip1/lip1 - ti/li ) / ( lip1/2 + li/2 ) ;

    if norm( dtds ) > 0
      Fiim1 = +dtds *  EI  / norm( cross( ti  , dtds) )  ;
      Fim1i = -dtds *  EI  / norm( cross( ti  , dtds) )  ;

      Fiip1 = +dtds *  EI  / norm( cross( tip1, dtds) )  ;
      Fip1i = -dtds *  EI  / norm( cross( tip1, dtds) )  ;


      aux3 = nodes2dofs( [indxim1 indxi indxip1] ,6) ;
      Avec = -EI /  ( lip1/2 + li/2 )^3 * [ 1 -2 1 ];

      % node im1
      aux = nodes2dofs( indxim1,6) ;
      fextAngSpr( aux(1:2:6) ) += Fim1i ;

      for indcoord=1:2:5
        KTAngSprAppr( aux(indcoord) , aux3(indcoord:6:end) ) +=  -Avec ;
      end
      % ---------------


      % node i
      aux = nodes2dofs( indxi  ,6) ;
      fextAngSpr( aux(1:2:6) ) += Fiim1 + Fiip1 ;

      for indcoord=1:2:5
        KTAngSprAppr( aux(indcoord) , aux3(indcoord:6:end) ) +=  2*Avec ;
      end
      % ---------------


      % node ip1
      aux = nodes2dofs( indxip1,6) ;
      fextAngSpr( aux(1:2:6) ) += Fip1i ;

      for indcoord=1:2:5
        KTAngSprAppr( aux(indcoord) , aux3(indcoord:6:end) ) +=  -Avec ;
      end
      % ---------------

    end
  end
end

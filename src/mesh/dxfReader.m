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

function [ nodesMat, conecMat, physicalNames ] = dxfReader( fileName )

  [c_Line, c_Poly, c_Cir, c_Arc, c_Poi] = f_LectDxf( fileName ) ;

  ndofpnode = 6 ;
  nelems = size(c_Line,1) ;

  Nodes = [] ;
  conecMat = cell(nelems,1) ;
  pos = 1 ;

  nnodesEntity = size(c_Poi,1) ;
	auxLength = length(c_Poi{1}) ;
  typesElem = [] ;

  for i = 1:nelems
    aux1 = [] ;
    aux2 = [] ;
    nod1 = c_Line{i}(1:3) ;
    nod2 = c_Line{i}(4:6) ;
    nnodes = size(Nodes,1) ;

    % MELCS
    % M: Material
    % E: Element type
    % L: Load
    % C: Cross section
    % S: Springs

		matEntity = str2num(c_Line{i,2}(1:2)) ;
    elemEntity = str2num(c_Line{i,2}(4:5)) ;
    loadEntity = str2num(c_Line{i,2}(7:8)) ;
    secEntity = str2num(c_Line{i,2}(10:11)) ;
    suppEntity = str2num(c_Line{i,2}(13:14)) ;

    entityVec = [ elemEntity loadEntity secEntity suppEntity ] ;

    typesElem = [ typesElem ; elemEntity ] ;

    if i == 1
			Nodes = [ Nodes ; nod1 ; nod2 ] ;
      conecMat{i,1} = [ matEntity entityVec pos pos+1 ] ;
      pos = pos+2 ;
    else
      for j = 1:nnodes
        if norm(Nodes(j,:)-nod1) < 1e-10
          aux1 = [ j ] ;
        elseif norm(Nodes(j,:)-nod2) < 1e-10
          aux2 = [ j ] ;
        end
      end


      if length(aux1) == 0 && length(aux2) == 0
          Nodes = [ Nodes ; nod1 ; nod2 ] ;
          conecMat{i,1} = [ matEntity entityVec pos pos+1 ] ;
          pos = pos+2 ;
      elseif length(aux1) == 0 && length(aux2) == 1
          Nodes = [ Nodes ; nod1 ] ;
          conecMat{i,1} = [ matEntity entityVec pos aux2 ] ;
          pos = pos+1 ;
      elseif length(aux1) == 1 &&  length(aux2) == 0
          Nodes = [ Nodes ; nod2 ] ;
          conecMat{i,1} = [ matEntity entityVec aux1 pos ] ;
          pos = pos+1 ;
      elseif length(aux1) == 1 && length(aux2) == 1

          conecMat{i,1} = [ matEntity entityVec aux1 aux2 ] ;

      end

    end

  end

  typesElem = unique(typesElem);

	nodesMat = Nodes ;

	lastElemType = typesElem(end) ;

  nnodes = size(Nodes,1) ;
  nnodesEntity = size(c_Poi,1) ;
  auxLength = length(c_Poi{1}) ;
  if auxLength > 0
    for i = 1:nnodesEntity
      nod = c_Poi{i}(1:3) ;
      aux = [] ;
      elemEntity = str2num(c_Poi{i,2}(4:5)) ;
      for j = 1:nnodes
        if norm(Nodes(j,:)-nod) < 1e-10
          nod = j ;

          suppEntity = str2num(c_Poi{i,2}(13:14)) ;
					loadEntity = str2num(c_Poi{i,2}(7:8)) ;
					entityVec = [ 0 elemEntity loadEntity 0 suppEntity  ] ;

          conecMat{nelems+1,1} = [ entityVec nod ] ;
          nelems=nelems+1 ;

        end
      end
    end

  end % for elements

end % function
